//! The truth manifest: what a study was asked to show and what it actually shows.
//!
//! Targets come from the resolved controls; realised values are computed from the generated
//! ADaM datasets (never from simulation state), exactly as an analyst would compute them.

use std::collections::BTreeMap;

use serde_json::{Value, json};
use trialforge_core::{Dataset, Date, Values};
use trialforge_designs::{Course, Resolved};

fn chars<'a>(ds: &'a Dataset, name: &str) -> Option<&'a [Option<String>]> {
    ds.column(name)?.as_char()
}

fn nums(ds: &Dataset, name: &str) -> Option<Vec<Option<f64>>> {
    match &ds.column(name)?.values {
        Values::Num(v) => Some(v.clone()),
        Values::Int(v) => Some(v.iter().map(|x| x.map(|n| n as f64)).collect()),
        _ => None,
    }
}

fn dates(ds: &Dataset, name: &str) -> Option<Vec<Option<Date>>> {
    match &ds.column(name)?.values {
        Values::Date(v) => Some(v.clone()),
        _ => None,
    }
}

fn round4(x: f64) -> f64 {
    trialforge_core::math::round_dp(x, 4)
}

fn mean(xs: &[f64]) -> Option<f64> {
    (!xs.is_empty()).then(|| xs.iter().sum::<f64>() / xs.len() as f64)
}

/// Kaplan–Meier median from (time, event) pairs.
pub fn km_median(mut data: Vec<(f64, bool)>) -> Option<f64> {
    data.sort_by(|a, b| a.0.total_cmp(&b.0).then(b.1.cmp(&a.1)));
    let mut at_risk = data.len() as f64;
    let mut surv = 1.0;
    let mut i = 0;
    while i < data.len() {
        let t = data[i].0;
        let (mut events, mut total) = (0.0, 0.0);
        while i < data.len() && data[i].0 == t {
            events += f64::from(u8::from(data[i].1));
            total += 1.0;
            i += 1;
        }
        if events > 0.0 {
            surv *= 1.0 - events / at_risk;
            if surv <= 0.5 {
                return Some(t);
            }
        }
        at_risk -= total;
    }
    None
}

/// Builds the manifest.
pub(crate) fn manifest(r: &Resolved, archetype: &str, sdtm: &[Dataset], adam: &[Dataset]) -> Value {
    let find = |name: &str| adam.iter().chain(sdtm).find(|d| d.name == name);
    let arms: Vec<String> = r.design.arms.iter().map(|a| a.arm.clone()).collect();
    let mut out = json!({
        "archetype": archetype,
        "studyid": r.design.studyid,
        "seed": r.design.seed,
        "n": r.design.n,
        "controls": r.controls,
        "arms": arms,
    });
    let adsl = find("ADSL");
    let by_arm = |values: Vec<(String, f64)>| -> BTreeMap<String, Vec<f64>> {
        let mut m: BTreeMap<String, Vec<f64>> = BTreeMap::new();
        for (a, v) in values {
            m.entry(a).or_default().push(v);
        }
        m
    };

    match &r.course {
        Course::Fixed(p) => {
            let endpoint = &p.endpoints[0];
            let eot = r
                .design
                .schedule
                .as_ref()
                .and_then(|s| {
                    s.visits
                        .iter()
                        .find(|v| v.kind == trialforge_model::VisitKind::EndOfTreatment)
                })
                .map(|v| v.visitnum as f64);
            let bds = ["ADLB", "ADVS"].iter().filter_map(|n| find(n)).find(|d| {
                chars(d, "PARAMCD")
                    .is_some_and(|c| c.iter().any(|x| x.as_deref() == Some(&endpoint.testcd)))
            });
            let mut realised = BTreeMap::new();
            if let (Some(d), Some(eot)) = (bds, eot) {
                let (pc, trtp, anl, avn) = (
                    chars(d, "PARAMCD"),
                    chars(d, "TRTP"),
                    chars(d, "ANL01FL"),
                    nums(d, "AVISITN"),
                );
                let chg = nums(d, "CHG");
                if let (Some(pc), Some(trtp), Some(anl), Some(avn), Some(chg)) =
                    (pc, trtp, anl, avn, chg)
                {
                    let vals = (0..d.n_rows())
                        .filter(|&i| {
                            pc[i].as_deref() == Some(&endpoint.testcd)
                                && anl[i].as_deref() == Some("Y")
                                && avn[i] == Some(eot)
                                && chg[i].is_some()
                        })
                        .filter_map(|i| Some((trtp[i].clone()?, chg[i]?)))
                        .collect();
                    for (arm, xs) in by_arm(vals) {
                        realised.insert(
                            arm,
                            json!({"mean_change": mean(&xs).map(round4), "n": xs.len()}),
                        );
                    }
                }
            }
            let targets: BTreeMap<&String, f64> = arms
                .iter()
                .zip(&endpoint.change_by_arm)
                .map(|(a, c)| (a, *c))
                .collect();
            let diff_target = endpoint
                .change_by_arm
                .last()
                .zip(endpoint.change_by_arm.first())
                .map(|(a, b)| a - b);
            let diff_realised = (|| {
                let a = realised.get(arms.last()?)?.get("mean_change")?.as_f64()?;
                let b = realised.get(arms.first()?)?.get("mean_change")?.as_f64()?;
                Some(round4(a - b))
            })();
            out["endpoint"] = json!({
                "paramcd": endpoint.testcd,
                "visit": "end of treatment",
                "target_mean_change": targets,
                "realised": realised,
                "target_difference": diff_target,
                "realised_difference": diff_realised,
            });

            if let Some(sl) = adsl {
                let (trt, dcs) = (chars(sl, "TRT01A"), chars(sl, "DCSREAS"));
                if let (Some(trt), Some(dcs)) = (trt, dcs) {
                    let weeks =
                        f64::from(r.design.schedule.as_ref().map_or(0, |s| s.treatment_days)) / 7.0;
                    let expected =
                        1.0 - trialforge_core::math::exp(-p.dropout_annual_rate * weeks / 52.18);
                    let mut per = BTreeMap::new();
                    for (arm, xs) in by_arm(
                        (0..sl.n_rows())
                            .filter_map(|i| {
                                Some((trt[i].clone()?, f64::from(u8::from(dcs[i].is_some()))))
                            })
                            .collect(),
                    ) {
                        per.insert(
                            arm,
                            json!({"proportion": mean(&xs).map(round4), "n": xs.len()}),
                        );
                    }
                    out["discontinuation"] = json!({
                        "target_dropout_excluding_ae_and_death": round4(expected),
                        "realised_any_reason": per,
                    });
                }
            }
            if let (Some(ae), Some(sl)) = (find("ADAE"), adsl) {
                let mut teae: BTreeMap<String, std::collections::BTreeSet<String>> =
                    BTreeMap::new();
                if let (Some(u), Some(t), Some(fl)) = (
                    chars(ae, "USUBJID"),
                    chars(ae, "TRTA"),
                    chars(ae, "TRTEMFL"),
                ) {
                    for i in 0..ae.n_rows() {
                        if let (Some(u), Some(t), Some("Y")) = (&u[i], &t[i], fl[i].as_deref()) {
                            teae.entry(t.clone()).or_default().insert(u.clone());
                        }
                    }
                }
                let denom = by_arm(chars(sl, "TRT01A").map_or(vec![], |t| {
                    t.iter().flatten().map(|a| (a.clone(), 1.0)).collect()
                }));
                let mut inc = BTreeMap::new();
                for (arm, d) in denom {
                    let k = teae.get(&arm).map_or(0, |s| s.len());
                    inc.insert(arm, json!({"subjects_with_teae": k, "n": d.len(), "proportion": round4(k as f64 / d.len() as f64)}));
                }
                out["teae_incidence"] = json!(inc);
            }
        }
        Course::Oncology(p) => {
            out["latent_median_pd_months"] = json!(
                arms.iter()
                    .zip(&p.median_pd_months)
                    .map(|(a, m)| (a.clone(), *m))
                    .collect::<BTreeMap<_, _>>()
            );
            if let Some(sl) = adsl {
                let (trt, rand, dth, lst) = (
                    chars(sl, "TRT01P"),
                    dates(sl, "RANDDT"),
                    dates(sl, "DTHDT"),
                    dates(sl, "LSTALVDT"),
                );
                if let (Some(trt), Some(rand), Some(dth), Some(lst)) = (trt, rand, dth, lst) {
                    let mut per: BTreeMap<String, Vec<(f64, bool)>> = BTreeMap::new();
                    for i in 0..sl.n_rows() {
                        let (Some(a), Some(r0)) = (&trt[i], rand[i]) else {
                            continue;
                        };
                        let end = dth[i].or(lst[i]);
                        let Some(end) = end else { continue };
                        let months = f64::from(end.days_since(r0) + 1)
                            / trialforge_processes::oncology::DAYS_PER_MONTH;
                        per.entry(a.clone())
                            .or_default()
                            .push((months, dth[i].is_some()));
                    }
                    let os: BTreeMap<String, Value> = per
                        .into_iter()
                        .map(|(a, d)| {
                            let deaths = d.iter().filter(|x| x.1).count();
                            (a, json!({"n": d.len(), "deaths": deaths, "km_median_os_months": km_median(d).map(round4)}))
                        })
                        .collect();
                    out["overall_survival"] = json!(os);
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::km_median;

    #[test]
    fn km_median_simple() {
        // Survival steps to 0.75, 0.5 at t = 2 → median 2.
        assert_eq!(
            km_median(vec![(1.0, true), (2.0, true), (3.0, false), (4.0, true)]),
            Some(2.0)
        );
        assert_eq!(km_median(vec![(1.0, false), (2.0, false)]), None);
    }
}
