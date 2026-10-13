//! The truth manifest: what a study was asked to show and what it actually shows.
//!
//! Targets come from the resolved controls; realised values are computed from the generated
//! ADaM datasets (never from simulation state), exactly as an analyst would compute them.

use std::collections::BTreeMap;

use serde_json::{Value, json};
use trialforge_core::stats::{Survival, km_median};
use trialforge_core::{Dataset, Values};
use trialforge_designs::{Course, Resolved};
use trialforge_processes::oncology::DAYS_PER_MONTH;

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

fn round4(x: f64) -> f64 {
    trialforge_core::math::round_dp(x, 4)
}

fn mean(xs: &[f64]) -> Option<f64> {
    (!xs.is_empty()).then(|| xs.iter().sum::<f64>() / xs.len() as f64)
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
        Course::Oncology(_, cal) => {
            let c = &r.controls;
            out["targets"] = json!({
                "os_median_months_control": c["os_median_months"],
                "os_hr": c["os_hr"],
                "pfs_median_months_control": c["pfs_median_months"],
                "pfs_hr": c["pfs_hr"],
                "confirmed_orr": arms.iter().zip(c["orr"].as_array().into_iter().flatten()).map(|(a, v)| (a.clone(), v.clone())).collect::<BTreeMap<_, _>>(),
            });
            out["calibration"] = json!({
                "confirmed_orr_on_calibration_cohort": arms.iter().zip(&cal.achieved_orr).map(|(a, v)| (a.clone(), round4(*v))).collect::<BTreeMap<_, _>>(),
                "orr_targets_reachable": cal.orr_reachable,
            });
            if let Some(tte) = find("ADTTE") {
                for (param, key) in [
                    ("OS", "overall_survival"),
                    ("PFS", "progression_free_survival"),
                ] {
                    if let Some(v) = tte_summary(tte, param, &arms) {
                        out[key] = v;
                    }
                }
            }
            if let Some(rs) = find("ADRS") {
                out["confirmed_orr"] = orr_summary(rs, &arms);
            }
        }
    }
    out
}

/// KM medians by arm and the stratified Cox HR (first arm vs the last) for one ADTTE parameter.
fn tte_summary(tte: &Dataset, param: &str, arms: &[String]) -> Option<Value> {
    let (pc, trt, strat) = (
        chars(tte, "PARAMCD")?,
        chars(tte, "TRT01P")?,
        chars(tte, "STRAT1R"),
    );
    let (aval, cnsr) = (nums(tte, "AVAL")?, nums(tte, "CNSR")?);
    let rows: Vec<usize> = (0..tte.n_rows())
        .filter(|&i| pc[i].as_deref() == Some(param))
        .collect();
    let obs = |i: usize| Survival {
        time: aval[i].unwrap_or(0.0) / DAYS_PER_MONTH,
        event: cnsr[i] == Some(0.0),
    };
    let mut per_arm = BTreeMap::new();
    for arm in arms {
        let d: Vec<Survival> = rows
            .iter()
            .copied()
            .filter(|&i| trt[i].as_deref() == Some(arm))
            .map(obs)
            .collect();
        let events = d.iter().filter(|x| x.event).count();
        per_arm.insert(
            arm.clone(),
            json!({"n": d.len(), "events": events, "km_median_months": km_median(&d).map(round4)}),
        );
    }
    let mut v = json!({"by_arm": per_arm});
    if arms.len() == 2 {
        let levels: Vec<String> = strat.map_or(vec![], |s| {
            let mut l: Vec<String> = rows.iter().filter_map(|&i| s[i].clone()).collect();
            l.sort();
            l.dedup();
            l
        });
        let data: Vec<Survival> = rows.iter().map(|&i| obs(i)).collect();
        let x: Vec<f64> = rows
            .iter()
            .map(|&i| f64::from(u8::from(trt[i].as_deref() == Some(&arms[0]))))
            .collect();
        let st: Vec<usize> = rows
            .iter()
            .map(|&i| {
                strat
                    .and_then(|s| s[i].as_ref())
                    .and_then(|v| levels.iter().position(|l| l == v))
                    .unwrap_or(0)
            })
            .collect();
        if let Some(fit) = trialforge_core::stats::cox_binary(&data, &x, &st) {
            let (lo, hi) = fit.hr_ci(0.95);
            v["stratified_cox_hr"] = json!({
                "comparison": format!("{} vs {}", arms[0], arms[1]),
                "hr": round4(fit.hr()),
                "ci95": [round4(lo), round4(hi)],
            });
        }
    }
    Some(v)
}

/// Confirmed ORR by arm from ADRS (`CBOR` in CR/PR), over the ITT population.
fn orr_summary(rs: &Dataset, arms: &[String]) -> Value {
    let (Some(pc), Some(trt), Some(av)) = (
        chars(rs, "PARAMCD"),
        chars(rs, "TRT01P"),
        chars(rs, "AVALC"),
    ) else {
        return Value::Null;
    };
    let mut out = BTreeMap::new();
    for arm in arms {
        let rows: Vec<usize> = (0..rs.n_rows())
            .filter(|&i| pc[i].as_deref() == Some("CBOR") && trt[i].as_deref() == Some(arm))
            .collect();
        let resp = rows
            .iter()
            .filter(|&&i| matches!(av[i].as_deref(), Some("CR" | "PR")))
            .count();
        out.insert(arm.clone(), json!({"n": rows.len(), "responders": resp, "proportion": round4(resp as f64 / rows.len().max(1) as f64)}));
    }
    json!(out)
}
