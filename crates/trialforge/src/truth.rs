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
        Course::Pk(_, pk) => {
            out["pk"] = pk_truth(r, pk, adam);
        }
        Course::Adas(_, adas) => {
            out["adas_cog"] = adas_truth(r, adas, adam);
        }
        Course::Outcomes(p) => {
            out["outcomes"] = outcomes_truth(r, p, adam);
        }
        Course::Crossover(_, pk, _) => {
            out["bioequivalence"] = be_truth(r, pk, adam);
        }
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
            let mut groups: BTreeMap<String, Vec<f64>> = BTreeMap::new();
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
                    groups = by_arm(vals);
                    for (arm, xs) in groups.clone() {
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
            if archetype == "dose_ranging_ph2" {
                out["dose_response"] = dose_response(r, &groups);
            }
            if let Some(margin) = r.controls.get("ni_margin").and_then(Value::as_f64) {
                out["non_inferiority"] = non_inferiority(&arms, &groups, margin);
            }
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

/// Phase 1 PK: the population model and what NCA (ADPP) recovers: geometric means by dose
/// and the power-model slope of exposure on dose (90% CI) against the proportionality target.
fn pk_truth(r: &Resolved, pk: &trialforge_processes::pk::PkParams, adam: &[Dataset]) -> Value {
    use trialforge_core::math;
    let mut out = json!({
        "analyte": pk.analyte,
        "model": {
            "compartments": pk.compartments,
            "cl_l_per_h": pk.cl,
            "v_l": pk.v,
            "ka_per_h": pk.ka,
            "half_life_h": round4(std::f64::consts::LN_2 * pk.v / pk.cl),
            "proportionality": pk.proportionality,
        },
    });
    let Some(adpp) = adam.iter().find(|d| d.name == "ADPP") else {
        return out;
    };
    let (Some(pc), Some(trt), Some(aval), Some(avisitn)) = (
        chars(adpp, "PARAMCD"),
        chars(adpp, "TRT01P"),
        nums(adpp, "AVAL"),
        nums(adpp, "AVISITN"),
    ) else {
        return out;
    };
    let dose_of: BTreeMap<&str, f64> = r
        .design
        .arms
        .iter()
        .map(|a| (a.arm.as_str(), a.treatment.dose))
        .collect();
    let multiple = pk.planned_doses > 1;
    // The profile the summary uses: the single dose, or the last of several.
    let last_day = |i: usize| {
        !multiple
            || avisitn[i].is_some_and(|d| {
                (0..adpp.n_rows())
                    .filter(|&j| trt[j] == trt[i])
                    .filter_map(|j| avisitn[j])
                    .all(|e| e <= d)
            })
    };
    let params: &[&str] = if multiple {
        &["CMAX", "AUCTAU", "CLFTAU", "LAMZHL", "ARAUC"]
    } else {
        &["CMAX", "AUCIFO", "CLFO", "LAMZHL"]
    };
    let mut values: BTreeMap<(&str, String), Vec<f64>> = BTreeMap::new();
    for i in 0..adpp.n_rows() {
        let (Some(p), Some(t), Some(v)) = (pc[i].as_deref(), trt[i].as_deref(), aval[i]) else {
            continue;
        };
        if params.contains(&p) && v > 0.0 && (p == "ARAUC" || last_day(i)) {
            values.entry((p, t.to_string())).or_default().push(v);
        }
    }
    let gm = |xs: &[f64]| math::exp(xs.iter().map(|x| math::ln(*x)).sum::<f64>() / xs.len() as f64);
    let mut by_dose = BTreeMap::new();
    for ((p, t), xs) in &values {
        let e: &mut Value = by_dose
            .entry(t.clone())
            .or_insert_with(|| json!({"dose_mg": dose_of.get(t.as_str())}));
        e[*p] = json!({"n": xs.len(), "geometric_mean": round4(gm(xs))});
    }
    out["realised_by_dose"] = json!(by_dose);
    // Power model ln(exposure) = a + b ln(dose) on individual values.
    let exposure = if multiple { "AUCTAU" } else { "AUCIFO" };
    for p in ["CMAX", exposure] {
        let pts: Vec<(f64, f64)> = values
            .iter()
            .filter(|((q, _), _)| *q == p)
            .flat_map(|((_, t), xs)| {
                let d = dose_of.get(t.as_str()).copied().unwrap_or(0.0);
                xs.iter().map(move |x| (math::ln(d), math::ln(*x)))
            })
            .collect();
        if let Some((b, lo, hi)) = slope_ci(&pts, 0.90) {
            out["power_model"][p] = json!({
                "slope": round4(b),
                "ci90": [round4(lo), round4(hi)],
                "target": pk.proportionality,
            });
        }
    }
    out
}

/// OLS slope of y on x with a two-sided `level` CI (t distribution).
fn slope_ci(pts: &[(f64, f64)], level: f64) -> Option<(f64, f64, f64)> {
    let n = pts.len() as f64;
    if pts.len() < 3 {
        return None;
    }
    let mx = pts.iter().map(|p| p.0).sum::<f64>() / n;
    let my = pts.iter().map(|p| p.1).sum::<f64>() / n;
    let sxx: f64 = pts.iter().map(|p| (p.0 - mx) * (p.0 - mx)).sum();
    if sxx <= 0.0 {
        return None;
    }
    let b = pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum::<f64>() / sxx;
    let a = my - b * mx;
    let rss: f64 = pts
        .iter()
        .map(|p| (p.1 - a - b * p.0) * (p.1 - a - b * p.0))
        .sum();
    let se = trialforge_core::math::sqrt(rss / (n - 2.0) / sxx);
    let t = t_quantile(0.5 + level / 2.0, n - 2.0);
    Some((b, b - t * se, b + t * se))
}

/// Student t quantile (Cornish–Fisher expansion; accurate to ~1e-4 for df ≥ 3).
fn t_quantile(p: f64, df: f64) -> f64 {
    let z = trialforge_core::math::norm_inv(p);
    let (z3, z5, z7) = (z * z * z, z * z * z * z * z, z * z * z * z * z * z * z);
    z + (z3 + z) / (4.0 * df)
        + (5.0 * z5 + 16.0 * z3 + 3.0 * z) / (96.0 * df * df)
        + (3.0 * z7 + 19.0 * z5 + 17.0 * z3 - 15.0 * z) / (384.0 * df * df * df)
}

/// Crossover bioequivalence from ADPP: per subject, the difference of mean log values on the
/// Test (or fed) and Reference (or fasted) periods; the GMR averages the sequence means (which
/// cancels period effects in these balanced designs), with a 90% CI from the pooled
/// within-sequence variance (df = subjects - sequences).
fn be_truth(r: &Resolved, pk: &trialforge_processes::pk::PkParams, adam: &[Dataset]) -> Value {
    use trialforge_core::math;
    let test_code = pk
        .relative_f
        .first()
        .map(|x| x.0.clone())
        .unwrap_or_default();
    let gmr = pk.relative_f.first().map_or(1.0, |x| x.1);
    let omega = math::sqrt(math::ln(1.0 + pk.iov_cv * pk.iov_cv));
    let mut out = json!({
        "analyte": pk.analyte,
        "target": {"gmr_auc": gmr, "ka_ratio": pk.relative_ka.first().map_or(1.0, |x| x.1), "intra_cv": pk.iov_cv},
        "limits": [0.8, 1.25],
    });
    let test: Vec<&str> = r
        .design
        .arms
        .iter()
        .flat_map(|a| a.periods.iter())
        .filter(|p| p.etcd == test_code)
        .map(|p| p.element.as_str())
        .collect();
    let (Some(adpp), Some(adsl)) = (
        adam.iter().find(|d| d.name == "ADPP"),
        adam.iter().find(|d| d.name == "ADSL"),
    ) else {
        return out;
    };
    let (Some(u), Some(pc), Some(trta), Some(aval)) = (
        chars(adpp, "USUBJID"),
        chars(adpp, "PARAMCD"),
        chars(adpp, "TRTA"),
        nums(adpp, "AVAL"),
    ) else {
        return out;
    };
    let seq: BTreeMap<&str, &str> = match (chars(adsl, "USUBJID"), chars(adsl, "TRTSEQA")) {
        (Some(us), Some(sq)) => us
            .iter()
            .zip(sq)
            .filter_map(|(a, b)| Some((a.as_deref()?, b.as_deref()?)))
            .collect(),
        _ => BTreeMap::new(),
    };
    for param in ["AUCLST", "AUCIFO", "CMAX"] {
        // subject → (log T values, log R values)
        let mut per: BTreeMap<&str, (Vec<f64>, Vec<f64>)> = BTreeMap::new();
        for i in 0..adpp.n_rows() {
            let (Some(x), Some(p), Some(t), Some(v)) = (
                u[i].as_deref(),
                pc[i].as_deref(),
                trta[i].as_deref(),
                aval[i],
            ) else {
                continue;
            };
            if p != param || v <= 0.0 {
                continue;
            }
            let e = per.entry(x).or_default();
            if test.contains(&t) {
                e.0.push(math::ln(v))
            } else {
                e.1.push(math::ln(v))
            }
        }
        let mut by_seq: BTreeMap<&str, Vec<f64>> = BTreeMap::new();
        let mut nt = 0.0;
        let mut nr = 0.0;
        for (x, (t, rr)) in &per {
            if t.is_empty() || rr.is_empty() {
                continue;
            }
            let Some(sq) = seq.get(x) else { continue };
            nt = t.len() as f64;
            nr = rr.len() as f64;
            let d = t.iter().sum::<f64>() / nt - rr.iter().sum::<f64>() / nr;
            by_seq.entry(sq).or_default().push(d);
        }
        let k = by_seq.len() as f64;
        let n: usize = by_seq.values().map(Vec::len).sum();
        if k < 2.0 || n as f64 <= k || by_seq.values().any(|v| v.is_empty()) {
            continue;
        }
        let means: Vec<f64> = by_seq
            .values()
            .map(|v| v.iter().sum::<f64>() / v.len() as f64)
            .collect();
        let theta = means.iter().sum::<f64>() / k;
        let ss: f64 = by_seq
            .values()
            .zip(&means)
            .map(|(v, m)| v.iter().map(|d| (d - m) * (d - m)).sum::<f64>())
            .sum();
        let df = n as f64 - k;
        let s2 = ss / df;
        let se =
            math::sqrt(s2 / (k * k) * by_seq.values().map(|v| 1.0 / v.len() as f64).sum::<f64>());
        let t = t_quantile(0.95, df);
        let (lo, hi) = (math::exp(theta - t * se), math::exp(theta + t * se));
        let sigma2_w = s2 / (1.0 / nt + 1.0 / nr);
        out[param] = json!({
            "n": n,
            "gmr": round4(math::exp(theta)),
            "ci90": [round4(lo), round4(hi)],
            "bioequivalent": lo >= 0.8 && hi <= 1.25,
            "intra_cv": round4(math::sqrt(math::exp(sigma2_w) - 1.0)),
        });
    }
    out["target"]["intra_sd_log"] = json!(round4(omega));
    out
}

/// Sigmoid Emax fit (Hill fixed at its target) to the realised arm means, weighted by n: the
/// ED50 is found on a log grid and refined by golden section; E0 and Emax are then linear.
fn dose_response(r: &Resolved, groups: &BTreeMap<String, Vec<f64>>) -> Value {
    use trialforge_core::math;
    let hill = r
        .controls
        .get("hill")
        .and_then(Value::as_f64)
        .unwrap_or(1.0);
    let ed50 = r
        .controls
        .get("ed50")
        .and_then(Value::as_f64)
        .unwrap_or(1.0);
    let pts: Vec<(f64, f64, f64)> = r
        .design
        .arms
        .iter()
        .filter_map(|a| {
            let xs = groups.get(&a.arm)?;
            (!xs.is_empty()).then(|| {
                (
                    a.treatment.dose,
                    xs.iter().sum::<f64>() / xs.len() as f64,
                    xs.len() as f64,
                )
            })
        })
        .collect();
    let mut out = json!({"target": {"ed50": ed50, "hill": hill}});
    let max_dose = pts.iter().map(|p| p.0).fold(0.0, f64::max);
    if pts.len() < 3 || max_dose <= 0.0 {
        return out;
    }
    let shape = |d: f64, e: f64| {
        if d <= 0.0 {
            0.0
        } else {
            let dh = math::pow(d, hill);
            dh / (math::pow(e, hill) + dh)
        }
    };
    // For a given ED50: weighted least squares of mean on [1, shape] → (E0, Emax, RSS).
    let fit = |e: f64| {
        let (mut sw, mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for &(d, y, w) in &pts {
            let x = shape(d, e);
            sw += w;
            sx += w * x;
            sy += w * y;
            sxx += w * x * x;
            sxy += w * x * y;
        }
        let det = sw * sxx - sx * sx;
        let emax = (sw * sxy - sx * sy) / det;
        let e0 = (sy - emax * sx) / sw;
        let rss: f64 = pts
            .iter()
            .map(|&(d, y, w)| {
                let r = y - e0 - emax * shape(d, e);
                w * r * r
            })
            .sum();
        (e0, emax, rss)
    };
    let (lo, hi) = (math::ln(max_dose / 1000.0), math::ln(max_dose * 10.0));
    let grid: Vec<f64> = (0..=400)
        .map(|k| lo + (hi - lo) * f64::from(k) / 400.0)
        .collect();
    let best = grid
        .iter()
        .copied()
        .min_by(|a, b| fit(math::exp(*a)).2.total_cmp(&fit(math::exp(*b)).2))
        .expect("non-empty grid");
    let step = (hi - lo) / 400.0;
    let (mut a, mut b) = (best - step, best + step);
    let phi = (math::sqrt(5.0) - 1.0) / 2.0;
    for _ in 0..60 {
        let (c, d) = (b - phi * (b - a), a + phi * (b - a));
        if fit(math::exp(c)).2 < fit(math::exp(d)).2 {
            b = d;
        } else {
            a = c;
        }
    }
    let e = math::exp((a + b) / 2.0);
    let (e0, emax, _) = fit(e);
    out["fitted"] = json!({"e0": round4(e0), "emax": round4(emax), "ed50": round4(e)});
    out["arms"] = json!(
        pts.iter()
            .map(|p| json!({"dose": p.0, "mean_change": round4(p.1), "n": p.2}))
            .collect::<Vec<_>>()
    );
    out
}

/// Non-inferiority of the last arm against the first (lower is better): the 95% CI of the
/// difference in mean change (pooled variance) against the margin.
fn non_inferiority(arms: &[String], groups: &BTreeMap<String, Vec<f64>>, margin: f64) -> Value {
    use trialforge_core::math;
    let (Some(ctrl), Some(test)) = (
        arms.first().and_then(|a| groups.get(a)),
        arms.last().and_then(|a| groups.get(a)),
    ) else {
        return json!({"margin": margin});
    };
    let stats = |xs: &[f64]| {
        let n = xs.len() as f64;
        let m = xs.iter().sum::<f64>() / n;
        (n, m, xs.iter().map(|x| (x - m) * (x - m)).sum::<f64>())
    };
    let ((n1, m1, ss1), (n2, m2, ss2)) = (stats(ctrl), stats(test));
    if n1 < 2.0 || n2 < 2.0 {
        return json!({"margin": margin});
    }
    let df = n1 + n2 - 2.0;
    let se = math::sqrt((ss1 + ss2) / df * (1.0 / n1 + 1.0 / n2));
    let diff = m2 - m1;
    let t = t_quantile(0.975, df);
    let (lo, hi) = (diff - t * se, diff + t * se);
    json!({
        "comparison": format!("{} minus {}", arms.last().expect("checked"), arms.first().expect("checked")),
        "margin": margin,
        "difference": round4(diff),
        "ci95": [round4(lo), round4(hi)],
        "non_inferior": hi < margin,
    })
}

/// ADAS-Cog: the mean change in the total at the end of treatment by arm (ADADAS, LOCF), against
/// the targets.
fn adas_truth(r: &Resolved, p: &trialforge_processes::adas::AdasParams, adam: &[Dataset]) -> Value {
    let targets: BTreeMap<&str, f64> = r
        .design
        .arms
        .iter()
        .zip(&p.change_by_arm)
        .map(|(a, c)| (a.arm.as_str(), *c))
        .collect();
    let mut out = json!({"target_change": targets, "target_baseline_mean": p.baseline_mean});
    let Some(d) = adam.iter().find(|d| d.name == "ADADAS") else {
        return out;
    };
    let (Some(pc), Some(trt), Some(anl), Some(avn), Some(chg), Some(base), Some(ablfl)) = (
        chars(d, "PARAMCD"),
        chars(d, "TRTP"),
        chars(d, "ANL01FL"),
        nums(d, "AVISITN"),
        nums(d, "CHG"),
        nums(d, "BASE"),
        chars(d, "ABLFL"),
    ) else {
        return out;
    };
    let last = avn
        .iter()
        .flatten()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let dtype = chars(d, "DTYPE");
    let mut observed: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let mut by: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let mut bases: Vec<f64> = Vec::new();
    for i in 0..d.n_rows() {
        if pc[i].as_deref() != Some("ACTOT") {
            continue;
        }
        if ablfl[i].as_deref() == Some("Y") {
            bases.extend(base[i]);
        }
        if anl[i].as_deref() == Some("Y") && avn[i] == Some(last) {
            if let (Some(t), Some(c)) = (trt[i].clone(), chg[i]) {
                if dtype.is_none_or(|x| x[i].is_none()) {
                    observed.entry(t.clone()).or_default().push(c);
                }
                by.entry(t).or_default().push(c);
            }
        }
    }
    let realised: BTreeMap<String, Value> = by
        .into_iter()
        .map(|(a, xs)| {
            (
                a,
                json!({"n": xs.len(), "mean_change": mean(&xs).map(round4)}),
            )
        })
        .collect();
    out["realised_change_locf"] = json!(realised);
    let observed: BTreeMap<String, Value> = observed
        .into_iter()
        .map(|(a, xs)| {
            (
                a,
                json!({"n": xs.len(), "mean_change": mean(&xs).map(round4)}),
            )
        })
        .collect();
    out["realised_change_observed"] = json!(observed);
    out["realised_baseline_mean"] = json!(mean(&bases).map(round4));
    out
}

/// Cardiovascular outcomes: the event-driven cutoff and, per endpoint, KM medians and the
/// stratified Cox HR (Drug A vs placebo) from ADTTE, against the targets.
fn outcomes_truth(
    r: &Resolved,
    p: &trialforge_processes::outcomes::OutcomesParams,
    adam: &[Dataset],
) -> Value {
    let arms: Vec<String> = r.design.arms.iter().map(|a| a.arm.clone()).collect();
    let hr = p.events.iter().find(|e| e.mace).map_or(1.0, |e| e.hr[0]);
    let mut out = json!({
        "target": {
            "first_mace_events": p.target_events,
            "hr": hr,
            "mace_rate_per_year": p.events.iter().filter(|e| e.mace).map(|e| e.rate).sum::<f64>(),
        },
    });
    let Some(tte) = adam.iter().find(|d| d.name == "ADTTE") else {
        return out;
    };
    if let (Some(pc), Some(cnsr)) = (chars(tte, "PARAMCD"), nums(tte, "CNSR")) {
        let events = (0..tte.n_rows())
            .filter(|&i| pc[i].as_deref() == Some("MACE") && cnsr[i] == Some(0.0))
            .count();
        out["realised_first_mace_events"] = json!(events);
    }
    for (param, key) in [
        ("MACE", "mace"),
        ("CVDEATH", "cv_death"),
        ("MI", "mi"),
        ("STROKE", "stroke"),
        ("DEATH", "all_cause_death"),
    ] {
        if let Some(v) = tte_summary(tte, param, &arms) {
            out[key] = v;
        }
    }
    out
}
