//! Phase 3 exit criterion: the 95% CI of the stratified Cox OS hazard ratio, estimated from
//! ADTTE, covers the target `os_hr` in 93–97% of replicate studies.
//!
//! Slow (1,000 studies); run explicitly:
//! `cargo test --release -p trialforge --test coverage -- --ignored --nocapture`.

use serde_json::json;
use trialforge::{StudyOptions, study};

fn covers(truth: &serde_json::Value, endpoint: &str, target: f64) -> bool {
    let ci = &truth[endpoint]["stratified_cox_hr"]["ci95"];
    let (lo, hi) = (ci[0].as_f64().unwrap(), ci[1].as_f64().unwrap());
    lo <= target && target <= hi
}

#[test]
#[ignore = "slow: 1,000 replicate studies"]
fn os_hazard_ratio_ci_coverage() {
    let reps: u64 = std::env::var("TF_COVERAGE_REPS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(1000);
    let (mut os, mut pfs) = (0, 0);
    for seed in 0..reps {
        let opts = StudyOptions {
            n: 400,
            seed: 1_000 + seed,
            studyid: None,
            controls: json!(null),
            defects: None,
        };
        let s = study("oncology_solid_tumor", &opts).unwrap();
        os += usize::from(covers(&s.truth, "overall_survival", 0.75));
        pfs += usize::from(covers(&s.truth, "progression_free_survival", 0.65));
    }
    let (os, pfs) = (os as f64 / reps as f64, pfs as f64 / reps as f64);
    println!("OS HR 95% CI coverage {os:.3}; PFS (scan-detected, informational) {pfs:.3}");
    assert!((0.93..=0.97).contains(&os), "OS coverage {os}");
}
