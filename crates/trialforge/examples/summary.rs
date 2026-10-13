//! Prints a study's datasets and truth manifest.
//!
//! `cargo run --example summary -- parallel_continuous 300 1 '{"indication":"hypertension"}'`

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let archetype = args.first().map_or("parallel_continuous", String::as_str);
    let n = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(300);
    let seed = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1);
    let mut opts = trialforge::StudyOptions::new(n, seed);
    if let Some(c) = args.get(3) {
        opts.controls = serde_json::from_str(c).expect("controls must be JSON");
    }
    let t0 = std::time::Instant::now();
    let study = trialforge::study(archetype, &opts).unwrap_or_else(|e| panic!("{e}"));
    let elapsed = t0.elapsed();
    for d in study.sdtm.iter().chain(&study.adam) {
        println!(
            "{:<7} {:>7} rows  {:>3} cols  {}",
            d.name,
            d.n_rows(),
            d.columns().len(),
            d.label
        );
    }
    println!("generated in {elapsed:?}");
    let report = trialforge::check::run(&study.sdtm, &study.adam);
    println!("{}", report.summary());
    if std::env::var("TF_TRUTH").is_err() {
        return;
    }
    println!("{}", serde_json::to_string_pretty(&study.truth).unwrap());
}
