//! `trialforge` command-line interface. Only `--version` exists so far;
//! `generate` and `verify` subcommands come once the engine can produce a study.

fn main() {
    match std::env::args().nth(1).as_deref() {
        Some("--version" | "-V") => println!("trialforge {}", trialforge::VERSION),
        _ => {
            eprintln!(
                "trialforge {}: no subcommands are implemented yet",
                trialforge::VERSION
            );
            std::process::exit(2);
        }
    }
}
