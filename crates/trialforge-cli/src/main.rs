//! `trialforge` command-line interface.
//!
//! ```text
//! trialforge generate spec.yaml -o out/            # formats and layout from the spec
//! trialforge generate spec.yaml -o out/ --submission
//! trialforge summary spec.yaml                     # datasets, checks and truth, no files
//! trialforge archetypes
//! trialforge controls parallel_continuous
//! ```

use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "trialforge",
    version,
    about = "Synthetic CDISC clinical-trial data from study specs"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Generate a study from a spec file (JSON or YAML) and write it.
    Generate {
        /// The study spec.
        spec: PathBuf,
        /// Output directory.
        #[arg(short, long, default_value = "trialforge-out")]
        out: PathBuf,
        /// Comma-separated formats (overrides the spec's `outputs.formats`).
        #[arg(long, value_delimiter = ',')]
        formats: Option<Vec<String>>,
        /// Write the submission layout (sdtm/ and adam/ with XPT and define.xml).
        #[arg(long)]
        submission: bool,
    },
    /// Generate a study from a spec and print its datasets and truth manifest, without writing.
    Summary {
        /// The study spec.
        spec: PathBuf,
    },
    /// List the archetypes.
    Archetypes,
    /// List an archetype's controls with their defaults.
    Controls {
        /// Archetype name.
        archetype: String,
    },
}

fn read(path: &PathBuf) -> Result<String, String> {
    std::fs::read_to_string(path).map_err(|e| format!("cannot read {}: {e}", path.display()))
}

fn run(cli: Cli) -> Result<(), String> {
    match cli.command {
        Command::Generate {
            spec,
            out,
            formats,
            submission,
        } => {
            let (study, spec) =
                trialforge::study_from_spec(&read(&spec)?).map_err(|e| e.to_string())?;
            let paths = if submission || (formats.is_none() && spec.outputs.layout == "submission")
            {
                study.write_submission(&out)
            } else {
                let formats = formats.unwrap_or(spec.outputs.formats);
                let formats: Vec<&str> = formats.iter().map(String::as_str).collect();
                study.write(&out, &formats)
            }
            .map_err(|e| e.to_string())?;
            println!(
                "{} ({}), seed {}: {} subjects, wrote {} files to {}",
                study.studyid,
                study.archetype,
                study.seed,
                study.dataset("DM").map_or(0, |d| d.n_rows()),
                paths.len(),
                out.display()
            );
        }
        Command::Summary { spec } => {
            let (study, _) =
                trialforge::study_from_spec(&read(&spec)?).map_err(|e| e.to_string())?;
            println!(
                "{} ({}), seed {}",
                study.studyid, study.archetype, study.seed
            );
            for d in study.sdtm.iter().chain(&study.adam) {
                println!("  {:<7} {:>8} rows  {}", d.name, d.n_rows(), d.label);
            }
            println!("verification: {}", trialforge::verify(&study).summary());
            println!(
                "{}",
                serde_json::to_string_pretty(&study.truth).expect("truth serializes")
            );
        }
        Command::Archetypes => {
            for (name, description) in trialforge::ARCHETYPES {
                println!("{name:<22} {description}");
            }
        }
        Command::Controls { archetype } => {
            let docs = trialforge::controls(&archetype).map_err(|e| e.to_string())?;
            for d in docs.as_array().expect("control docs are a list") {
                println!(
                    "{:<20} default {:<18} {}",
                    d["name"].as_str().unwrap_or(""),
                    d["default"].to_string(),
                    d["description"].as_str().unwrap_or("")
                );
            }
        }
    }
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
