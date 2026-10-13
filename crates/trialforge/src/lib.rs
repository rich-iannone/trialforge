//! trialforge: synthetic CDISC clinical-trial data from declarative study specs.
//!
//! This facade crate re-exports the workspace crates. The Python and R bindings and the CLI
//! depend on it alone, so the internal crate layout can change without touching them.

pub use trialforge_adam as adam;
pub use trialforge_calibrate as calibrate;
pub use trialforge_cdisc as cdisc;
pub use trialforge_check as check;
pub use trialforge_core as core;
pub use trialforge_designs as designs;
pub use trialforge_hooks as hooks;
pub use trialforge_io as io;
pub use trialforge_model as model;
pub use trialforge_processes as processes;
pub use trialforge_sdtm as sdtm;

/// The version of the trialforge core.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
