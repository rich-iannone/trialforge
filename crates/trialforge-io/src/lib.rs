//! Output writers: SAS XPORT v5 ([`xpt`]), CSV ([`csv`]), Arrow ([`arrow`], feature `arrow`)
//! and Parquet ([`parquet`], feature `parquet`). Dataset-JSON and define.xml come in Phase 2.

#[cfg(feature = "arrow")]
pub mod arrow;
pub mod csv;
#[cfg(feature = "parquet")]
pub mod parquet;
pub mod xpt;
