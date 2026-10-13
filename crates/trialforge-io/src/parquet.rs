//! Parquet writer (feature `parquet`): one file per dataset, uncompressed with dictionary
//! encoding, labels kept as Arrow field metadata. Output is deterministic for a given
//! `parquet` crate version (the writer records its version in the file).

use std::path::Path;

use parquet::arrow::ArrowWriter;
use parquet::file::properties::WriterProperties;
use trialforge_core::Dataset;

/// Errors from writing Parquet.
#[derive(Debug)]
pub enum ParquetError {
    /// The dataset could not be converted or encoded.
    Encode(String),
    /// Writing failed.
    Io(std::io::Error),
}

impl core::fmt::Display for ParquetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ParquetError::Encode(m) => f.write_str(m),
            ParquetError::Io(e) => e.fmt(f),
        }
    }
}

impl std::error::Error for ParquetError {}

/// Encodes `dataset` as Parquet bytes.
pub fn to_bytes(dataset: &Dataset) -> Result<Vec<u8>, ParquetError> {
    let batch =
        crate::arrow::to_record_batch(dataset).map_err(|e| ParquetError::Encode(e.to_string()))?;
    let props = WriterProperties::builder()
        .set_created_by("trialforge".to_string())
        .build();
    let mut out = Vec::new();
    let mut w = ArrowWriter::try_new(&mut out, batch.schema(), Some(props))
        .map_err(|e| ParquetError::Encode(e.to_string()))?;
    w.write(&batch)
        .map_err(|e| ParquetError::Encode(e.to_string()))?;
    w.close().map_err(|e| ParquetError::Encode(e.to_string()))?;
    Ok(out)
}

/// Writes `dataset` to `path` as Parquet.
pub fn write_file(dataset: &Dataset, path: &Path) -> Result<(), ParquetError> {
    let bytes = to_bytes(dataset)?;
    std::fs::write(path, bytes).map_err(ParquetError::Io)
}
