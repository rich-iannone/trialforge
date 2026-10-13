//! Dataset → Arrow `RecordBatch` (feature `arrow`).
//!
//! Types: `Char` → `Utf8`, `Num` → `Float64`, `Int` → `Int64`, `Date` → `Date32`. Labels travel
//! as `label` metadata on each field and on the schema (with the dataset `name`).

use std::collections::HashMap;
use std::sync::Arc;

use arrow_array::{ArrayRef, Date32Array, Float64Array, Int64Array, RecordBatch, StringArray};
use arrow_schema::{ArrowError, DataType, Field, Schema};
use trialforge_core::{Dataset, Values};

/// The dataset as a single Arrow record batch.
pub fn to_record_batch(ds: &Dataset) -> Result<RecordBatch, ArrowError> {
    let mut fields = Vec::with_capacity(ds.columns().len());
    let mut arrays: Vec<ArrayRef> = Vec::with_capacity(ds.columns().len());
    for c in ds.columns() {
        let (dtype, array): (DataType, ArrayRef) = match &c.values {
            Values::Char(v) => (
                DataType::Utf8,
                Arc::new(StringArray::from_iter(v.iter().map(Option::as_deref))),
            ),
            Values::Num(v) => (DataType::Float64, Arc::new(Float64Array::from(v.clone()))),
            Values::Int(v) => (DataType::Int64, Arc::new(Int64Array::from(v.clone()))),
            Values::Date(v) => (
                DataType::Date32,
                Arc::new(Date32Array::from(
                    v.iter().map(|d| d.map(|d| d.days())).collect::<Vec<_>>(),
                )),
            ),
        };
        let metadata = HashMap::from([("label".to_string(), c.label.clone())]);
        fields.push(Field::new(&c.name, dtype, true).with_metadata(metadata));
        arrays.push(array);
    }
    let schema = Arc::new(Schema::new_with_metadata(
        fields,
        HashMap::from([
            ("name".to_string(), ds.name.clone()),
            ("label".to_string(), ds.label.clone()),
        ]),
    ));
    RecordBatch::try_new(schema, arrays)
}
