//! CDISC Dataset-JSON 1.1 writer: one JSON file per dataset.
//!
//! Types: `Char` → `string` (with `length`), `Int` → `integer`, `Num` → `double`, `Date` →
//! `date` with `targetDataType: "integer"` and `displayFormat: "DATE9."` (a SAS numeric date
//! written as an ISO 8601 string). Missing values are `null`. OIDs match the define.xml written
//! by [`crate::define`], and `metaDataRef` points at it. Output is deterministic: timestamps are
//! the study date passed in.

use std::fmt::Write as _;

use trialforge_cdisc::spec;
use trialforge_core::{Dataset, Values};

use crate::define::{DefineInfo, Standard};

/// The Dataset-JSON version written.
pub const VERSION: &str = "1.1.0";

fn json(s: &str) -> String {
    serde_json::to_string(s).expect("strings serialize")
}

/// Encodes `dataset` as Dataset-JSON. `info` supplies the study, standard and date.
pub fn to_string(dataset: &Dataset, info: &DefineInfo<'_>, source_version: &str) -> String {
    let keys = spec(&dataset.name)
        .map(|s| s.keys.clone())
        .unwrap_or_default();
    let label = match info.standard {
        Standard::Sdtm => "SDTM",
        Standard::Adam => "ADaM",
    };
    let stamp = format!("{}T00:00:00", info.created);
    let mut out = String::new();
    let _ = write!(
        out,
        "{{\"datasetJSONCreationDateTime\":{},\"datasetJSONVersion\":{},\"fileOID\":{},\"dbLastModifiedDateTime\":{},\
         \"originator\":\"trialforge\",\"sourceSystem\":{{\"name\":\"trialforge\",\"version\":{}}},\"studyOID\":{},\
         \"metaDataVersionOID\":{},\"metaDataRef\":\"define.xml\",\"itemGroupOID\":{},\"records\":{},\"name\":{},\"label\":{},\"columns\":[",
        json(&stamp),
        json(VERSION),
        json(&format!("TF.{}.{}", info.studyid, dataset.name)),
        json(&stamp),
        json(source_version),
        json(&format!("ST.{}", info.studyid)),
        json(&format!("MDV.{}.{label}", info.studyid)),
        json(&format!("IG.{}", dataset.name)),
        dataset.n_rows(),
        json(&dataset.name),
        json(&dataset.label),
    );
    for (k, c) in dataset.columns().iter().enumerate() {
        if k > 0 {
            out.push(',');
        }
        let _ = write!(
            out,
            "{{\"itemOID\":{},\"name\":{},\"label\":{}",
            json(&format!("IT.{}.{}", dataset.name, c.name)),
            json(&c.name),
            json(&c.label)
        );
        match &c.values {
            Values::Char(v) => {
                let len = v
                    .iter()
                    .flatten()
                    .map(String::len)
                    .max()
                    .unwrap_or(0)
                    .max(1);
                let _ = write!(out, ",\"dataType\":\"string\",\"length\":{len}");
            }
            Values::Int(_) => out.push_str(",\"dataType\":\"integer\""),
            Values::Num(_) => out.push_str(",\"dataType\":\"double\""),
            Values::Date(_) => {
                out.push_str(",\"dataType\":\"date\",\"targetDataType\":\"integer\",\"displayFormat\":\"DATE9.\"");
            }
        }
        if let Some(p) = keys.iter().position(|key| *key == c.name) {
            let _ = write!(out, ",\"keySequence\":{}", p + 1);
        }
        out.push('}');
    }
    out.push_str("],\"rows\":[");
    for i in 0..dataset.n_rows() {
        if i > 0 {
            out.push(',');
        }
        out.push('[');
        for (k, c) in dataset.columns().iter().enumerate() {
            if k > 0 {
                out.push(',');
            }
            match &c.values {
                Values::Char(v) => {
                    out.push_str(&v[i].as_deref().map_or_else(|| "null".into(), json))
                }
                Values::Int(v) => {
                    out.push_str(&v[i].map_or_else(|| "null".into(), |n| n.to_string()))
                }
                Values::Num(v) => out.push_str(&v[i].filter(|x| x.is_finite()).map_or_else(
                    || "null".into(),
                    |x| serde_json::to_string(&x).expect("finite floats serialize"),
                )),
                Values::Date(v) => {
                    out.push_str(&v[i].map_or_else(|| "null".into(), |d| json(&d.to_string())))
                }
            }
        }
        out.push(']');
    }
    out.push_str("]}\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use trialforge_core::{Column, Date};

    #[test]
    fn writes_valid_json_with_types_and_nulls() {
        let ds = Dataset::new(
            "ADSL",
            "Subject-Level Analysis Dataset",
            vec![
                Column::new(
                    "USUBJID",
                    "Unique Subject Identifier",
                    Values::Char(vec![Some("X-1".into()), None]),
                ),
                Column::new("AGE", "Age", Values::Int(vec![Some(61), None])),
                Column::new("BMIBL", "Baseline BMI", Values::Num(vec![Some(27.5), None])),
                Column::new(
                    "TRTSDT",
                    "Start",
                    Values::Date(vec![Date::from_ymd(2024, 1, 2), None]),
                ),
            ],
        )
        .unwrap();
        let info = DefineInfo {
            studyid: "S1",
            description: "d",
            created: Date::from_ymd(2025, 1, 1).unwrap(),
            standard: Standard::Adam,
        };
        let v: serde_json::Value = serde_json::from_str(&to_string(&ds, &info, "0.1.0")).unwrap();
        assert_eq!(v["records"], 2);
        assert_eq!(v["datasetJSONVersion"], "1.1.0");
        assert_eq!(v["columns"][0]["keySequence"], 2);
        assert_eq!(v["columns"][3]["dataType"], "date");
        assert_eq!(
            v["rows"][0],
            serde_json::json!(["X-1", 61, 27.5, "2024-01-02"])
        );
        assert_eq!(v["rows"][1], serde_json::json!([null, null, null, null]));
    }
}
