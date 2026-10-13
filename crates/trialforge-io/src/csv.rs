//! CSV writer (RFC 4180): a header of variable names, `\n` line endings, fields quoted when they
//! contain a comma, quote or line break. Missing values are empty; dates are ISO 8601; numbers
//! use the shortest representation that round-trips.

use std::path::Path;

use trialforge_core::{Dataset, Values};

fn field(out: &mut String, s: &str) {
    if s.contains([',', '"', '\n', '\r']) {
        out.push('"');
        out.push_str(&s.replace('"', "\"\""));
        out.push('"');
    } else {
        out.push_str(s);
    }
}

/// Encodes `dataset` as CSV.
pub fn to_string(dataset: &Dataset) -> String {
    let cols = dataset.columns();
    let mut out = String::new();
    for (j, c) in cols.iter().enumerate() {
        if j > 0 {
            out.push(',');
        }
        field(&mut out, &c.name);
    }
    out.push('\n');
    for i in 0..dataset.n_rows() {
        for (j, c) in cols.iter().enumerate() {
            if j > 0 {
                out.push(',');
            }
            match &c.values {
                Values::Char(v) => field(&mut out, v[i].as_deref().unwrap_or("")),
                Values::Num(v) => {
                    if let Some(x) = v[i] {
                        out.push_str(&x.to_string());
                    }
                }
                Values::Int(v) => {
                    if let Some(x) = v[i] {
                        out.push_str(&x.to_string());
                    }
                }
                Values::Date(v) => {
                    if let Some(x) = v[i] {
                        out.push_str(&x.to_string());
                    }
                }
            }
        }
        out.push('\n');
    }
    out
}

/// Writes `dataset` to `path` as CSV.
pub fn write_file(dataset: &Dataset, path: &Path) -> std::io::Result<()> {
    std::fs::write(path, to_string(dataset))
}

#[cfg(test)]
mod tests {
    use super::*;
    use trialforge_core::{Column, Date};

    #[test]
    fn quoting_and_missing_values() {
        let ds = Dataset::new(
            "X",
            "x",
            vec![
                Column::new(
                    "A",
                    "a",
                    Values::Char(vec![Some("plain".into()), Some("has, comma".into()), None]),
                ),
                Column::new("B", "b", Values::Num(vec![Some(0.1), None, Some(-2.5)])),
                Column::new(
                    "C",
                    "c",
                    Values::Date(vec![Date::from_ymd(2024, 1, 2), None, None]),
                ),
            ],
        )
        .unwrap();
        assert_eq!(
            to_string(&ds),
            "A,B,C\nplain,0.1,2024-01-02\n\"has, comma\",,\n,-2.5,\n"
        );
    }
}
