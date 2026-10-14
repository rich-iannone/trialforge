//! SAS transport (XPORT) version 5 writer, after the SAS TS-140 layout.
//!
//! One dataset per file, as FDA submissions expect. Output is byte-for-byte deterministic: the
//! header "created" and "modified" timestamps are a caller-supplied date (trialforge passes the
//! data cutoff) rather than the wall clock.
//!
//! Type mapping: `Char` → character (blank-padded; missing is all blanks); `Num`, `Int` → 8-byte
//! IBM floating point (missing is `.`); `Date` → 8-byte IBM float of the SAS date with a
//! `DATE9.` format.

use std::io::{self, Write};
use std::path::Path;

use trialforge_core::{Dataset, Date, Values};

const RECORD: usize = 80;
const NAMESTR_LEN: usize = 140;
const MAX_CHAR_LEN: usize = 200;

/// Errors from writing XPT.
#[derive(Debug)]
pub enum XptError {
    /// The dataset breaks a transport v5 limit (name or label length, character length, …).
    Invalid(String),
    /// Writing failed.
    Io(io::Error),
}

impl core::fmt::Display for XptError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            XptError::Invalid(m) => f.write_str(m),
            XptError::Io(e) => e.fmt(f),
        }
    }
}

impl std::error::Error for XptError {}

impl From<io::Error> for XptError {
    fn from(e: io::Error) -> Self {
        XptError::Io(e)
    }
}

/// Writes `dataset` to `path` as XPT v5.
pub fn write_file(dataset: &Dataset, path: &Path, timestamp: Date) -> Result<(), XptError> {
    let bytes = to_bytes(dataset, timestamp)?;
    std::fs::File::create(path)?.write_all(&bytes)?;
    Ok(())
}

/// Encodes `dataset` as XPT v5 bytes.
pub fn to_bytes(dataset: &Dataset, timestamp: Date) -> Result<Vec<u8>, XptError> {
    let layout = Layout::of(dataset)?;
    let stamp = timestamp.sas_header_datetime();
    let mut out = Vec::new();

    // Library header.
    record(
        &mut out,
        "HEADER RECORD*******LIBRARY HEADER RECORD!!!!!!!000000000000000000000000000000  ",
    );
    record(
        &mut out,
        &format!(
            "{:<8}{:<8}{:<8}{:<8}{:<8}{:24}{stamp}",
            "SAS", "SAS", "SASLIB", "9.4", "RUST", ""
        ),
    );
    record(&mut out, &format!("{stamp}{:64}", ""));

    // Member header and descriptor.
    record(
        &mut out,
        "HEADER RECORD*******MEMBER  HEADER RECORD!!!!!!!000000000000000001600000000140  ",
    );
    record(
        &mut out,
        "HEADER RECORD*******DSCRPTR HEADER RECORD!!!!!!!000000000000000000000000000000  ",
    );
    record(
        &mut out,
        &format!(
            "{:<8}{:<8}{:<8}{:<8}{:<8}{:24}{stamp}",
            "SAS", dataset.name, "SASDATA", "9.4", "RUST", ""
        ),
    );
    record(
        &mut out,
        &format!("{stamp}{:16}{:<40}{:8}", "", dataset.label, ""),
    );

    // Variable descriptors (namestr records), padded to a whole number of records.
    record(
        &mut out,
        &format!(
            "HEADER RECORD*******NAMESTR HEADER RECORD!!!!!!!000000{:04}00000000000000000000  ",
            layout.vars.len()
        ),
    );
    let start = out.len();
    for (i, v) in layout.vars.iter().enumerate() {
        out.extend_from_slice(&namestr(v, i));
    }
    pad(&mut out, start);

    // Observations.
    record(
        &mut out,
        "HEADER RECORD*******OBS     HEADER RECORD!!!!!!!000000000000000000000000000000  ",
    );
    let start = out.len();
    for row in 0..dataset.n_rows() {
        for (v, col) in layout.vars.iter().zip(dataset.columns()) {
            match &col.values {
                Values::Char(x) => {
                    let s = x[row].as_deref().unwrap_or("");
                    out.extend_from_slice(s.as_bytes());
                    out.extend(std::iter::repeat_n(b' ', v.length - s.len()));
                }
                Values::Num(x) => out.extend_from_slice(&ibm(x[row])),
                Values::Int(x) => out.extend_from_slice(&ibm(x[row].map(|n| n as f64))),
                Values::Date(x) => {
                    out.extend_from_slice(&ibm(x[row].map(|d| f64::from(d.sas_days()))))
                }
                Values::DateTime(x) => {
                    out.extend_from_slice(&ibm(x[row].map(|d| d.sas_seconds() as f64)))
                }
            }
        }
    }
    pad(&mut out, start);
    Ok(out)
}

/// A variable's transport attributes.
struct Var<'a> {
    name: &'a str,
    label: &'a str,
    is_char: bool,
    length: usize,
    format: Option<(&'static str, i16)>,
    position: usize,
}

struct Layout<'a> {
    vars: Vec<Var<'a>>,
}

impl<'a> Layout<'a> {
    fn of(ds: &'a Dataset) -> Result<Self, XptError> {
        let invalid = |m: String| Err(XptError::Invalid(m));
        if !valid_name(&ds.name) {
            return invalid(format!(
                "dataset name {:?} must be 1–8 ASCII letters, digits or _",
                ds.name
            ));
        }
        if !ds.label.is_ascii() || ds.label.len() > 40 {
            return invalid(format!(
                "{} label must be ASCII and at most 40 characters",
                ds.name
            ));
        }
        if ds.columns().len() > 9999 {
            return invalid(format!("{} has more than 9999 variables", ds.name));
        }
        let mut position = 0;
        let mut vars = Vec::with_capacity(ds.columns().len());
        for c in ds.columns() {
            if !valid_name(&c.name) {
                return invalid(format!(
                    "{}.{} must be 1–8 ASCII letters, digits or _",
                    ds.name, c.name
                ));
            }
            if !c.label.is_ascii() || c.label.len() > 40 {
                return invalid(format!(
                    "{}.{} label must be ASCII and at most 40 characters",
                    ds.name, c.name
                ));
            }
            let (is_char, length, format) = match &c.values {
                Values::Char(x) => {
                    if x.iter().flatten().any(|s| !s.is_ascii()) {
                        return invalid(format!("{}.{} has non-ASCII values", ds.name, c.name));
                    }
                    let len = x
                        .iter()
                        .flatten()
                        .map(String::len)
                        .max()
                        .unwrap_or(0)
                        .max(1);
                    if len > MAX_CHAR_LEN {
                        return invalid(format!(
                            "{}.{} has values longer than {MAX_CHAR_LEN}",
                            ds.name, c.name
                        ));
                    }
                    (true, len, None)
                }
                Values::Num(_) | Values::Int(_) => (false, 8, None),
                Values::Date(_) => (false, 8, Some(("DATE", 9))),
                Values::DateTime(_) => (false, 8, Some(("DATETIME", 20))),
            };
            vars.push(Var {
                name: &c.name,
                label: &c.label,
                is_char,
                length,
                format,
                position,
            });
            position += length;
        }
        Ok(Self { vars })
    }
}

fn valid_name(s: &str) -> bool {
    !s.is_empty() && s.len() <= 8 && s.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_')
}

fn namestr(v: &Var<'_>, index: usize) -> [u8; NAMESTR_LEN] {
    let mut b = [0u8; NAMESTR_LEN];
    let mut put = |at: usize, bytes: &[u8]| b[at..at + bytes.len()].copy_from_slice(bytes);
    put(0, &(if v.is_char { 2i16 } else { 1 }).to_be_bytes()); // ntype
    put(2, &0i16.to_be_bytes()); // nhfun
    put(4, &(v.length as i16).to_be_bytes()); // nlng
    put(6, &((index + 1) as i16).to_be_bytes()); // nvar0
    put(8, format!("{:<8}", v.name).as_bytes()); // nname
    put(16, format!("{:<40}", v.label).as_bytes()); // nlabel
    let (form, width) = v.format.unwrap_or(("", 0));
    put(56, format!("{form:<8}").as_bytes()); // nform
    put(64, &width.to_be_bytes()); // nfl
    put(66, &0i16.to_be_bytes()); // nfd
    put(68, &0i16.to_be_bytes()); // nfj
    // 70..72 nfill stays zero.
    put(72, format!("{:<8}", "").as_bytes()); // niform
    put(80, &0i16.to_be_bytes()); // nifl
    put(82, &0i16.to_be_bytes()); // nifd
    put(84, &(v.position as i32).to_be_bytes()); // npos
    // 88..140 rest stays zero.
    b
}

fn record(out: &mut Vec<u8>, s: &str) {
    debug_assert_eq!(s.len(), RECORD, "header record must be 80 bytes: {s:?}");
    out.extend_from_slice(s.as_bytes());
}

/// Blank-pads everything written since `start` to a whole number of 80-byte records.
fn pad(out: &mut Vec<u8>, start: usize) {
    let rem = (out.len() - start) % RECORD;
    if rem != 0 {
        out.extend(std::iter::repeat_n(b' ', RECORD - rem));
    }
}

/// An IEEE 754 double as an 8-byte IBM System/360 hexadecimal float. `None` is SAS missing.
fn ibm(x: Option<f64>) -> [u8; 8] {
    let Some(x) = x.filter(|x| x.is_finite()) else {
        return [b'.', 0, 0, 0, 0, 0, 0, 0];
    };
    if x == 0.0 {
        return [0; 8];
    }
    let bits = x.to_bits();
    let sign = (bits >> 63) as u8;
    let exp = ((bits >> 52) & 0x7ff) as i32;
    if exp == 0 {
        return [0; 8]; // subnormals are far below IBM's range; flush to zero
    }
    // value = m56 * 2^-56 * 2^p, with m56 the 53-bit significand moved into 56 bits.
    let mut m56 = ((bits & 0x000f_ffff_ffff_ffff) | (1 << 52)) << 3;
    let p = exp - 1022;
    let q = (p + 3).div_euclid(4); // ceil(p / 4): IBM exponents are powers of 16
    m56 >>= 4 * q - p; // shift of 0–3 bits; the low 3 bits of m56 are zero, so nothing is lost
    let e = q + 64;
    assert!((0..=127).contains(&e), "{x} is outside the IBM float range");
    let mut b = [0u8; 8];
    b[0] = (sign << 7) | e as u8;
    b[1..].copy_from_slice(&m56.to_be_bytes()[1..]);
    b
}

#[cfg(test)]
mod tests {
    use super::*;
    use trialforge_core::Column;

    #[test]
    fn ibm_known_values() {
        assert_eq!(ibm(Some(1.0)), [0x41, 0x10, 0, 0, 0, 0, 0, 0]);
        assert_eq!(ibm(Some(-1.0)), [0xC1, 0x10, 0, 0, 0, 0, 0, 0]);
        assert_eq!(ibm(Some(100.0)), [0x42, 0x64, 0, 0, 0, 0, 0, 0]);
        assert_eq!(ibm(Some(0.5)), [0x40, 0x80, 0, 0, 0, 0, 0, 0]);
        assert_eq!(
            ibm(Some(0.1)),
            [0x40, 0x19, 0x99, 0x99, 0x99, 0x99, 0x99, 0x9A]
        );
        assert_eq!(ibm(None), [0x2E, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(ibm(Some(0.0)), [0; 8]);
    }

    #[test]
    fn layout_is_record_aligned() {
        let ds = Dataset::new(
            "DM",
            "Demographics",
            vec![
                Column::new(
                    "USUBJID",
                    "Unique Subject Identifier",
                    Values::Char(vec![Some("A-1".into()), None]),
                ),
                Column::new("AGE", "Age", Values::Int(vec![Some(61), None])),
                Column::new(
                    "RFDT",
                    "Reference Date",
                    Values::Date(vec![Date::from_ymd(2024, 1, 2), None]),
                ),
            ],
        )
        .unwrap();
        let bytes = to_bytes(&ds, Date::from_ymd(2025, 1, 1).unwrap()).unwrap();
        assert_eq!(bytes.len() % RECORD, 0);
        // 8 header records, 3 namestrs (420 bytes → 6 records), OBS header, 2 rows × 19 bytes (1 record).
        assert_eq!(bytes.len(), RECORD * (8 + 6 + 1 + 1));
        assert_eq!(&bytes[..20], b"HEADER RECORD*******");
    }

    #[test]
    fn rejects_overlong_labels() {
        let ds = Dataset::new(
            "DM",
            "Demographics",
            vec![Column::new("X", "x".repeat(41), Values::Int(vec![Some(1)]))],
        )
        .unwrap();
        assert!(matches!(
            to_bytes(&ds, Date::from_days(0)),
            Err(XptError::Invalid(_))
        ));
    }
}
