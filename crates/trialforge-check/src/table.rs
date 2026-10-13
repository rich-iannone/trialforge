//! Loose, string-friendly access to datasets for the rules.

use std::collections::HashMap;

use trialforge_core::{Dataset, Date, Values};

/// A dataset with lookups that return `None` instead of failing (structure is L0's job).
pub(crate) struct T<'a> {
    pub ds: &'a Dataset,
}

impl<'a> T<'a> {
    pub fn new(ds: &'a Dataset) -> Self {
        Self { ds }
    }

    pub fn name(&self) -> &str {
        &self.ds.name
    }

    pub fn n(&self) -> usize {
        self.ds.n_rows()
    }

    pub fn has(&self, var: &str) -> bool {
        self.ds.column(var).is_some()
    }

    /// Character values; `None` if absent or not character.
    pub fn c(&self, var: &str) -> Option<&'a [Option<String>]> {
        self.ds.column(var)?.as_char()
    }

    /// Numeric values as f64 (Num, Int, and Date as days).
    pub fn f(&self, var: &str) -> Option<Vec<Option<f64>>> {
        Some(match &self.ds.column(var)?.values {
            Values::Num(v) => v.clone(),
            Values::Int(v) => v.iter().map(|x| x.map(|n| n as f64)).collect(),
            Values::Date(v) => v.iter().map(|x| x.map(|d| f64::from(d.days()))).collect(),
            Values::Char(_) => return None,
        })
    }

    /// Integer values.
    pub fn i(&self, var: &str) -> Option<Vec<Option<i64>>> {
        Some(match &self.ds.column(var)?.values {
            Values::Int(v) => v.clone(),
            Values::Num(v) => v.iter().map(|x| x.map(|n| n as i64)).collect(),
            _ => return None,
        })
    }

    /// Date values (ADaM dates, or SDTM ISO dates parsed when complete).
    pub fn d(&self, var: &str) -> Option<Vec<Option<Date>>> {
        Some(match &self.ds.column(var)?.values {
            Values::Date(v) => v.clone(),
            Values::Char(v) => v
                .iter()
                .map(|x| x.as_deref().and_then(Date::parse_iso))
                .collect(),
            _ => return None,
        })
    }

    /// Row index by USUBJID (first occurrence).
    pub fn by_subject(&self) -> HashMap<&'a str, usize> {
        let mut m = HashMap::new();
        if let Some(u) = self.c("USUBJID") {
            for (i, x) in u.iter().enumerate() {
                if let Some(x) = x.as_deref() {
                    m.entry(x).or_insert(i);
                }
            }
        }
        m
    }
}

/// The dataset called `name`.
pub(crate) fn get<'a>(sets: &'a [Dataset], name: &str) -> Option<T<'a>> {
    sets.iter().find(|d| d.name == name).map(T::new)
}

/// Whether `s` is a valid ISO 8601 date or partial date (`YYYY`, `YYYY-MM`, `YYYY-MM-DD`).
pub(crate) fn iso_ok(s: &str) -> bool {
    let b = s.as_bytes();
    let digits =
        |r: std::ops::Range<usize>| r.clone().all(|i| b.get(i).is_some_and(u8::is_ascii_digit));
    match b.len() {
        4 => digits(0..4),
        7 => {
            digits(0..4)
                && b[4] == b'-'
                && digits(5..7)
                && (1..=12).contains(&s[5..7].parse::<u32>().unwrap_or(0))
        }
        10 => Date::parse_iso(s).is_some(),
        _ => false,
    }
}
