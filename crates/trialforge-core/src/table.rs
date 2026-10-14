//! In-memory datasets.
//!
//! A [`Dataset`] is a named, labelled set of equal-length [`Column`]s. It is deliberately
//! small: the crates that emit and derive data build `Dataset`s, and the writers and language
//! bindings convert them to XPT, Arrow, R data frames and so on.

use crate::date::{Date, DateTime};

/// The values of one column. `None` is a missing value.
#[derive(Clone, Debug, PartialEq)]
pub enum Values {
    /// Character data (SAS `Char`).
    Char(Vec<Option<String>>),
    /// Floating-point numbers (SAS `Num`).
    Num(Vec<Option<f64>>),
    /// Whole numbers, such as `AGE`, `--SEQ` and `--DY` (SAS `Num`).
    Int(Vec<Option<i64>>),
    /// Calendar dates, such as ADaM `*DT` variables (SAS `Num` with a `DATE9.` format).
    Date(Vec<Option<Date>>),
    /// Date-times, such as ADaM `*DTM` variables (SAS `Num` with a `DATETIME20.` format).
    DateTime(Vec<Option<DateTime>>),
}

impl Values {
    /// The number of values.
    pub fn len(&self) -> usize {
        match self {
            Values::Char(v) => v.len(),
            Values::Num(v) => v.len(),
            Values::Int(v) => v.len(),
            Values::Date(v) => v.len(),
            Values::DateTime(v) => v.len(),
        }
    }

    /// Whether there are no values.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// A named, labelled column.
#[derive(Clone, Debug, PartialEq)]
pub struct Column {
    /// Variable name, e.g. `USUBJID`.
    pub name: String,
    /// Variable label, e.g. `Unique Subject Identifier`.
    pub label: String,
    /// The values.
    pub values: Values,
}

impl Column {
    /// A new column.
    pub fn new(name: impl Into<String>, label: impl Into<String>, values: Values) -> Self {
        Self {
            name: name.into(),
            label: label.into(),
            values,
        }
    }

    /// The character values, if this is a character column.
    pub fn as_char(&self) -> Option<&[Option<String>]> {
        match &self.values {
            Values::Char(v) => Some(v),
            _ => None,
        }
    }

    /// The integer values, if this is an integer column.
    pub fn as_int(&self) -> Option<&[Option<i64>]> {
        match &self.values {
            Values::Int(v) => Some(v),
            _ => None,
        }
    }
}

/// A named, labelled dataset (an SDTM domain or an ADaM dataset).
#[derive(Clone, Debug, PartialEq)]
pub struct Dataset {
    /// Dataset name, e.g. `DM` or `ADSL`.
    pub name: String,
    /// Dataset label, e.g. `Demographics`.
    pub label: String,
    columns: Vec<Column>,
    n_rows: usize,
}

impl Dataset {
    /// A dataset from columns, which must all have the same length and distinct names.
    pub fn new(
        name: impl Into<String>,
        label: impl Into<String>,
        columns: Vec<Column>,
    ) -> Result<Self, TableError> {
        let name = name.into();
        let n_rows = columns.first().map_or(0, |c| c.values.len());
        for (i, c) in columns.iter().enumerate() {
            if c.values.len() != n_rows {
                return Err(TableError::LengthMismatch {
                    dataset: name,
                    column: c.name.clone(),
                    expected: n_rows,
                    found: c.values.len(),
                });
            }
            if columns[..i].iter().any(|o| o.name == c.name) {
                return Err(TableError::DuplicateColumn {
                    dataset: name,
                    column: c.name.clone(),
                });
            }
        }
        Ok(Self {
            name,
            label: label.into(),
            columns,
            n_rows,
        })
    }

    /// The columns, in order.
    pub fn columns(&self) -> &[Column] {
        &self.columns
    }

    /// The column called `name`.
    pub fn column(&self, name: &str) -> Option<&Column> {
        self.columns.iter().find(|c| c.name == name)
    }

    /// The number of rows.
    pub fn n_rows(&self) -> usize {
        self.n_rows
    }
}

/// Errors from building a [`Dataset`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TableError {
    /// A column's length differs from the first column's.
    LengthMismatch {
        /// The dataset being built.
        dataset: String,
        /// The offending column.
        column: String,
        /// The first column's length.
        expected: usize,
        /// This column's length.
        found: usize,
    },
    /// Two columns share a name.
    DuplicateColumn {
        /// The dataset being built.
        dataset: String,
        /// The repeated name.
        column: String,
    },
}

impl core::fmt::Display for TableError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            TableError::LengthMismatch {
                dataset,
                column,
                expected,
                found,
            } => write!(
                f,
                "{dataset}.{column} has {found} values but the dataset has {expected} rows"
            ),
            TableError::DuplicateColumn { dataset, column } => {
                write!(f, "{dataset} has more than one column named {column}")
            }
        }
    }
}

impl std::error::Error for TableError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_ragged_and_duplicate_columns() {
        let a = Column::new("A", "a", Values::Int(vec![Some(1), Some(2)]));
        let b = Column::new("B", "b", Values::Int(vec![Some(1)]));
        assert!(matches!(
            Dataset::new("X", "x", vec![a.clone(), b]),
            Err(TableError::LengthMismatch { .. })
        ));
        assert!(matches!(
            Dataset::new("X", "x", vec![a.clone(), a]),
            Err(TableError::DuplicateColumn { .. })
        ));
    }
}
