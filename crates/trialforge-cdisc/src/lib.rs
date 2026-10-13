//! The CDISC metadata bundle: SDTMIG 3.4 and ADaMIG 1.3 dataset and variable specs, a
//! version-pinned controlled terminology (CT) snapshot, and the open coding dictionaries
//! ([`dictionaries`]).
//!
//! The bundle is plain TSV under `assets/`, embedded at compile time so it travels with the
//! crate (including when vendored for R). It is trialforge's source of truth; a test
//! cross-checks its labels against Pointblank's templates.
//!
//! - `assets/sdtmig-3.4/{datasets,variables}.tsv`, `assets/adamig-1.3/…`: dataset label, class,
//!   structure and keys; per-variable label, type (`Char`/`Num`), core (`Req`/`Exp`/`Perm`, and
//!   `Cond` for ADaM) and codelist.
//! - `assets/ct/2024-09-27.tsv`: a curated subset of NCI EVS CDISC CT 2024-09-27 covering every
//!   codelist-controlled value trialforge emits. `source` records whether a term came from
//!   Pointblank's bundled subset or was curated here.

pub mod dictionaries;

/// Days after the last dose during which a new AE still counts as treatment-emergent (SUPPAE
/// `AETRTEM`, ADAE `TRTEMFL`).
pub const TEAE_WINDOW_DAYS: i32 = 30;

use std::collections::HashMap;
use std::sync::OnceLock;

use trialforge_core::{Column, Dataset, TableError, Values};

/// A variable's SAS type.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VarType {
    /// Character.
    Char,
    /// Numeric (including dates in ADaM).
    Num,
}

/// A variable's core status.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Core {
    /// Required: must be present and populated.
    Req,
    /// Expected: must be present; may be null.
    Exp,
    /// Permissible.
    Perm,
    /// Conditionally required (ADaM).
    Cond,
}

/// A variable's metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct VariableSpec {
    /// Variable name.
    pub name: String,
    /// Variable label (at most 40 characters).
    pub label: String,
    /// SAS type.
    pub var_type: VarType,
    /// Core status.
    pub core: Core,
    /// Controlling codelist, if any (a key into [`ct`]).
    pub codelist: Option<String>,
    /// Define-XML origin type: `Collected`, `Derived`, `Assigned`, `Protocol` or
    /// `Predecessor`.
    pub origin: String,
    /// Origin source (`Investigator`, `Sponsor`, `Vendor`, `Subject`) or, for predecessors, the
    /// source variable (e.g. `DM.AGE`). Empty if not applicable.
    pub source: String,
    /// How a derived variable is computed (Define-XML method). Empty if not derived.
    pub method: String,
}

/// A dataset's metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct DatasetSpec {
    /// Dataset name.
    pub name: String,
    /// Dataset label.
    pub label: String,
    /// Observation class, e.g. `FINDINGS` or `BASIC DATA STRUCTURE`.
    pub class: String,
    /// Structure, e.g. `One record per subject`.
    pub structure: String,
    /// Natural key variables.
    pub keys: Vec<String>,
    /// Variables in standard order.
    pub variables: Vec<VariableSpec>,
}

/// An implementation guide's datasets.
#[derive(Clone, Debug, PartialEq)]
pub struct Standard {
    /// `SDTMIG` or `ADaMIG`.
    pub name: &'static str,
    /// Version, e.g. `3.4`.
    pub version: &'static str,
    /// Datasets in bundle order.
    pub datasets: Vec<DatasetSpec>,
}

/// A codelist.
#[derive(Clone, Debug, PartialEq)]
pub struct Codelist {
    /// Whether sponsors may add terms.
    pub extensible: bool,
    /// Submission values, in CT order.
    pub terms: Vec<String>,
}

impl Codelist {
    /// Whether `value` is a term of this codelist.
    pub fn contains(&self, value: &str) -> bool {
        self.terms.iter().any(|t| t == value)
    }
}

/// A controlled terminology package.
#[derive(Clone, Debug, PartialEq)]
pub struct Terminology {
    /// Package date, e.g. `2024-09-27`.
    pub version: &'static str,
    /// Codelists by submission name (e.g. `SEX`, `UNIT`).
    pub codelists: HashMap<String, Codelist>,
}

impl Terminology {
    /// The codelist called `name`.
    pub fn codelist(&self, name: &str) -> Option<&Codelist> {
        self.codelists.get(name)
    }
}

/// Rows of a tab-separated file with a header row, as maps from column name to value. Missing
/// trailing cells are empty strings. Lines starting with `#` are comments.
pub fn read_tsv(text: &'static str) -> Vec<HashMap<&'static str, &'static str>> {
    let mut lines = text
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'));
    let header: Vec<&str> = lines.next().expect("TSV header").split('\t').collect();
    lines
        .map(|l| {
            let cells: Vec<&str> = l.split('\t').collect();
            assert!(
                cells.len() <= header.len(),
                "too many cells in TSV row: {l}"
            );
            header
                .iter()
                .zip(cells.into_iter().chain(std::iter::repeat("")))
                .map(|(h, c)| (*h, c))
                .collect()
        })
        .collect()
}

fn load_standard(
    name: &'static str,
    version: &'static str,
    datasets: &'static str,
    variables: &'static str,
) -> Standard {
    let vars = read_tsv(variables);
    let datasets = read_tsv(datasets)
        .into_iter()
        .map(|d| DatasetSpec {
            name: d["dataset"].to_string(),
            label: d["label"].to_string(),
            class: d["class"].to_string(),
            structure: d["structure"].to_string(),
            keys: d["keys"].split(',').map(str::to_string).collect(),
            variables: vars
                .iter()
                .filter(|v| v["dataset"] == d["dataset"])
                .map(|v| VariableSpec {
                    name: v["variable"].to_string(),
                    label: v["label"].to_string(),
                    var_type: match v["type"] {
                        "Char" => VarType::Char,
                        "Num" => VarType::Num,
                        t => panic!("bad type {t} for {}", v["variable"]),
                    },
                    core: match v["core"] {
                        "Req" => Core::Req,
                        "Exp" => Core::Exp,
                        "Perm" => Core::Perm,
                        "Cond" => Core::Cond,
                        c => panic!("bad core {c} for {}", v["variable"]),
                    },
                    codelist: Some(v["codelist"])
                        .filter(|c| !c.is_empty())
                        .map(str::to_string),
                    origin: v["origin"].to_string(),
                    source: v["source"].to_string(),
                    method: v["method"].to_string(),
                })
                .collect(),
        })
        .collect();
    Standard {
        name,
        version,
        datasets,
    }
}

/// SDTMIG 3.4.
pub fn sdtmig() -> &'static Standard {
    static S: OnceLock<Standard> = OnceLock::new();
    S.get_or_init(|| {
        load_standard(
            "SDTMIG",
            "3.4",
            include_str!("../assets/sdtmig-3.4/datasets.tsv"),
            include_str!("../assets/sdtmig-3.4/variables.tsv"),
        )
    })
}

/// ADaMIG 1.3.
pub fn adamig() -> &'static Standard {
    static S: OnceLock<Standard> = OnceLock::new();
    S.get_or_init(|| {
        load_standard(
            "ADaMIG",
            "1.3",
            include_str!("../assets/adamig-1.3/datasets.tsv"),
            include_str!("../assets/adamig-1.3/variables.tsv"),
        )
    })
}

/// The CDISC CT snapshot (2024-09-27).
pub fn ct() -> &'static Terminology {
    static T: OnceLock<Terminology> = OnceLock::new();
    T.get_or_init(|| {
        let mut codelists: HashMap<String, Codelist> = HashMap::new();
        for r in read_tsv(include_str!("../assets/ct/2024-09-27.tsv")) {
            let cl = codelists
                .entry(r["codelist"].to_string())
                .or_insert_with(|| Codelist {
                    extensible: r["extensible"] == "Y",
                    terms: Vec::new(),
                });
            cl.terms.push(r["term"].to_string());
        }
        Terminology {
            version: "2024-09-27",
            codelists,
        }
    })
}

/// The spec for dataset `name`, from SDTMIG or ADaMIG.
pub fn spec(name: &str) -> Option<&'static DatasetSpec> {
    sdtmig()
        .datasets
        .iter()
        .chain(&adamig().datasets)
        .find(|d| d.name == name)
}

/// Errors from building a dataset against its spec.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MetadataError {
    /// A dataset the bundle does not define.
    UnknownDataset(String),
    /// A variable the spec does not define.
    UnknownVariable {
        /// The dataset.
        dataset: String,
        /// The unknown variable.
        variable: String,
    },
    /// The columns could not form a dataset.
    Table(TableError),
}

impl core::fmt::Display for MetadataError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            MetadataError::UnknownDataset(d) => {
                write!(f, "{d} is not a dataset in the metadata bundle")
            }
            MetadataError::UnknownVariable { dataset, variable } => {
                write!(
                    f,
                    "{variable} is not a {dataset} variable in the metadata bundle"
                )
            }
            MetadataError::Table(e) => e.fmt(f),
        }
    }
}

impl std::error::Error for MetadataError {}

/// Builds dataset `name` from (variable, values) pairs: labels come from the bundle and columns
/// are put in standard order. Bundle variables without values are left out; values for
/// variables the bundle does not define are an error.
pub fn build(name: &str, values: Vec<(&str, Values)>) -> Result<Dataset, MetadataError> {
    spec(name)
        .ok_or_else(|| MetadataError::UnknownDataset(name.to_string()))?
        .build(values)
}

impl DatasetSpec {
    /// The variable called `name`.
    pub fn variable(&self, name: &str) -> Option<&VariableSpec> {
        self.variables.iter().find(|v| v.name == name)
    }

    /// See [`build`].
    pub fn build(&self, mut values: Vec<(&str, Values)>) -> Result<Dataset, MetadataError> {
        if let Some((unknown, _)) = values.iter().find(|(n, _)| self.variable(n).is_none()) {
            return Err(MetadataError::UnknownVariable {
                dataset: self.name.clone(),
                variable: (*unknown).to_string(),
            });
        }
        let mut columns = Vec::with_capacity(values.len());
        for v in &self.variables {
            if let Some(pos) = values.iter().position(|(n, _)| *n == v.name) {
                let (_, vals) = values.swap_remove(pos);
                columns.push(Column::new(v.name.clone(), v.label.clone(), vals));
            }
        }
        Dataset::new(self.name.clone(), self.label.clone(), columns).map_err(MetadataError::Table)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bundle_fits_sas_transport_v5() {
        for spec in sdtmig().datasets.iter().chain(&adamig().datasets) {
            assert!(
                spec.name.len() <= 8 && spec.label.len() <= 40,
                "{}",
                spec.name
            );
            assert!(!spec.variables.is_empty(), "{} has variables", spec.name);
            for k in &spec.keys {
                assert!(
                    spec.variable(k).is_some(),
                    "{} key {k} is a variable",
                    spec.name
                );
            }
            for v in &spec.variables {
                assert!(
                    v.name.len() <= 8 && v.name.is_ascii(),
                    "{}.{}",
                    spec.name,
                    v.name
                );
                assert!(
                    v.label.len() <= 40 && v.label.is_ascii(),
                    "{}.{}: {}",
                    spec.name,
                    v.name,
                    v.label
                );
            }
        }
    }

    #[test]
    fn every_variable_has_a_define_origin() {
        let types = [
            "Collected",
            "Derived",
            "Assigned",
            "Protocol",
            "Predecessor",
        ];
        for spec in sdtmig().datasets.iter().chain(&adamig().datasets) {
            for v in &spec.variables {
                assert!(
                    types.contains(&v.origin.as_str()),
                    "{}.{}: {:?}",
                    spec.name,
                    v.name,
                    v.origin
                );
                assert_eq!(
                    v.origin == "Derived",
                    !v.method.is_empty() && v.origin != "Assigned",
                    "{}.{} method",
                    spec.name,
                    v.name
                );
                if v.origin == "Predecessor" {
                    let (ds, var) = v
                        .source
                        .split_once('.')
                        .unwrap_or_else(|| panic!("{}.{}", spec.name, v.name));
                    let src = super::spec(ds)
                        .unwrap_or_else(|| panic!("{}.{} → {ds}", spec.name, v.name));
                    assert!(
                        src.variable(var).is_some(),
                        "{}.{} → {}",
                        spec.name,
                        v.name,
                        v.source
                    );
                }
            }
        }
    }

    #[test]
    fn every_codelist_reference_resolves() {
        for spec in sdtmig().datasets.iter().chain(&adamig().datasets) {
            for v in &spec.variables {
                if let Some(cl) = &v.codelist {
                    assert!(
                        ct().codelist(cl).is_some(),
                        "{}.{} → missing codelist {cl}",
                        spec.name,
                        v.name
                    );
                }
            }
        }
    }

    #[test]
    fn paired_codelists_line_up() {
        for (cd, name) in [
            ("LBTESTCD", "LBTEST"),
            ("VSTESTCD", "VSTEST"),
            ("TSPARMCD", "TSPARM"),
        ] {
            let (a, b) = (ct().codelist(cd).unwrap(), ct().codelist(name).unwrap());
            assert_eq!(a.terms.len(), b.terms.len(), "{cd}/{name}");
        }
    }

    #[test]
    fn build_orders_columns_and_rejects_unknowns() {
        let ds = build(
            "DM",
            vec![
                ("USUBJID", Values::Char(vec![Some("X-1".into())])),
                ("STUDYID", Values::Char(vec![Some("X".into())])),
            ],
        )
        .unwrap();
        let names: Vec<&str> = ds.columns().iter().map(|c| c.name.as_str()).collect();
        assert_eq!(names, ["STUDYID", "USUBJID"]);
        assert_eq!(ds.columns()[1].label, "Unique Subject Identifier");
        assert!(matches!(
            build("DM", vec![("NOPE", Values::Int(vec![]))]),
            Err(MetadataError::UnknownVariable { .. })
        ));
        assert!(matches!(
            build("XX", vec![]),
            Err(MetadataError::UnknownDataset(_))
        ));
    }
}
