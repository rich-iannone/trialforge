//! The CDISC metadata bundle: SDTM and ADaM variable specs and, later, version-pinned
//! controlled terminology.
//!
//! The spike carries the variables of DM, EX, DS (SDTMIG 3.4) and ADSL (ADaMIG 1.3) as static
//! tables. Labels follow the implementation guides (and Pointblank's templates). Phase 1
//! replaces these tables with a JSON bundle shared with Pointblank.

use trialforge_core::{Column, Dataset, TableError, Values};

/// A dataset's metadata: name, label and variables in standard order.
#[derive(Clone, Copy, Debug)]
pub struct DatasetSpec {
    /// Dataset name.
    pub name: &'static str,
    /// Dataset label.
    pub label: &'static str,
    /// (variable name, label), in standard order.
    pub variables: &'static [(&'static str, &'static str)],
}

/// SDTM DM (Demographics).
pub const DM: DatasetSpec = DatasetSpec {
    name: "DM",
    label: "Demographics",
    variables: &[
        ("STUDYID", "Study Identifier"),
        ("DOMAIN", "Domain Abbreviation"),
        ("USUBJID", "Unique Subject Identifier"),
        ("SUBJID", "Subject Identifier for the Study"),
        ("RFSTDTC", "Subject Reference Start Date/Time"),
        ("RFENDTC", "Subject Reference End Date/Time"),
        ("RFXSTDTC", "Date/Time of First Study Treatment"),
        ("RFXENDTC", "Date/Time of Last Study Treatment"),
        ("RFICDTC", "Date/Time of Informed Consent"),
        ("RFPENDTC", "Date/Time of End of Participation"),
        ("DTHDTC", "Date/Time of Death"),
        ("DTHFL", "Subject Death Flag"),
        ("SITEID", "Study Site Identifier"),
        ("AGE", "Age"),
        ("AGEU", "Age Units"),
        ("SEX", "Sex"),
        ("RACE", "Race"),
        ("ETHNIC", "Ethnicity"),
        ("ARMCD", "Planned Arm Code"),
        ("ARM", "Description of Planned Arm"),
        ("ACTARMCD", "Actual Arm Code"),
        ("ACTARM", "Description of Actual Arm"),
        ("COUNTRY", "Country"),
        ("DMDTC", "Date/Time of Collection"),
        ("DMDY", "Study Day of Collection"),
    ],
};

/// SDTM EX (Exposure).
pub const EX: DatasetSpec = DatasetSpec {
    name: "EX",
    label: "Exposure",
    variables: &[
        ("STUDYID", "Study Identifier"),
        ("DOMAIN", "Domain Abbreviation"),
        ("USUBJID", "Unique Subject Identifier"),
        ("EXSEQ", "Sequence Number"),
        ("EXTRT", "Name of Treatment"),
        ("EXDOSE", "Dose"),
        ("EXDOSU", "Dose Units"),
        ("EXDOSFRM", "Dose Form"),
        ("EXDOSFRQ", "Dosing Frequency per Interval"),
        ("EXROUTE", "Route of Administration"),
        ("EPOCH", "Epoch"),
        ("EXSTDTC", "Start Date/Time of Treatment"),
        ("EXENDTC", "End Date/Time of Treatment"),
        ("EXSTDY", "Study Day of Start of Treatment"),
        ("EXENDY", "Study Day of End of Treatment"),
    ],
};

/// SDTM DS (Disposition).
pub const DS: DatasetSpec = DatasetSpec {
    name: "DS",
    label: "Disposition",
    variables: &[
        ("STUDYID", "Study Identifier"),
        ("DOMAIN", "Domain Abbreviation"),
        ("USUBJID", "Unique Subject Identifier"),
        ("DSSEQ", "Sequence Number"),
        ("DSTERM", "Reported Term for the Disposition Event"),
        ("DSDECOD", "Standardized Disposition Term"),
        ("DSCAT", "Category for Disposition Event"),
        ("EPOCH", "Epoch"),
        ("DSSTDTC", "Start Date/Time of Disposition Event"),
        ("DSSTDY", "Study Day of Start of Disposition Event"),
    ],
};

/// ADaM ADSL (Subject-Level Analysis Dataset).
pub const ADSL: DatasetSpec = DatasetSpec {
    name: "ADSL",
    label: "Subject-Level Analysis Dataset",
    variables: &[
        ("STUDYID", "Study Identifier"),
        ("USUBJID", "Unique Subject Identifier"),
        ("SUBJID", "Subject Identifier for the Study"),
        ("SITEID", "Study Site Identifier"),
        ("COUNTRY", "Country"),
        ("AGE", "Age"),
        ("AGEU", "Age Units"),
        ("AGEGR1", "Pooled Age Group 1"),
        ("SEX", "Sex"),
        ("RACE", "Race"),
        ("ETHNIC", "Ethnicity"),
        ("ARMCD", "Planned Arm Code"),
        ("ARM", "Description of Planned Arm"),
        ("ACTARMCD", "Actual Arm Code"),
        ("ACTARM", "Description of Actual Arm"),
        ("TRT01P", "Planned Treatment for Period 01"),
        ("TRT01A", "Actual Treatment for Period 01"),
        ("STRAT1R", "Stratification Factor 1 (Randomized)"),
        ("ITTFL", "Intent-To-Treat Population Flag"),
        ("SAFFL", "Safety Population Flag"),
        ("RANDDT", "Date of Randomization"),
        ("TRTSDT", "Date of First Exposure to Treatment"),
        ("TRTEDT", "Date of Last Exposure to Treatment"),
        ("EOTSTT", "End of Treatment Status"),
        ("DCTREAS", "Reason for Discontinuation of Treatment"),
        ("EOSSTT", "End of Study Status"),
        ("DCSREAS", "Reason for Discontinuation from Study"),
        ("EOSDT", "End of Study Date"),
        ("DTHFL", "Subject Death Flag"),
        ("DTHDT", "Date of Death"),
        ("LSTALVDT", "Date Last Known Alive"),
        ("DCUTDT", "Data Cutoff Date"),
    ],
};

/// Every dataset spec in the bundle.
pub const ALL: &[DatasetSpec] = &[DM, EX, DS, ADSL];

/// Errors from building a dataset against its spec.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MetadataError {
    /// A variable the spec does not define.
    UnknownVariable {
        /// The dataset.
        dataset: &'static str,
        /// The unknown variable.
        variable: String,
    },
    /// The columns could not form a dataset.
    Table(TableError),
}

impl core::fmt::Display for MetadataError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
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

impl DatasetSpec {
    /// Builds a dataset from (name, values) pairs: labels come from the spec and columns are
    /// put in standard order. Spec variables without values are left out (permissible
    /// variables are optional); values for variables the spec does not define are an error.
    pub fn build(&self, mut values: Vec<(&str, Values)>) -> Result<Dataset, MetadataError> {
        if let Some((unknown, _)) = values
            .iter()
            .find(|(n, _)| !self.variables.iter().any(|(v, _)| v == n))
        {
            return Err(MetadataError::UnknownVariable {
                dataset: self.name,
                variable: (*unknown).to_string(),
            });
        }
        let mut columns = Vec::with_capacity(values.len());
        for (name, label) in self.variables {
            if let Some(pos) = values.iter().position(|(n, _)| n == name) {
                let (_, v) = values.swap_remove(pos);
                columns.push(Column::new(*name, *label, v));
            }
        }
        Dataset::new(self.name, self.label, columns).map_err(MetadataError::Table)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_and_labels_fit_sas_transport_v5() {
        for spec in ALL {
            assert!(
                spec.name.len() <= 8 && spec.label.len() <= 40,
                "{}",
                spec.name
            );
            for (name, label) in spec.variables {
                assert!(name.len() <= 8, "{}.{name}", spec.name);
                assert!(label.len() <= 40, "{}.{name}: {label}", spec.name);
            }
        }
    }

    #[test]
    fn build_orders_columns_and_rejects_unknowns() {
        let ds = DM
            .build(vec![
                ("USUBJID", Values::Char(vec![Some("X-1".into())])),
                ("STUDYID", Values::Char(vec![Some("X".into())])),
            ])
            .unwrap();
        let names: Vec<&str> = ds.columns().iter().map(|c| c.name.as_str()).collect();
        assert_eq!(names, ["STUDYID", "USUBJID"]);
        assert_eq!(ds.columns()[1].label, "Unique Subject Identifier");
        assert!(matches!(
            DM.build(vec![("NOPE", Values::Int(vec![]))]),
            Err(MetadataError::UnknownVariable { .. })
        ));
    }
}
