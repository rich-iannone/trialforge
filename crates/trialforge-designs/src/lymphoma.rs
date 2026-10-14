//! `heme_lymphoma`: a randomised, open-label study of Drug A versus investigator's choice in
//! relapsed or refractory diffuse large B-cell lymphoma, stratified by IPI risk. Tumour
//! assessments every 8 weeks follow the Lugano 2014 classification: nodal and extranodal
//! lesions measured by perpendicular diameters (TR), an FDG-PET Deauville score, and CT- and
//! PET-CT-based responses (RS); ADRS, PFS and OS as for solid tumours.
//!
//! The survival model, targets and calibration are the solid-tumour archetype's (see
//! `oncology`); the response target is the objective response rate on PET-CT (complete or
//! partial metabolic response, no confirmation under Lugano).

use serde::{Deserialize, Serialize};
use trialforge_calibrate::oncology::Calibration;
use trialforge_model::Design;
use trialforge_processes::oncology::OncologyParams;

use crate::oncology::{Variant, resolve_with};
use crate::{ControlDoc, StudyOptions};

/// Controls for `heme_lymphoma`: those of `oncology_solid_tumor`, with lymphoma defaults.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Controls(pub crate::oncology::Controls);

impl Default for Controls {
    fn default() -> Self {
        Self(crate::oncology::Controls {
            os_median_months: 12.0,
            os_hr: 0.7,
            pfs_median_months: 4.0,
            pfs_hr: 0.6,
            orr: [0.6, 0.35],
            ..crate::oncology::Controls::default()
        })
    }
}

pub(crate) fn docs() -> Vec<ControlDoc> {
    crate::oncology::docs()
        .into_iter()
        .map(|mut d| {
            if d.name == "orr" {
                d.description = "Objective response rate on PET-CT (Lugano 2014: CMR or PMR, no confirmation): [Drug A, Investigator's Choice]";
            }
            if d.name.ends_with("months") && d.name.starts_with("os") || d.name == "pfs_median_months" {
                d.description = match d.name {
                    "os_median_months" => "Median overall survival of the investigator's choice arm, in months",
                    _ => "Median progression-free survival of the investigator's choice arm, in months",
                };
            }
            d
        })
        .collect()
}

pub(crate) fn resolve(opts: &StudyOptions, c: &Controls) -> (Design, OncologyParams, Calibration) {
    resolve_with(opts, &c.0, Variant::Lymphoma)
}
