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
    // The solid-tumor descriptions that name the comparator arm, reworded for this one.
    crate::oncology::docs()
        .into_iter()
        .map(|mut d| {
            d.description = match d.name {
                "os_median_months" => "Median overall survival of the Investigator's Choice arm, in months",
                "os_hr" => "Overall-survival hazard ratio, Drug A vs Investigator's Choice (exact proportional hazards)",
                "pfs_median_months" => "Median progression-free survival of the Investigator's Choice arm, in months (of the latent time to \
                     progression or death; observed KM medians fall on the 8-weekly assessment dates)",
                "pfs_hr" => "Progression hazard ratio, Drug A vs Investigator's Choice",
                "orr" => "Objective response rate on PET-CT (Lugano 2014: CMR or PMR, no confirmation): [Drug A, Investigator's Choice]",
                "crossover_prob" => "Probability that an Investigator's Choice subject who progresses crosses over to Drug A (crossover on progression)",
                "crossover_os_hr" => "Hazard ratio of death after crossing over, relative to staying on investigator's choice; null uses os_hr",
                _ => d.description,
            };
            d
        })
        .collect()
}

pub(crate) fn resolve(opts: &StudyOptions, c: &Controls) -> (Design, OncologyParams, Calibration) {
    resolve_with(opts, &c.0, Variant::Lymphoma)
}
