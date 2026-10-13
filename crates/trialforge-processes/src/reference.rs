//! Reference data the processes draw on (`assets/*.tsv`): the lab panel, vital signs,
//! background AE rates, medical-history prevalence and which medicines treat what.
//!
//! These are plausible defaults for an adult trial population, not calibrated epidemiology.
//! Archetypes override what matters for their indication (e.g. HbA1c in diabetes).

use std::sync::OnceLock;

use trialforge_cdisc::read_tsv;

/// Whether an analyte varies on the natural or the log scale.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scale {
    /// Additive variation.
    Linear,
    /// Multiplicative (log-normal) variation.
    Log,
}

/// A lab analyte.
#[derive(Clone, Debug, PartialEq)]
pub struct Analyte {
    /// `LBTESTCD`.
    pub testcd: String,
    /// `LBTEST`.
    pub test: String,
    /// `LBCAT`.
    pub category: String,
    /// `LBSPEC`.
    pub specimen: String,
    /// Standard (SI) unit, `LBSTRESU`.
    pub si_unit: String,
    /// Original (conventional) unit, `LBORRESU`.
    pub conv_unit: String,
    /// SI = conventional × factor + offset.
    pub factor: f64,
    /// See `factor`.
    pub offset: f64,
    /// Decimals reported in SI.
    pub si_dp: u32,
    /// Decimals reported in conventional units.
    pub conv_dp: u32,
    /// Scale of variation.
    pub scale: Scale,
    /// Population mean (SI) for (female, male).
    pub mean: (f64, f64),
    /// Between-subject SD on the model scale (SI units, or log units for `Scale::Log`).
    pub sd_between: f64,
    /// Within-subject (visit-to-visit) SD on the model scale.
    pub sd_within: f64,
    /// Reference range (SI) for females.
    pub range_f: (f64, f64),
    /// Reference range (SI) for males.
    pub range_m: (f64, f64),
    /// Correlation group: analytes sharing a group share part of their subject effect.
    pub group: Option<String>,
}

/// A vital sign.
#[derive(Clone, Debug, PartialEq)]
pub struct VitalSign {
    /// `VSTESTCD`.
    pub testcd: String,
    /// `VSTEST`.
    pub test: String,
    /// Unit (original and standard).
    pub unit: String,
    /// Decimals reported.
    pub dp: u32,
    /// Population mean for (female, male).
    pub mean: (f64, f64),
    /// Between-subject SD.
    pub sd_between: f64,
    /// Within-subject SD.
    pub sd_within: f64,
    /// `VSPOS`, if recorded.
    pub position: Option<String>,
    /// Correlation group.
    pub group: Option<String>,
}

/// Background incidence and character of an adverse event.
#[derive(Clone, Debug, PartialEq)]
pub struct AeBackground {
    /// Preferred Term.
    pub pt: String,
    /// Events per subject-year.
    pub rate_per_year: f64,
    /// Median duration in days.
    pub median_days: f64,
    /// Probability an event is serious.
    pub p_serious: f64,
    /// Severity probabilities (mild, moderate, severe).
    pub severity: [f64; 3],
}

/// Prevalence of a medical-history condition.
#[derive(Clone, Debug, PartialEq)]
pub struct HistoryCondition {
    /// Preferred Term.
    pub pt: String,
    /// `MHCAT`.
    pub category: String,
    /// Prevalence at age 55 for (female, male).
    pub prevalence: (f64, f64),
    /// Change in log-odds per decade of age.
    pub log_odds_per_decade: f64,
    /// Whether the condition is usually ongoing.
    pub chronic: bool,
    /// Onset is up to this many years before consent.
    pub max_years: f64,
}

/// Whether a medicine treats a chronic condition or an acute event.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TreatmentKind {
    /// Long-term treatment of a medical-history condition.
    Chronic,
    /// Short treatment of an adverse event.
    Acute,
}

/// The medicines used for an indication.
#[derive(Clone, Debug, PartialEq)]
pub struct Indication {
    /// Preferred Term treated.
    pub pt: String,
    /// Chronic or acute.
    pub kind: TreatmentKind,
    /// Probability a subject with the indication takes one of the medicines.
    pub p_treat: f64,
    /// (ingredient, weight).
    pub ingredients: Vec<(String, f64)>,
}

fn num(s: &str) -> f64 {
    s.parse().unwrap_or_else(|_| panic!("not a number: {s:?}"))
}

fn opt(s: &str) -> Option<String> {
    (!s.is_empty()).then(|| s.to_string())
}

/// The lab panel.
pub fn lab_panel() -> &'static [Analyte] {
    static L: OnceLock<Vec<Analyte>> = OnceLock::new();
    L.get_or_init(|| {
        read_tsv(include_str!("../assets/lab_panel.tsv"))
            .into_iter()
            .map(|r| Analyte {
                testcd: r["testcd"].into(),
                test: r["test"].into(),
                category: r["category"].into(),
                specimen: r["spec"].into(),
                si_unit: r["si_unit"].into(),
                conv_unit: r["conv_unit"].into(),
                factor: num(r["factor"]),
                offset: num(r["offset"]),
                si_dp: num(r["si_dp"]) as u32,
                conv_dp: num(r["conv_dp"]) as u32,
                scale: if r["scale"] == "log" {
                    Scale::Log
                } else {
                    Scale::Linear
                },
                mean: (num(r["mean_f"]), num(r["mean_m"])),
                sd_between: num(r["sd_between"]),
                sd_within: num(r["sd_within"]),
                range_f: (num(r["lo_f"]), num(r["hi_f"])),
                range_m: (num(r["lo_m"]), num(r["hi_m"])),
                group: opt(r["group"]),
            })
            .collect()
    })
}

/// The vital signs.
pub fn vital_signs() -> &'static [VitalSign] {
    static V: OnceLock<Vec<VitalSign>> = OnceLock::new();
    V.get_or_init(|| {
        read_tsv(include_str!("../assets/vital_signs.tsv"))
            .into_iter()
            .map(|r| VitalSign {
                testcd: r["testcd"].into(),
                test: r["test"].into(),
                unit: r["unit"].into(),
                dp: num(r["dp"]) as u32,
                mean: (num(r["mean_f"]), num(r["mean_m"])),
                sd_between: num(r["sd_between"]),
                sd_within: num(r["sd_within"]),
                position: opt(r["position"]),
                group: opt(r["group"]),
            })
            .collect()
    })
}

/// Background AE rates.
pub fn ae_background() -> &'static [AeBackground] {
    static A: OnceLock<Vec<AeBackground>> = OnceLock::new();
    A.get_or_init(|| {
        read_tsv(include_str!("../assets/ae_background.tsv"))
            .into_iter()
            .map(|r| AeBackground {
                pt: r["pt"].into(),
                rate_per_year: num(r["rate_per_year"]),
                median_days: num(r["median_days"]),
                p_serious: num(r["p_serious"]),
                severity: [num(r["p_mild"]), num(r["p_moderate"]), num(r["p_severe"])],
            })
            .collect()
    })
}

/// Medical-history prevalence.
pub fn history_conditions() -> &'static [HistoryCondition] {
    static H: OnceLock<Vec<HistoryCondition>> = OnceLock::new();
    H.get_or_init(|| {
        read_tsv(include_str!("../assets/mh_prevalence.tsv"))
            .into_iter()
            .map(|r| HistoryCondition {
                pt: r["pt"].into(),
                category: r["category"].into(),
                prevalence: (num(r["prev_f"]), num(r["prev_m"])),
                log_odds_per_decade: num(r["log_odds_per_decade"]),
                chronic: r["chronic"] == "Y",
                max_years: num(r["max_years"]),
            })
            .collect()
    })
}

/// Indications and their medicines.
pub fn indications() -> &'static [Indication] {
    static I: OnceLock<Vec<Indication>> = OnceLock::new();
    I.get_or_init(|| {
        read_tsv(include_str!("../assets/cm_indications.tsv"))
            .into_iter()
            .map(|r| Indication {
                pt: r["indication"].into(),
                kind: if r["kind"] == "chronic" {
                    TreatmentKind::Chronic
                } else {
                    TreatmentKind::Acute
                },
                p_treat: num(r["p_treat"]),
                ingredients: r["ingredients"]
                    .split(';')
                    .zip(r["weights"].split(';'))
                    .map(|(i, w)| (i.to_string(), num(w)))
                    .collect(),
            })
            .collect()
    })
}

/// The medicines for `pt` of the given kind.
pub fn indication(pt: &str, kind: TreatmentKind) -> Option<&'static Indication> {
    indications().iter().find(|i| i.pt == pt && i.kind == kind)
}

#[cfg(test)]
mod tests {
    use super::*;
    use trialforge_cdisc::{ct, dictionaries};

    #[test]
    fn lab_panel_is_controlled() {
        let ct = ct();
        for a in lab_panel() {
            let cd = ct.codelist("LBTESTCD").unwrap();
            let pos = cd
                .terms
                .iter()
                .position(|t| *t == a.testcd)
                .unwrap_or_else(|| panic!("{}", a.testcd));
            assert_eq!(
                ct.codelist("LBTEST").unwrap().terms[pos],
                a.test,
                "{} name",
                a.testcd
            );
            for u in [&a.si_unit, &a.conv_unit] {
                assert!(
                    ct.codelist("UNIT").unwrap().contains(u),
                    "{}: unit {u}",
                    a.testcd
                );
            }
            assert!(
                ct.codelist("SPECTYPE").unwrap().contains(&a.specimen),
                "{}",
                a.testcd
            );
            assert!(
                a.range_f.0 < a.range_f.1 && a.range_m.0 < a.range_m.1,
                "{}",
                a.testcd
            );
        }
    }

    #[test]
    fn vital_signs_are_controlled() {
        let ct = ct();
        for v in vital_signs() {
            let cd = ct.codelist("VSTESTCD").unwrap();
            let pos = cd
                .terms
                .iter()
                .position(|t| *t == v.testcd)
                .unwrap_or_else(|| panic!("{}", v.testcd));
            assert_eq!(ct.codelist("VSTEST").unwrap().terms[pos], v.test);
            assert!(ct.codelist("UNIT").unwrap().contains(&v.unit), "{}", v.unit);
            if let Some(p) = &v.position {
                assert!(ct.codelist("POSITION").unwrap().contains(p));
            }
        }
    }

    #[test]
    fn terms_and_medicines_resolve() {
        for a in ae_background() {
            assert!(dictionaries::ae_term(&a.pt).is_some(), "AE {}", a.pt);
            assert!(
                (a.severity.iter().sum::<f64>() - 1.0).abs() < 1e-9,
                "{} severity sums to 1",
                a.pt
            );
        }
        for h in history_conditions() {
            assert!(dictionaries::ae_term(&h.pt).is_some(), "MH {}", h.pt);
        }
        for i in indications() {
            assert!(
                dictionaries::ae_term(&i.pt).is_some(),
                "indication {}",
                i.pt
            );
            for (ing, _) in &i.ingredients {
                assert!(dictionaries::ingredient(ing).is_some(), "{} → {ing}", i.pt);
            }
        }
    }
}
