//! Open coding dictionaries.
//!
//! - [`ae_terms`]: an open, MedDRA-shaped dictionary of MedDRA System Organ Class names and
//!   common Preferred Terms (British spelling, primary SOC). It has no LLT/HLT/HLGT levels and
//!   no MedDRA codes, and it is not MedDRA. Used to code AE and MH.
//! - [`ingredients`]: RxNorm ingredient names (NLM, public domain) with the WHO ATC level-2
//!   class name and a typical dose. Used to code CM.
//!
//! Licensed dictionaries (MedDRA, WHODrug) are bring-your-own and come in a later phase.

use std::sync::OnceLock;

use crate::read_tsv;

/// Name and version recorded wherever a coding dictionary is cited (e.g. define.xml).
pub const AE_DICTIONARY: (&str, &str) = ("TRIALFORGE OPEN MEDDRA-SHAPED TERMS", "0.1");
/// Name and version of the con-med dictionary.
pub const CM_DICTIONARY: (&str, &str) = ("TRIALFORGE OPEN RXNORM INGREDIENTS", "0.1");

/// A Preferred Term and its primary System Organ Class.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AeTerm {
    /// Preferred Term, e.g. `Headache`.
    pub pt: String,
    /// Primary System Organ Class, e.g. `Nervous system disorders`.
    pub soc: String,
}

/// A medication ingredient with its class and a typical dose.
#[derive(Clone, Debug, PartialEq)]
pub struct Ingredient {
    /// RxNorm ingredient name, upper case (`CMDECOD`).
    pub name: String,
    /// ATC level-2 class name (`CMCLAS`).
    pub atc_class: String,
    /// Typical dose per administration (`CMDOSE`).
    pub dose: f64,
    /// Dose unit (`CMDOSU`, CT `UNIT`).
    pub unit: String,
    /// Route (`CMROUTE`, CT `ROUTE`).
    pub route: String,
    /// Frequency (`CMDOSFRQ`, CT `FREQ`).
    pub frequency: String,
    /// Dose form (CT `FRM`).
    pub form: String,
}

/// Every Preferred Term in the open AE dictionary.
pub fn ae_terms() -> &'static [AeTerm] {
    static T: OnceLock<Vec<AeTerm>> = OnceLock::new();
    T.get_or_init(|| {
        read_tsv(include_str!("../assets/dictionaries/ae_terms.tsv"))
            .into_iter()
            .map(|r| AeTerm {
                pt: r["pt"].to_string(),
                soc: r["soc"].to_string(),
            })
            .collect()
    })
}

/// The term for Preferred Term `pt`.
pub fn ae_term(pt: &str) -> Option<&'static AeTerm> {
    ae_terms().iter().find(|t| t.pt == pt)
}

/// The 27 MedDRA System Organ Class names.
pub fn socs() -> &'static [String] {
    static S: OnceLock<Vec<String>> = OnceLock::new();
    S.get_or_init(|| {
        read_tsv(include_str!("../assets/dictionaries/meddra_soc.tsv"))
            .into_iter()
            .map(|r| r["soc"].to_string())
            .collect()
    })
}

/// Every ingredient in the open con-med dictionary.
pub fn ingredients() -> &'static [Ingredient] {
    static I: OnceLock<Vec<Ingredient>> = OnceLock::new();
    I.get_or_init(|| {
        read_tsv(include_str!(
            "../assets/dictionaries/rxnorm_ingredients.tsv"
        ))
        .into_iter()
        .map(|r| Ingredient {
            name: r["ingredient"].to_string(),
            atc_class: r["atc_class"].to_string(),
            dose: r["dose"].parse().expect("numeric dose"),
            unit: r["unit"].to_string(),
            route: r["route"].to_string(),
            frequency: r["frequency"].to_string(),
            form: r["form"].to_string(),
        })
        .collect()
    })
}

/// The ingredient called `name`.
pub fn ingredient(name: &str) -> Option<&'static Ingredient> {
    ingredients().iter().find(|i| i.name == name)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ct;

    #[test]
    fn ae_terms_map_to_real_socs() {
        assert_eq!(socs().len(), 27);
        assert!(ae_terms().len() >= 150);
        for t in ae_terms() {
            assert!(socs().contains(&t.soc), "{} → {}", t.pt, t.soc);
            assert!(t.pt.len() <= 200 && t.pt.is_ascii());
        }
        assert_eq!(ae_term("Headache").unwrap().soc, "Nervous system disorders");
    }

    #[test]
    fn ingredient_fields_are_controlled_terms() {
        for i in ingredients() {
            assert!(
                ct().codelist("UNIT").unwrap().contains(&i.unit),
                "{}: unit {}",
                i.name,
                i.unit
            );
            assert!(
                ct().codelist("ROUTE").unwrap().contains(&i.route),
                "{}: route {}",
                i.name,
                i.route
            );
            assert!(
                ct().codelist("FREQ").unwrap().contains(&i.frequency),
                "{}: freq {}",
                i.name,
                i.frequency
            );
            assert!(
                ct().codelist("FRM").unwrap().contains(&i.form),
                "{}: form {}",
                i.name,
                i.form
            );
            assert!(i.dose > 0.0);
        }
    }
}
