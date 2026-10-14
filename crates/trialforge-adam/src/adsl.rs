//! ADSL from DM, EX, DS and VS (and TA for treatment numbering).

use std::collections::HashMap;

use trialforge_cdisc::build;
use trialforge_core::{Dataset, Date, Values};

use crate::{DeriveError, Frame, dt, find, require};

/// What ADSL takes from outside SDTM.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AdslInputs {
    /// `USUBJID` → randomised stratum (`STRAT1R`), from the randomisation system.
    pub strata: HashMap<String, String>,
    /// The data cutoff (`DCUTDT`), for event-driven studies.
    pub cutoff: Option<Date>,
}

#[derive(Default)]
struct Facts {
    first_dose: Option<Date>,
    last_dose: Option<Date>,
    randomized: Option<Date>,
    eot: Option<String>,
    eos: Option<(String, Option<Date>)>,
    height: Option<f64>,
    weight: Option<f64>,
}

/// Derives ADSL.
pub fn adsl(sdtm: &[Dataset], inputs: &AdslInputs) -> Result<Dataset, DeriveError> {
    let dm = Frame::new(require(sdtm, "DM")?);
    let usubjid = dm.chars("USUBJID")?;
    let mut facts: HashMap<&str, Facts> = HashMap::new();

    // EX: first and last exposure (admiral: derive_vars_merged on EXSTDTC / EXENDTC).
    if let Some(ex) = find(sdtm, "EX").map(Frame::new) {
        for ((u, st), en) in ex
            .chars("USUBJID")?
            .iter()
            .zip(ex.chars("EXSTDTC")?)
            .zip(ex.chars("EXENDTC")?)
        {
            let Some(u) = u.as_deref() else { continue };
            let f = facts.entry(u).or_default();
            if let Some(d) = dt(st) {
                f.first_dose = Some(f.first_dose.map_or(d, |x| x.min(d)));
            }
            if let Some(d) = dt(en) {
                f.last_dose = Some(f.last_dose.map_or(d, |x| x.max(d)));
            }
        }
    }

    // DS: randomisation, end of treatment, end of study (by DSSCAT).
    let ds = Frame::new(require(sdtm, "DS")?);
    let (du, ddecod, dscat, ddtc) = (
        ds.chars("USUBJID")?,
        ds.chars("DSDECOD")?,
        ds.chars("DSSCAT")?,
        ds.chars("DSSTDTC")?,
    );
    for i in 0..ds.n() {
        let Some(u) = du[i].as_deref() else { continue };
        let f = facts.entry(u).or_default();
        match (ddecod[i].as_deref(), dscat[i].as_deref()) {
            (Some("RANDOMIZED"), _) => f.randomized = dt(&ddtc[i]),
            (Some(d), Some("END OF TREATMENT")) => f.eot = Some(d.to_string()),
            (Some(d), Some("END OF STUDY")) => f.eos = Some((d.to_string(), dt(&ddtc[i]))),
            _ => {}
        }
    }

    // VS: baseline height and weight (the last pre-exposure records).
    if let Some(vs) = find(sdtm, "VS").map(Frame::new) {
        let (vu, vt, vv, vfl) = (
            vs.chars("USUBJID")?,
            vs.chars("VSTESTCD")?,
            vs.nums("VSSTRESN")?,
            vs.chars("VSLOBXFL")?,
        );
        for i in 0..vs.n() {
            let Some(u) = vu[i].as_deref() else { continue };
            if vfl[i].as_deref() != Some("Y") {
                continue;
            }
            let f = facts.entry(u).or_default();
            match vt[i].as_deref() {
                Some("HEIGHT") => f.height = vv[i],
                Some("WEIGHT") => f.weight = vv[i],
                _ => {}
            }
        }
    }

    // TA: treatment numbering follows arm order of first appearance.
    let mut arm_order: Vec<String> = Vec::new();
    if let Some(ta) = find(sdtm, "TA").map(Frame::new) {
        for a in ta.chars("ARMCD")?.iter().flatten() {
            if !arm_order.contains(a) {
                arm_order.push(a.clone());
            }
        }
    }
    let arm_n = |cd: &Option<String>| {
        cd.as_ref()
            .and_then(|c| arm_order.iter().position(|a| a == c))
            .map(|p| p as f64 + 1.0)
    };

    let none = Facts::default();
    let per: Vec<&Facts> = usubjid
        .iter()
        .map(|u| u.as_deref().and_then(|u| facts.get(u)).unwrap_or(&none))
        .collect();
    let copy = |name: &str| -> Result<Values, DeriveError> {
        Ok(dm
            .ds
            .column(name)
            .ok_or_else(|| DeriveError::Input(format!("DM has no {name}")))?
            .values
            .clone())
    };
    let age = dm.nums("AGE")?;
    let (armcd, actarmcd) = (dm.chars("ARMCD")?, dm.chars("ACTARMCD")?);
    let (arm, actarm) = (dm.chars("ARM")?, dm.chars("ACTARM")?);
    let (dthdtc, rfpendtc) = (dm.chars("DTHDTC")?, dm.chars("RFPENDTC")?);
    let n = dm.n();
    let flag = |b: bool| Some(if b { "Y" } else { "N" }.to_string());
    let status = |d: Option<&str>| {
        Some(
            match d {
                None => "ONGOING",
                Some("COMPLETED") => "COMPLETED",
                Some(_) => "DISCONTINUED",
            }
            .to_string(),
        )
    };
    let reason = |d: Option<&str>| d.filter(|d| *d != "COMPLETED").map(str::to_string);
    let bmi = |f: &Facts| match (f.height, f.weight) {
        (Some(h), Some(w)) if h > 0.0 => Some(trialforge_core::math::round_dp(
            w / (h / 100.0) / (h / 100.0),
            1,
        )),
        _ => None,
    };

    let mut cols: Vec<(&str, Values)> = vec![
        ("STUDYID", copy("STUDYID")?),
        ("USUBJID", copy("USUBJID")?),
        ("SUBJID", copy("SUBJID")?),
        ("SITEID", copy("SITEID")?),
        ("COUNTRY", copy("COUNTRY")?),
        ("AGE", copy("AGE")?),
        ("AGEU", copy("AGEU")?),
        (
            "AGEGR1",
            Values::Char(
                age.iter()
                    .map(|a| a.map(|a| if a < 65.0 { "<65" } else { ">=65" }.to_string()))
                    .collect(),
            ),
        ),
        (
            "AGEGR1N",
            Values::Num(
                age.iter()
                    .map(|a| a.map(|a| if a < 65.0 { 1.0 } else { 2.0 }))
                    .collect(),
            ),
        ),
        ("SEX", copy("SEX")?),
        ("RACE", copy("RACE")?),
        ("ETHNIC", copy("ETHNIC")?),
        ("ARMCD", copy("ARMCD")?),
        ("ARM", copy("ARM")?),
        ("ACTARMCD", copy("ACTARMCD")?),
        ("ACTARM", copy("ACTARM")?),
        ("TRT01P", Values::Char(arm.to_vec())),
        ("TRT01PN", Values::Num(armcd.iter().map(arm_n).collect())),
        ("TRT01A", Values::Char(actarm.to_vec())),
        ("TRT01AN", Values::Num(actarmcd.iter().map(arm_n).collect())),
        (
            "ITTFL",
            Values::Char(per.iter().map(|f| flag(f.randomized.is_some())).collect()),
        ),
        (
            "SAFFL",
            Values::Char(per.iter().map(|f| flag(f.first_dose.is_some())).collect()),
        ),
        (
            "COMPLFL",
            Values::Char(
                per.iter()
                    .map(|f| flag(f.eos.as_ref().is_some_and(|(d, _)| d == "COMPLETED")))
                    .collect(),
            ),
        ),
        (
            "RANDDT",
            Values::Date(per.iter().map(|f| f.randomized).collect()),
        ),
        (
            "TRTSDT",
            Values::Date(per.iter().map(|f| f.first_dose).collect()),
        ),
        (
            "TRTEDT",
            Values::Date(per.iter().map(|f| f.last_dose).collect()),
        ),
        (
            "TRTDURD",
            Values::Num(
                per.iter()
                    .map(|f| match (f.first_dose, f.last_dose) {
                        (Some(a), Some(b)) => Some(f64::from(b.days_since(a) + 1)),
                        _ => None,
                    })
                    .collect(),
            ),
        ),
        (
            "EOTSTT",
            Values::Char(per.iter().map(|f| status(f.eot.as_deref())).collect()),
        ),
        (
            "DCTREAS",
            Values::Char(per.iter().map(|f| reason(f.eot.as_deref())).collect()),
        ),
        (
            "EOSSTT",
            Values::Char(
                per.iter()
                    .map(|f| status(f.eos.as_ref().map(|(d, _)| d.as_str())))
                    .collect(),
            ),
        ),
        (
            "DCSREAS",
            Values::Char(
                per.iter()
                    .map(|f| reason(f.eos.as_ref().map(|(d, _)| d.as_str())))
                    .collect(),
            ),
        ),
        (
            "EOSDT",
            Values::Date(
                per.iter()
                    .map(|f| f.eos.as_ref().and_then(|(_, d)| *d))
                    .collect(),
            ),
        ),
        ("DTHFL", copy("DTHFL")?),
        ("DTHDT", Values::Date(dthdtc.iter().map(dt).collect())),
        // Last known alive: the later of the last dose and the end of participation.
        (
            "LSTALVDT",
            Values::Date(
                per.iter()
                    .zip(rfpendtc)
                    .map(|(f, p)| f.last_dose.max(dt(p)))
                    .collect(),
            ),
        ),
    ];
    if per.iter().any(|f| f.height.is_some()) {
        cols.push((
            "HEIGHTBL",
            Values::Num(per.iter().map(|f| f.height).collect()),
        ));
        cols.push((
            "WEIGHTBL",
            Values::Num(per.iter().map(|f| f.weight).collect()),
        ));
        cols.push(("BMIBL", Values::Num(per.iter().map(|f| bmi(f)).collect())));
    }
    if !inputs.strata.is_empty() {
        cols.push((
            "STRAT1R",
            Values::Char(
                usubjid
                    .iter()
                    .map(|u| u.as_ref().and_then(|u| inputs.strata.get(u).cloned()))
                    .collect(),
            ),
        ));
    }
    if let Some(c) = inputs.cutoff {
        cols.push(("DCUTDT", Values::Date(vec![Some(c); n])));
    }
    if let (Some(ta), Some(se)) = (find(sdtm, "TA"), find(sdtm, "SE")) {
        crossover_periods(&mut cols, &dm, ta, se)?;
    }
    Ok(build("ADSL", cols)?)
}

/// Period variables of a crossover (TA has `TREATMENT k` epochs): `TRTkkP/PN/A/AN` from TA
/// elements, `APkkSDT`/`APkkEDT` from the subject's SE, `TRTSEQP`/`TRTSEQA` from DM.
fn crossover_periods(
    cols: &mut Vec<(&str, Values)>,
    dm: &Frame<'_>,
    ta: &Dataset,
    se: &Dataset,
) -> Result<(), DeriveError> {
    let ta = Frame::new(ta);
    let (ta_arm, ta_epoch, ta_etcd, ta_el) = (
        ta.chars("ARMCD")?,
        ta.chars("EPOCH")?,
        ta.chars("ETCD")?,
        ta.chars("ELEMENT")?,
    );
    let period = |e: &Option<String>| {
        e.as_deref()
            .and_then(|e| e.strip_prefix("TREATMENT "))
            .and_then(|k| k.parse::<i64>().ok())
    };
    let n_periods = (0..ta.n())
        .filter_map(|i| period(&ta_epoch[i]))
        .max()
        .unwrap_or(0);
    if n_periods == 0 {
        return Ok(());
    }
    // Treatment elements numbered in order of first appearance in TA.
    let mut etcds: Vec<&str> = Vec::new();
    let mut element: HashMap<(&str, i64), (String, f64)> = HashMap::new();
    for i in 0..ta.n() {
        if let (Some(k), Some(arm), Some(cd)) = (
            period(&ta_epoch[i]),
            ta_arm[i].as_deref(),
            ta_etcd[i].as_deref(),
        ) {
            if !etcds.contains(&cd) {
                etcds.push(cd);
            }
            let num = etcds.iter().position(|c| *c == cd).expect("added") as f64 + 1.0;
            element.insert((arm, k), (ta_el[i].clone().unwrap_or_default(), num));
        }
    }
    let se = Frame::new(se);
    let (se_u, se_epoch, se_st, se_en) = (
        se.chars("USUBJID")?,
        se.chars("EPOCH")?,
        se.chars("SESTDTC")?,
        se.chars("SEENDTC")?,
    );
    // (start, end) of each subject's element in each period.
    type Span = (Option<Date>, Option<Date>);
    let mut spans: HashMap<(&str, i64), Span> = HashMap::new();
    for i in 0..se.n() {
        if let (Some(u), Some(k)) = (se_u[i].as_deref(), period(&se_epoch[i])) {
            spans.insert((u, k), (dt(&se_st[i]), dt(&se_en[i])));
        }
    }
    let usubjid = dm.chars("USUBJID")?;
    let (armcd, actarmcd) = (dm.chars("ARMCD")?, dm.chars("ACTARMCD")?);
    let (arm, actarm) = (dm.chars("ARM")?, dm.chars("ACTARM")?);
    let name =
        |cd: &Option<String>, k: i64| cd.as_deref().and_then(|c| element.get(&(c, k)).cloned());
    let names: [[&str; 4]; 4] = [
        ["TRT01P", "TRT01PN", "TRT01A", "TRT01AN"],
        ["TRT02P", "TRT02PN", "TRT02A", "TRT02AN"],
        ["TRT03P", "TRT03PN", "TRT03A", "TRT03AN"],
        ["TRT04P", "TRT04PN", "TRT04A", "TRT04AN"],
    ];
    let starts = ["AP01SDT", "AP02SDT", "AP03SDT", "AP04SDT"];
    let ends = ["AP01EDT", "AP02EDT", "AP03EDT", "AP04EDT"];
    for k in 1..=n_periods.min(4) {
        let v = names[(k - 1) as usize];
        let col = |cd: &[Option<String>]| {
            Values::Char(cd.iter().map(|c| name(c, k).map(|x| x.0)).collect())
        };
        let num = |cd: &[Option<String>]| {
            Values::Num(cd.iter().map(|c| name(c, k).map(|x| x.1)).collect())
        };
        let set = |cols: &mut Vec<(&str, Values)>, key: &'static str, values: Values| match cols
            .iter_mut()
            .find(|(n, _)| *n == key)
        {
            Some(slot) => slot.1 = values,
            None => cols.push((key, values)),
        };
        let span = |u: &Option<String>| u.as_deref().and_then(|u| spans.get(&(u, k)).copied());
        // Actual treatment only for periods the subject entered.
        let entered: Vec<bool> = usubjid.iter().map(|u| span(u).is_some()).collect();
        let only_entered = |v: Values| match v {
            Values::Char(x) => Values::Char(
                x.into_iter()
                    .zip(&entered)
                    .map(|(x, e)| x.filter(|_| *e))
                    .collect(),
            ),
            Values::Num(x) => Values::Num(
                x.into_iter()
                    .zip(&entered)
                    .map(|(x, e)| x.filter(|_| *e))
                    .collect(),
            ),
            other => other,
        };
        set(cols, v[0], col(armcd));
        set(cols, v[1], num(armcd));
        set(cols, v[2], only_entered(col(actarmcd)));
        set(cols, v[3], only_entered(num(actarmcd)));
        cols.push((
            starts[(k - 1) as usize],
            Values::Date(usubjid.iter().map(|u| span(u).and_then(|s| s.0)).collect()),
        ));
        cols.push((
            ends[(k - 1) as usize],
            Values::Date(usubjid.iter().map(|u| span(u).and_then(|s| s.1)).collect()),
        ));
    }
    cols.push(("TRTSEQP", Values::Char(arm.to_vec())));
    cols.push(("TRTSEQA", Values::Char(actarm.to_vec())));
    Ok(())
}
