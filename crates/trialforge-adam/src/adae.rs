//! ADAE (OCCDS) from AE and ADSL; `DLTFL` from the SUPPAE `AEDLT` qualifier when present.

use std::collections::HashSet;

use trialforge_cdisc::{TEAE_WINDOW_DAYS, build};
use trialforge_core::{Dataset, Values};

use crate::{DeriveError, Frame, dt, index, rel_day};

/// Derives ADAE.
pub fn adae(
    ae: &Dataset,
    suppae: Option<&Dataset>,
    adsl: &Dataset,
) -> Result<Dataset, DeriveError> {
    let ae = Frame::new(ae);
    let sl = Frame::new(adsl);
    let sl_idx = index(sl.chars("USUBJID")?);
    let usubjid = ae.chars("USUBJID")?;
    let row: Vec<Option<usize>> = usubjid
        .iter()
        .map(|u| u.as_deref().and_then(|u| sl_idx.get(u).copied()))
        .collect();

    let (trtsdt, trtedt) = (sl.dates("TRTSDT")?, sl.dates("TRTEDT")?);
    let from_sl_char = |name: &str| -> Result<Values, DeriveError> {
        let v = sl.chars(name)?;
        Ok(Values::Char(
            row.iter().map(|r| r.and_then(|r| v[r].clone())).collect(),
        ))
    };
    let from_sl_num = |name: &str| -> Result<Values, DeriveError> {
        let v = sl.nums(name)?;
        Ok(Values::Num(
            row.iter().map(|r| r.and_then(|r| v[r])).collect(),
        ))
    };
    let copy =
        |name: &str| -> Result<Values, DeriveError> { Ok(Values::Char(ae.chars(name)?.to_vec())) };

    let astdt: Vec<_> = ae.chars("AESTDTC")?.iter().map(dt).collect();
    let aendt: Vec<_> = ae.chars("AEENDTC")?.iter().map(dt).collect();
    let sdt: Vec<_> = row.iter().map(|r| r.and_then(|r| trtsdt[r])).collect();
    let edt: Vec<_> = row.iter().map(|r| r.and_then(|r| trtedt[r])).collect();
    let n = ae.n();

    let trtem: Vec<bool> = (0..n)
        .map(|i| match (astdt[i], sdt[i], edt[i]) {
            (Some(a), Some(s), Some(e)) => a >= s && a <= e.add_days(TEAE_WINDOW_DAYS),
            _ => false,
        })
        .collect();

    // First occurrences among treatment-emergent events, in start order.
    let (decod, bodsys) = (ae.chars("AEDECOD")?, ae.chars("AEBODSYS")?);
    let aeseq = ae.ints("AESEQ")?;
    let mut order: Vec<usize> = (0..n).filter(|&i| trtem[i]).collect();
    order.sort_by(|&a, &b| {
        (&usubjid[a], astdt[a], aeseq[a]).cmp(&(&usubjid[b], astdt[b], aeseq[b]))
    });
    let (mut occ, mut occ_soc, mut occ_pt) = (vec![false; n], vec![false; n], vec![false; n]);
    let (mut seen, mut seen_soc, mut seen_pt) = (HashSet::new(), HashSet::new(), HashSet::new());
    for i in order {
        occ[i] = seen.insert(&usubjid[i]);
        occ_soc[i] = seen_soc.insert((&usubjid[i], &bodsys[i]));
        occ_pt[i] = seen_pt.insert((&usubjid[i], &decod[i]));
    }
    let y = |b: &bool| b.then(|| "Y".to_string());
    let sev = ae.chars("AESEV")?;

    let mut cols: Vec<(&str, Values)> = vec![
        ("STUDYID", copy("STUDYID")?),
        ("USUBJID", copy("USUBJID")?),
        ("SUBJID", from_sl_char("SUBJID")?),
        ("SITEID", from_sl_char("SITEID")?),
        ("AGE", from_sl_num("AGE")?),
        ("SEX", from_sl_char("SEX")?),
        ("RACE", from_sl_char("RACE")?),
        ("SAFFL", from_sl_char("SAFFL")?),
        ("TRTA", from_sl_char("TRT01A")?),
        ("TRTAN", from_sl_num("TRT01AN")?),
        ("TRTSDT", Values::Date(sdt.clone())),
        ("TRTEDT", Values::Date(edt)),
        ("AESEQ", Values::Int(aeseq.clone())),
        ("AETERM", copy("AETERM")?),
        ("AEDECOD", copy("AEDECOD")?),
        ("AEBODSYS", copy("AEBODSYS")?),
        ("AESOC", copy("AESOC")?),
        ("AESEV", copy("AESEV")?),
        ("AESER", copy("AESER")?),
        ("AEREL", copy("AEREL")?),
        ("AEACN", copy("AEACN")?),
        ("AEOUT", copy("AEOUT")?),
        ("ASEV", copy("AESEV")?),
        (
            "ASEVN",
            Values::Num(
                sev.iter()
                    .map(|s| match s.as_deref() {
                        Some("MILD") => Some(1.0),
                        Some("MODERATE") => Some(2.0),
                        Some("SEVERE") => Some(3.0),
                        _ => None,
                    })
                    .collect(),
            ),
        ),
        ("AREL", copy("AEREL")?),
        ("ASTDT", Values::Date(astdt.clone())),
        ("AENDT", Values::Date(aendt.clone())),
        (
            "ASTDY",
            Values::Int((0..n).map(|i| rel_day(astdt[i], sdt[i])).collect()),
        ),
        (
            "AENDY",
            Values::Int((0..n).map(|i| rel_day(aendt[i], sdt[i])).collect()),
        ),
        (
            "ADURN",
            Values::Num(
                (0..n)
                    .map(|i| match (astdt[i], aendt[i]) {
                        (Some(a), Some(b)) => Some(f64::from(b.days_since(a) + 1)),
                        _ => None,
                    })
                    .collect(),
            ),
        ),
        (
            "ADURU",
            Values::Char(
                (0..n)
                    .map(|i| aendt[i].and(astdt[i]).map(|_| "DAYS".to_string()))
                    .collect(),
            ),
        ),
        ("TRTEMFL", Values::Char(trtem.iter().map(y).collect())),
        ("AOCCFL", Values::Char(occ.iter().map(y).collect())),
        ("AOCCSFL", Values::Char(occ_soc.iter().map(y).collect())),
        ("AOCCPFL", Values::Char(occ_pt.iter().map(y).collect())),
    ];
    // DLTFL (admiral derive_vars_merged() on SUPPAE AEDLT, as metatools::combine_supp()).
    if let Some(sup) = suppae.map(Frame::new) {
        let (su, qnam, idvarval, qval) = (
            sup.chars("USUBJID")?,
            sup.chars("QNAM")?,
            sup.chars("IDVARVAL")?,
            sup.chars("QVAL")?,
        );
        let mut dlt: std::collections::HashMap<(&str, i64), &str> =
            std::collections::HashMap::new();
        for i in 0..sup.n() {
            if let (Some(x), Some("AEDLT"), Some(v), Some(q)) = (
                su[i].as_deref(),
                qnam[i].as_deref(),
                idvarval[i].as_deref(),
                qval[i].as_deref(),
            ) {
                if let Ok(seq) = v.parse::<i64>() {
                    dlt.insert((x, seq), q);
                }
            }
        }
        if !dlt.is_empty() {
            cols.push((
                "DLTFL",
                Values::Char(
                    (0..n)
                        .map(|i| {
                            let x = usubjid[i].as_deref()?;
                            dlt.get(&(x, aeseq[i]?)).map(|q| q.to_string())
                        })
                        .collect(),
                ),
            ));
        }
    }
    let periods = crate::Periods::new(adsl)?;
    let keys: Vec<(usize, Option<trialforge_core::Date>)> = row
        .iter()
        .zip(&astdt)
        .map(|(r, d)| (r.unwrap_or(0), *d))
        .collect();
    crate::by_period(&mut cols, &periods, &keys);
    Ok(build("ADAE", cols)?)
}
