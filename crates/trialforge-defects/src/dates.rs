//! Date defects.

use std::collections::HashMap;

use trialforge_core::{Date, Draws, Stream, study_day};

use crate::{Amount, Cell, Ctx, Defect, Edit, Expected, chars, ints, parse_date, pick};

/// The largest step back from AESTDTC, in days.
const MAX_DAYS_BACK: i32 = 14;

/// `ae_end_before_start`: AEENDTC moves 1–14 days before AESTDTC, never before informed
/// consent (so only the order is wrong), and AEENDY is recomputed from it. ADAE keeps the
/// clean dates, so its ASTDT/AENDT (and ADURN) no longer match AE.
pub(crate) fn ae_end_before_start(cx: &mut Ctx<'_>, stream: Stream, amount: Amount) -> Vec<Defect> {
    let (Some(ae), Some(dm)) = (cx.sdtm("AE"), cx.sdtm("DM")) else {
        return Vec::new();
    };
    // USUBJID → (RFICDTC, RFSTDTC).
    let refs: HashMap<&str, (Option<Date>, Option<Date>)> = match (
        chars(dm, "USUBJID"),
        chars(dm, "RFICDTC"),
        chars(dm, "RFSTDTC"),
    ) {
        (Some(u), Some(ic), Some(st)) => (0..dm.n_rows())
            .filter_map(|i| {
                let date = |v: &Option<String>| v.as_deref().and_then(parse_date);
                Some((u[i].as_deref()?, (date(&ic[i]), date(&st[i]))))
            })
            .collect(),
        _ => return Vec::new(),
    };
    let (Some(u), Some(seq), Some(st), Some(en)) = (
        chars(ae, "USUBJID"),
        ints(ae, "AESEQ"),
        chars(ae, "AESTDTC"),
        chars(ae, "AEENDTC"),
    ) else {
        return Vec::new();
    };
    let has_endy = ae.column("AEENDY").is_some();

    // (row, start, days from consent to start).
    let candidates: Vec<(usize, Date, i32)> = (0..ae.n_rows())
        .filter(|&i| !cx.is_used("AE", i))
        .filter_map(|i| {
            let start = parse_date(st[i].as_deref()?)?;
            parse_date(en[i].as_deref()?)?;
            let consent = refs.get(u[i].as_deref()?)?.0?;
            let slack = start.days_since(consent);
            (slack >= 1).then_some((i, start, slack))
        })
        .collect();
    let chosen = pick(
        &mut Draws::new(stream.label("select")),
        candidates.len(),
        amount,
        |_| true,
    );

    let adae = cx.adam("ADAE");
    let adae_row = |x: &str, s: i64| -> Option<usize> {
        let ad = adae?;
        let (au, aseq) = (chars(ad, "USUBJID")?, ints(ad, "AESEQ")?);
        (0..ad.n_rows()).find(|&j| au[j].as_deref() == Some(x) && aseq[j] == Some(s))
    };

    let mut edits = Vec::new();
    let mut defects = Vec::new();
    for c in chosen {
        let (i, start, slack) = candidates[c];
        let mut d = Draws::new(stream.label("value").index(i as u64));
        let back = d.int_range(1, i64::from(slack.min(MAX_DAYS_BACK))) as i32;
        let end = start.add_days(-back);
        let original = en[i].clone().expect("a candidate has an end date");
        // Keep any time of day.
        let injected = format!("{end}{}", &original[10..]);
        let x = u[i].clone().expect("a candidate has a subject");
        edits.push(Edit {
            variable: "AEENDTC".into(),
            row: i,
            value: Cell::Char(Some(injected.clone())),
        });
        if has_endy && let Some(rfst) = refs[x.as_str()].1 {
            edits.push(Edit {
                variable: "AEENDY".into(),
                row: i,
                value: Cell::Int(Some(i64::from(study_day(end, rfst)))),
            });
        }
        let mut expected = vec![Expected {
            rule: "L1.ORDER",
            dataset: "AE".into(),
            row: Some(i),
        }];
        if let Some(s) = seq[i]
            && let Some(j) = adae_row(&x, s)
        {
            expected.push(Expected {
                rule: "L2.ADAE",
                dataset: "ADAE".into(),
                row: Some(j),
            });
        }
        defects.push(Defect {
            id: "ae_end_before_start",
            dataset: "AE".into(),
            row: i,
            key: seq[i].map(|s| format!("AESEQ={s}")),
            usubjid: Some(x),
            variable: "AEENDTC".into(),
            original: Some(original),
            injected: Some(injected),
            paired_row: None,
            expected,
        });
    }
    for d in &defects {
        cx.used.insert((d.dataset.clone(), d.row));
    }
    cx.edit_sdtm("AE", &edits);
    defects
}
