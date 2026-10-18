//! L1 `L1.LUGANO`: for lymphoma (RS `RSCAT = LUGANO CLASSIFICATION`), every product of
//! perpendicular diameters, sum of products and CT- and PET-CT-based overall response is
//! recomputed from the TR measurements, lesion states and Deauville scores.

use std::collections::{BTreeMap, HashMap, HashSet};

use trialforge_core::Dataset;

use crate::table::get;
use crate::{Level, Report};

const L: Level = Level::L1;
const RULE: &str = "L1.LUGANO";

/// A reported value: its row and value.
type Reported = Option<(usize, Option<f64>)>;

#[derive(Default)]
struct Visit {
    /// Per target lesion: (longest, perpendicular, PPD row and value).
    target: BTreeMap<String, (Option<f64>, Option<f64>, Reported)>,
    spd: Reported,
    deauville: Option<f64>,
    non_target: Vec<String>,
    new_lesion: bool,
}

/// Whether the study's RS follows Lugano.
pub(crate) fn applies(sdtm: &[Dataset]) -> bool {
    get(sdtm, "RS")
        .and_then(|rs| {
            rs.c("RSCAT").map(|c| {
                c.iter()
                    .any(|x| x.as_deref() == Some("LUGANO CLASSIFICATION"))
            })
        })
        .unwrap_or(false)
}

pub(crate) fn l1(sdtm: &[Dataset], r: &mut Report) {
    let (Some(tr), Some(rs), Some(tu)) = (get(sdtm, "TR"), get(sdtm, "RS"), get(sdtm, "TU")) else {
        return;
    };
    let mut nodes: HashSet<(&str, &str)> = HashSet::new();
    if let (Some(u), Some(l), Some(loc)) = (tu.c("USUBJID"), tu.c("TULNKID"), tu.c("TULOC")) {
        for i in 0..tu.n() {
            if let (Some(x), Some(id), Some("LYMPH NODE")) =
                (u[i].as_deref(), l[i].as_deref(), loc[i].as_deref())
            {
                nodes.insert((x, id));
            }
        }
    }
    let (Some(u), Some(grp), Some(test), Some(lnk), Some(orres), Some(val), Some(vn)) = (
        tr.c("USUBJID"),
        tr.c("TRGRPID"),
        tr.c("TRTESTCD"),
        tr.c("TRLNKID"),
        tr.c("TRORRES"),
        tr.f("TRSTRESN"),
        tr.i("VISITNUM"),
    ) else {
        return;
    };
    let mut subjects: HashMap<&str, BTreeMap<i64, Visit>> = HashMap::new();
    for i in 0..tr.n() {
        let (Some(x), Some(v)) = (u[i].as_deref(), vn[i]) else {
            continue;
        };
        let visit = subjects.entry(x).or_default().entry(v).or_default();
        let id = lnk[i].clone().unwrap_or_default();
        match (grp[i].as_deref(), test[i].as_deref()) {
            (Some("TARGET"), Some("LDIAM")) => visit.target.entry(id).or_default().0 = val[i],
            (Some("TARGET"), Some("LPERP")) => visit.target.entry(id).or_default().1 = val[i],
            (Some("TARGET"), Some("PPD")) => {
                visit.target.entry(id).or_default().2 = Some((i, val[i]))
            }
            (_, Some("SUMPPD")) => visit.spd = Some((i, val[i])),
            (_, Some("FDP5PS")) => visit.deauville = val[i],
            (Some("NON-TARGET"), Some("TUMSTATE")) => {
                visit.non_target.push(orres[i].clone().unwrap_or_default())
            }
            (Some("NEW"), Some("TUMSTATE")) => visit.new_lesion = true,
            _ => {}
        }
    }
    let mut responses: HashMap<(&str, i64, &str), (usize, &str)> = HashMap::new();
    if let (Some(u), Some(vn), Some(sc), Some(res)) = (
        rs.c("USUBJID"),
        rs.i("VISITNUM"),
        rs.c("RSSCAT"),
        rs.c("RSSTRESC"),
    ) {
        for i in 0..rs.n() {
            if let (Some(x), Some(v), Some(s), Some(o)) =
                (u[i].as_deref(), vn[i], sc[i].as_deref(), res[i].as_deref())
            {
                responses.insert((x, v, s), (i, o));
            }
        }
    }
    let mut matched = 0;
    for (x, visits) in &subjects {
        let mut baseline: Option<f64> = None;
        let mut nadir: HashMap<String, (f64, f64)> = HashMap::new();
        for (&vn, v) in visits {
            let products: Vec<Option<f64>> = v
                .target
                .values()
                .map(|t| t.0.zip(t.1).map(|(l, w)| l * w))
                .collect();
            for (t, p) in v.target.values().zip(&products) {
                if let Some((row, got)) = t.2
                    && got != *p
                {
                    r.add(
                        L,
                        RULE,
                        "TR",
                        Some("TRSTRESN"),
                        Some(row),
                        "PPD is not the product of the diameters".into(),
                    );
                }
            }
            let complete = products.iter().all(Option::is_some);
            let spd: f64 = products.iter().flatten().sum();
            match v.spd {
                Some((row, got)) if !complete || got.is_none_or(|g| (g - spd).abs() > 1e-6) => {
                    r.add(
                        L,
                        RULE,
                        "TR",
                        Some("TRSTRESN"),
                        Some(row),
                        "SUMPPD is not the sum of the products".into(),
                    );
                }
                None if complete && !v.target.is_empty() => {
                    r.add(
                        L,
                        RULE,
                        "TR",
                        Some("TRTESTCD"),
                        None,
                        format!("{x} visit {vn}: no SUMPPD although every lesion was measured"),
                    );
                }
                _ => {}
            }
            let Some(base) = baseline else {
                baseline = Some(spd);
                for (id, t) in &v.target {
                    if let (Some(l), Some(w)) = (t.0, t.1) {
                        nadir.insert(id.clone(), (l * w, l));
                    }
                }
                continue;
            };
            // CT-based response.
            let lesion_pd = v
                .target
                .iter()
                .any(|(id, t)| match (t.0, t.1, nadir.get(id)) {
                    (Some(l), Some(w), Some(&(np, nl))) => {
                        l > 15.0 && l * w >= 1.5 * np && l - nl >= 5.0
                    }
                    _ => false,
                });
            let ntl_pd = v.non_target.iter().any(|s| s == "UNEQUIVOCAL PROGRESSION");
            let all_small = v.target.iter().all(|(id, t)| {
                let l = t.0.unwrap_or(f64::INFINITY);
                if nodes.contains(&(*x, id.as_str())) {
                    l <= 15.0
                } else {
                    l == 0.0
                }
            });
            let ct = if v.new_lesion || ntl_pd || lesion_pd {
                "PD"
            } else if !complete {
                "NE"
            } else if all_small && v.non_target.iter().all(|s| s == "ABSENT") {
                "CR"
            } else if spd <= 0.5 * base {
                "PR"
            } else {
                "SD"
            };
            let pet = match (ct, v.deauville) {
                ("PD", _) => "PMD",
                (_, None) => "NE",
                (_, Some(d)) if d <= 3.0 => "CMR",
                ("NE", _) => "NE",
                _ if spd <= 0.5 * base => "PMR",
                _ => "SMD",
            };
            for (t, l) in v.target.iter() {
                if let (Some(len), Some(w)) = (l.0, l.1) {
                    let e = nadir.entry(t.clone()).or_insert((len * w, len));
                    e.0 = e.0.min(len * w);
                    e.1 = e.1.min(len);
                }
            }
            for (scat, want) in [("CT", ct), ("PET-CT", pet)] {
                let got = responses.get(&(*x, vn, scat));
                matched += usize::from(got.is_some());
                if got.map(|g| g.1) != Some(want) {
                    r.add(
                        L,
                        RULE,
                        "RS",
                        Some("RSSTRESC"),
                        got.map(|g| g.0),
                        format!(
                            "{x} visit {vn} {scat}: Lugano gives {want}, RS has {:?}",
                            got.map(|g| g.1)
                        ),
                    );
                }
            }
        }
    }
    if matched != responses.len() {
        r.add(
            L,
            RULE,
            "RS",
            None,
            None,
            format!(
                "{} RS records have no TR assessment",
                responses.len() - matched
            ),
        );
    }
}
