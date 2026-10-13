//! Define-XML 2.1 writer.
//!
//! One `define.xml` per standard (SDTM or ADaM), built from the metadata bundle (labels, core,
//! origins, methods, codelists) and from the data itself (lengths, data types, the codelist terms
//! actually used, and value-level metadata per test or parameter). Output is deterministic: the
//! creation timestamp is the study date passed in, not the wall clock.
//!
//! Contents, in schema order:
//! - `def:Standards`: the IG and the CDISC CT package.
//! - `def:ValueListDef` / `def:WhereClauseDef`: value-level metadata for `--ORRES` and `--STRESN`
//!   by `--TESTCD` (LB, VS) and `AVAL` by `PARAMCD` (BDS).
//! - `ItemGroupDef` per dataset (with `def:Class` and `def:leaf` pointing at the XPT file),
//!   `ItemDef` per variable (with `def:Origin`), `CodeList` per controlled codelist (only terms
//!   that occur) plus a sponsor `PARAMCD` codelist with decodes, `MethodDef` per derivation.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use trialforge_cdisc::{Core, DatasetSpec, VarType, ct, spec};
use trialforge_core::{Dataset, Date, Values};

/// Which standard a define.xml describes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Standard {
    /// SDTM tabulation datasets.
    Sdtm,
    /// ADaM analysis datasets.
    Adam,
}

/// Study-level information for the define.xml header.
#[derive(Clone, Debug, PartialEq)]
pub struct DefineInfo<'a> {
    /// `STUDYID`.
    pub studyid: &'a str,
    /// Study description (e.g. the trial title).
    pub description: &'a str,
    /// Creation date (deterministic; trialforge passes the study's data date).
    pub created: Date,
    /// Which standard.
    pub standard: Standard,
}

/// Errors from building a define.xml.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DefineError(pub String);

impl core::fmt::Display for DefineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for DefineError {}

fn esc(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            c => out.push(c),
        }
    }
    out
}

/// OID-safe text: letters, digits, `.`, `_` and `-` only.
fn oid(s: &str) -> String {
    s.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '.' | '_' | '-') {
                c
            } else {
                '_'
            }
        })
        .collect()
}

fn description(out: &mut String, indent: &str, text: &str) {
    let _ = writeln!(
        out,
        "{indent}<Description><TranslatedText xml:lang=\"en\">{}</TranslatedText></Description>",
        esc(text)
    );
}

/// The ODM data type, length and significant digits of a set of values.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Shape {
    data_type: &'static str,
    length: Option<usize>,
    digits: Option<usize>,
    display_format: Option<&'static str>,
}

fn decimals(x: f64) -> usize {
    let s = x.to_string();
    s.split_once('.').map_or(0, |(_, f)| f.len())
}

fn shape_of<'a>(var: &str, values: impl Iterator<Item = ValueRef<'a>>) -> Shape {
    let mut kind = None;
    let (mut len, mut digits, mut partial, mut all_ints) = (0usize, 0usize, false, true);
    for v in values {
        match v {
            ValueRef::Char(s) => {
                kind = Some("char");
                len = len.max(s.len());
                if s.len() < 10 {
                    partial = true;
                }
            }
            ValueRef::Int(n) => {
                kind = Some("int");
                len = len.max(n.to_string().len());
            }
            ValueRef::Num(x) => {
                kind = Some("num");
                if x.fract() != 0.0 {
                    all_ints = false;
                }
                len = len.max(x.to_string().len());
                digits = digits.max(decimals(x));
            }
            ValueRef::Date => kind = Some("date"),
        }
    }
    let is_dtc = var.ends_with("DTC");
    match kind {
        Some("char") if is_dtc => Shape {
            data_type: if partial { "partialDate" } else { "date" },
            length: None,
            digits: None,
            display_format: None,
        },
        Some("int") => Shape {
            data_type: "integer",
            length: Some(len.max(1)),
            digits: None,
            display_format: None,
        },
        Some("num") if all_ints => Shape {
            data_type: "integer",
            length: Some(len.max(1)),
            digits: None,
            display_format: None,
        },
        Some("num") => Shape {
            data_type: "float",
            length: Some(len.max(1)),
            digits: Some(digits),
            display_format: None,
        },
        Some("date") => Shape {
            data_type: "integer",
            length: Some(8),
            digits: None,
            display_format: Some("DATE9."),
        },
        _ if is_dtc => Shape {
            data_type: "date",
            length: None,
            digits: None,
            display_format: None,
        },
        _ => Shape {
            data_type: "text",
            length: Some(len.max(1)),
            digits: None,
            display_format: None,
        },
    }
}

#[derive(Clone, Copy)]
enum ValueRef<'a> {
    Char(&'a str),
    Int(i64),
    Num(f64),
    Date,
}

fn values_at<'a>(
    v: &'a Values,
    rows: impl Iterator<Item = usize> + 'a,
) -> impl Iterator<Item = ValueRef<'a>> + 'a {
    rows.filter_map(move |i| match v {
        Values::Char(x) => x[i].as_deref().map(ValueRef::Char),
        Values::Int(x) => x[i].map(ValueRef::Int),
        Values::Num(x) => x[i].map(ValueRef::Num),
        Values::Date(x) => x[i].map(|_| ValueRef::Date),
    })
}

fn shape_attrs(s: &Shape) -> String {
    let mut a = format!(" DataType=\"{}\"", s.data_type);
    if let Some(l) = s.length {
        let _ = write!(a, " Length=\"{l}\"");
    }
    if let Some(d) = s.digits {
        let _ = write!(a, " SignificantDigits=\"{d}\"");
    }
    if let Some(f) = s.display_format {
        let _ = write!(a, " def:DisplayFormat=\"{f}\"");
    }
    a
}

fn origin(out: &mut String, indent: &str, origin: &str, source: &str) {
    match origin {
        "Predecessor" => {
            let _ = writeln!(out, "{indent}<def:Origin Type=\"Predecessor\">");
            description(out, &format!("{indent}  "), source);
            let _ = writeln!(out, "{indent}</def:Origin>");
        }
        "Collected" | "Assigned" if !source.is_empty() => {
            let _ = writeln!(
                out,
                "{indent}<def:Origin Type=\"{origin}\" Source=\"{}\"/>",
                esc(source)
            );
        }
        other => {
            let _ = writeln!(out, "{indent}<def:Origin Type=\"{}\"/>", esc(other));
        }
    }
}

/// Value-level metadata for one variable: (`whereclause variable`, values in order).
struct ValueList<'a> {
    variable: &'a str,
    by: &'a str,
    groups: Vec<(String, Vec<usize>)>,
}

fn value_lists<'a>(ds: &'a Dataset, standard: Standard) -> Vec<ValueList<'a>> {
    let (by, targets): (&str, Vec<&str>) = match (standard, ds.name.as_str()) {
        (Standard::Sdtm, "LB") => ("LBTESTCD", vec!["LBORRES", "LBSTRESN"]),
        (Standard::Sdtm, "VS") => ("VSTESTCD", vec!["VSORRES", "VSSTRESN"]),
        (Standard::Adam, _) if ds.column("PARAMCD").is_some() => ("PARAMCD", vec!["AVAL"]),
        _ => return Vec::new(),
    };
    let Some(keys) = ds.column(by).and_then(|c| c.as_char()) else {
        return Vec::new();
    };
    let mut order: Vec<String> = Vec::new();
    let mut rows: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (i, k) in keys.iter().enumerate() {
        let Some(k) = k else { continue };
        if !rows.contains_key(k) {
            order.push(k.clone());
        }
        rows.entry(k.clone()).or_default().push(i);
    }
    targets
        .into_iter()
        .filter(|t| ds.column(t).is_some())
        .map(|t| ValueList {
            variable: t,
            by,
            groups: order.iter().map(|k| (k.clone(), rows[k].clone())).collect(),
        })
        .collect()
}

/// Builds the define.xml for `datasets` (all of one standard).
pub fn to_string(datasets: &[Dataset], info: &DefineInfo<'_>) -> Result<String, DefineError> {
    let (ig_name, ig_version, purpose, label) = match info.standard {
        Standard::Sdtm => ("SDTMIG", "3.4", "Tabulation", "SDTM"),
        Standard::Adam => ("ADaMIG", "1.3", "Analysis", "ADaM"),
    };
    let specs: Vec<(&Dataset, &DatasetSpec)> = datasets
        .iter()
        .map(|d| {
            spec(&d.name)
                .map(|s| (d, s))
                .ok_or_else(|| DefineError(format!("{} is not in the metadata bundle", d.name)))
        })
        .collect::<Result<_, _>>()?;
    let ct = ct();
    let study = oid(info.studyid);
    let mut out = String::new();

    let _ = writeln!(out, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
    let _ = writeln!(
        out,
        "<ODM xmlns=\"http://www.cdisc.org/ns/odm/v1.3\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" \
         xmlns:def=\"http://www.cdisc.org/ns/def/v2.1\" ODMVersion=\"1.3.2\" FileType=\"Snapshot\" \
         FileOID=\"DEF.{study}.{label}\" CreationDateTime=\"{}T00:00:00\" Originator=\"trialforge\" \
         SourceSystem=\"trialforge\" def:Context=\"Submission\">",
        info.created
    );
    let _ = writeln!(out, "  <Study OID=\"ST.{study}\">");
    let _ = writeln!(out, "    <GlobalVariables>");
    let _ = writeln!(out, "      <StudyName>{}</StudyName>", esc(info.studyid));
    let _ = writeln!(
        out,
        "      <StudyDescription>{}</StudyDescription>",
        esc(info.description)
    );
    let _ = writeln!(
        out,
        "      <ProtocolName>{}</ProtocolName>",
        esc(info.studyid)
    );
    let _ = writeln!(out, "    </GlobalVariables>");
    let _ = writeln!(
        out,
        "    <MetaDataVersion OID=\"MDV.{study}.{label}\" Name=\"{} {label}\" Description=\"{label} metadata for {}\" def:DefineVersion=\"2.1.0\">",
        esc(info.studyid),
        esc(info.studyid)
    );
    let _ = writeln!(out, "      <def:Standards>");
    let _ = writeln!(
        out,
        "        <def:Standard OID=\"STD.IG\" Name=\"{ig_name}\" Type=\"IG\" Version=\"{ig_version}\" Status=\"Final\"/>"
    );
    let _ = writeln!(
        out,
        "        <def:Standard OID=\"STD.CT\" Name=\"CDISC/NCI\" Type=\"CT\" PublishingSet=\"SDTM\" Version=\"{}\" Status=\"Final\"/>",
        ct.version
    );
    let _ = writeln!(out, "      </def:Standards>");

    // Value-level metadata.
    let mut where_clauses: BTreeMap<String, (String, String)> = BTreeMap::new(); // oid → (item oid, value)
    let mut vlm_items: Vec<(String, String, Shape, String, String, String)> = Vec::new(); // (oid, name, shape, desc, origin, source)
    let mut vlm_refs: BTreeMap<(String, String), String> = BTreeMap::new(); // (dataset, var) → VL oid
    for (d, s) in &specs {
        for vl in value_lists(d, info.standard) {
            let col = d.column(vl.variable).expect("filtered");
            let var_spec = s.variable(vl.variable).expect("bundle variable");
            let vl_oid = format!("VL.{}.{}", d.name, vl.variable);
            let _ = writeln!(out, "      <def:ValueListDef OID=\"{vl_oid}\">");
            for (k, (value, rows)) in vl.groups.iter().enumerate() {
                let item_oid = format!("IT.{}.{}.{}", d.name, vl.variable, oid(value));
                let wc_oid = format!("WC.{}.{}.{}", d.name, vl.by, oid(value));
                let _ = writeln!(
                    out,
                    "        <ItemRef ItemOID=\"{item_oid}\" OrderNumber=\"{}\" Mandatory=\"No\">",
                    k + 1
                );
                let _ = writeln!(
                    out,
                    "          <def:WhereClauseRef WhereClauseOID=\"{wc_oid}\"/>"
                );
                let _ = writeln!(out, "        </ItemRef>");
                where_clauses.insert(wc_oid, (format!("IT.{}.{}", d.name, vl.by), value.clone()));
                let shape = shape_of(vl.variable, values_at(&col.values, rows.iter().copied()));
                vlm_items.push((
                    item_oid,
                    vl.variable.to_string(),
                    shape,
                    format!("{} where {} = {value}", var_spec.label, vl.by),
                    var_spec.origin.clone(),
                    var_spec.source.clone(),
                ));
            }
            let _ = writeln!(out, "      </def:ValueListDef>");
            vlm_refs.insert((d.name.clone(), vl.variable.to_string()), vl_oid);
        }
    }
    for (wc, (item, value)) in &where_clauses {
        let _ = writeln!(out, "      <def:WhereClauseDef OID=\"{wc}\">");
        let _ = writeln!(
            out,
            "        <RangeCheck SoftHard=\"Soft\" def:ItemOID=\"{item}\" Comparator=\"EQ\">"
        );
        let _ = writeln!(out, "          <CheckValue>{}</CheckValue>", esc(value));
        let _ = writeln!(out, "        </RangeCheck>");
        let _ = writeln!(out, "      </def:WhereClauseDef>");
    }

    // Item groups.
    let mut methods: Vec<(String, String, String)> = Vec::new(); // (oid, name, text)
    let mut codelist_use: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    let mut paramcd_lists: Vec<(String, Vec<(String, String)>)> = Vec::new();
    for (d, s) in &specs {
        let reference = s.class == "TRIAL DESIGN";
        let domain = if info.standard == Standard::Sdtm && !d.name.starts_with("SUPP") {
            format!(" Domain=\"{}\"", d.name)
        } else {
            String::new()
        };
        let _ = writeln!(
            out,
            "      <ItemGroupDef OID=\"IG.{n}\"{domain} Name=\"{n}\" SASDatasetName=\"{n}\" Repeating=\"{rep}\" IsReferenceData=\"{refd}\" Purpose=\"{purpose}\" def:Structure=\"{st}\" def:StandardOID=\"STD.IG\" def:ArchiveLocationID=\"LF.{n}\">",
            n = d.name,
            rep = if s.structure == "One record per subject" {
                "No"
            } else {
                "Yes"
            },
            refd = if reference { "Yes" } else { "No" },
            st = esc(&s.structure),
        );
        description(&mut out, "        ", &d.label);
        for (k, c) in d.columns().iter().enumerate() {
            let v = s.variable(&c.name).expect("built from the bundle");
            let key = s
                .keys
                .iter()
                .position(|x| *x == c.name)
                .map(|p| format!(" KeySequence=\"{}\"", p + 1))
                .unwrap_or_default();
            let mandatory = if v.core == Core::Req { "Yes" } else { "No" };
            let method = if v.origin == "Derived" {
                let m = format!("MT.{}.{}", d.name, c.name);
                methods.push((
                    m.clone(),
                    format!("Algorithm to derive {}.{}", d.name, c.name),
                    v.method.clone(),
                ));
                format!(" MethodOID=\"{m}\"")
            } else {
                String::new()
            };
            let _ = writeln!(
                out,
                "        <ItemRef ItemOID=\"IT.{}.{}\" OrderNumber=\"{}\" Mandatory=\"{mandatory}\"{key}{method}/>",
                d.name,
                c.name,
                k + 1
            );
            if let (Some(cl), Values::Char(vals)) = (&v.codelist, &c.values) {
                codelist_use
                    .entry(cl.clone())
                    .or_default()
                    .extend(vals.iter().flatten().cloned());
            }
        }
        let _ = writeln!(out, "        <def:Class Name=\"{}\"/>", esc(&s.class));
        let file = format!("{}.xpt", d.name.to_ascii_lowercase());
        let _ = writeln!(
            out,
            "        <def:leaf ID=\"LF.{}\" xlink:href=\"{file}\">",
            d.name
        );
        let _ = writeln!(out, "          <def:title>{file}</def:title>");
        let _ = writeln!(out, "        </def:leaf>");
        let _ = writeln!(out, "      </ItemGroupDef>");

        // A sponsor PARAMCD codelist with decodes, for BDS datasets.
        if let (Some(pc), Some(param)) = (
            d.column("PARAMCD").and_then(|c| c.as_char()),
            d.column("PARAM").and_then(|c| c.as_char()),
        ) {
            let mut items: Vec<(String, String)> = Vec::new();
            for (a, b) in pc.iter().zip(param) {
                if let (Some(a), Some(b)) = (a, b) {
                    if !items.iter().any(|(x, _)| x == a) {
                        items.push((a.clone(), b.clone()));
                    }
                }
            }
            paramcd_lists.push((d.name.clone(), items));
        }
    }

    // Item definitions.
    for (d, s) in &specs {
        for c in d.columns() {
            let v = s.variable(&c.name).expect("built from the bundle");
            let mut shape = shape_of(&c.name, values_at(&c.values, 0..d.n_rows()));
            if v.var_type == VarType::Char && c.name.ends_with("DTC") && d.n_rows() == 0 {
                shape.data_type = "date";
            }
            let _ = writeln!(
                out,
                "      <ItemDef OID=\"IT.{}.{}\" Name=\"{}\" SASFieldName=\"{}\"{}>",
                d.name,
                c.name,
                c.name,
                c.name,
                shape_attrs(&shape)
            );
            description(&mut out, "        ", &v.label);
            if let Some(cl) = &v.codelist {
                let _ = writeln!(out, "        <CodeListRef CodeListOID=\"CL.{}\"/>", oid(cl));
            } else if c.name == "PARAMCD" {
                let _ = writeln!(
                    out,
                    "        <CodeListRef CodeListOID=\"CL.{}.PARAMCD\"/>",
                    d.name
                );
            }
            origin(&mut out, "        ", &v.origin, &v.source);
            if let Some(vl) = vlm_refs.get(&(d.name.clone(), c.name.clone())) {
                let _ = writeln!(out, "        <def:ValueListRef ValueListOID=\"{vl}\"/>");
            }
            let _ = writeln!(out, "      </ItemDef>");
        }
    }
    for (item_oid, name, shape, desc, orig, source) in &vlm_items {
        let _ = writeln!(
            out,
            "      <ItemDef OID=\"{item_oid}\" Name=\"{name}\" SASFieldName=\"{name}\"{}>",
            shape_attrs(shape)
        );
        description(&mut out, "        ", desc);
        origin(&mut out, "        ", orig, source);
        let _ = writeln!(out, "      </ItemDef>");
    }

    // Code lists: CDISC CT (only the terms used, in CT order), then sponsor PARAMCD lists.
    for (name, used) in &codelist_use {
        let Some(list) = ct.codelist(name) else {
            continue;
        };
        let _ = writeln!(
            out,
            "      <CodeList OID=\"CL.{}\" Name=\"{}\" DataType=\"text\" def:StandardOID=\"STD.CT\">",
            oid(name),
            esc(name)
        );
        for (k, term) in list.terms.iter().filter(|t| used.contains(*t)).enumerate() {
            let _ = writeln!(
                out,
                "        <EnumeratedItem CodedValue=\"{}\" OrderNumber=\"{}\"/>",
                esc(term),
                k + 1
            );
        }
        let _ = writeln!(out, "      </CodeList>");
    }
    for (ds, items) in &paramcd_lists {
        let _ = writeln!(
            out,
            "      <CodeList OID=\"CL.{ds}.PARAMCD\" Name=\"{ds} Parameter Code\" DataType=\"text\">"
        );
        for (k, (code, decode)) in items.iter().enumerate() {
            let _ = writeln!(
                out,
                "        <CodeListItem CodedValue=\"{}\" OrderNumber=\"{}\">",
                esc(code),
                k + 1
            );
            let _ = writeln!(
                out,
                "          <Decode><TranslatedText xml:lang=\"en\">{}</TranslatedText></Decode>",
                esc(decode)
            );
            let _ = writeln!(out, "        </CodeListItem>");
        }
        let _ = writeln!(out, "      </CodeList>");
    }

    for (m, name, text) in &methods {
        let _ = writeln!(
            out,
            "      <MethodDef OID=\"{m}\" Name=\"{}\" Type=\"Computation\">",
            esc(name)
        );
        description(&mut out, "        ", text);
        let _ = writeln!(out, "      </MethodDef>");
    }
    let _ = writeln!(out, "    </MetaDataVersion>");
    let _ = writeln!(out, "  </Study>");
    let _ = writeln!(out, "</ODM>");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escapes_and_oids() {
        assert_eq!(esc("a<b & \"c\""), "a&lt;b &amp; &quot;c&quot;");
        assert_eq!(oid("10^9/L"), "10_9_L");
    }

    #[test]
    fn shapes() {
        let s = shape_of("X", [ValueRef::Num(1.5), ValueRef::Num(12.25)].into_iter());
        assert_eq!(
            (s.data_type, s.length, s.digits),
            ("float", Some(5), Some(2))
        );
        let s = shape_of(
            "MHSTDTC",
            [ValueRef::Char("2015-06"), ValueRef::Char("2016-01-02")].into_iter(),
        );
        assert_eq!(s.data_type, "partialDate");
        let s = shape_of("AGE", [ValueRef::Int(7), ValueRef::Int(81)].into_iter());
        assert_eq!((s.data_type, s.length), ("integer", Some(2)));
    }
}
