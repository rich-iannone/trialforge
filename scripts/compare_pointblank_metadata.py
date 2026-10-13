"""Compare trialforge's CDISC metadata bundle with Pointblank's SDTM/ADaM templates.

Run from the repo root with an interpreter that has Pointblank installed:

    python scripts/compare_pointblank_metadata.py

Reports variables whose label or type differs where both define them. Known, intended
differences are listed in KNOWN and do not fail the run.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

from pointblank.metadata._adam_templates import get_adam_dataset, list_adam_datasets
from pointblank.metadata._sdtm_templates import get_sdtm_domain, list_sdtm_domains

ASSETS = Path(__file__).resolve().parents[1] / "crates" / "trialforge-cdisc" / "assets"

# (dataset, variable): why trialforge differs on purpose.
KNOWN: dict[tuple[str, str], str] = {
    ("ADSL", "DCTREAS"): "Pointblank's label is 41 characters (over the XPT limit); we use admiral's wording",
    ("LB", "LBORNRLO"): "SDTMIG 3.4 wording ('... in Orig Unit')",
    ("LB", "LBORNRHI"): "SDTMIG 3.4 wording ('... in Orig Unit')",
    ("VS", "VSTESTCD"): "SDTMIG 3.4 label",
    ("VS", "VSTEST"): "SDTMIG 3.4 label",
    ("VS", "VSPOS"): "SDTMIG 3.4 label",
    ("VS", "VSDTC"): "SDTMIG 3.4 label",
    ("VS", "VSDY"): "SDTMIG 3.4 label",
    ("DS", "DSSTDY"): "SDTMIG 3.4 label",
    ("MH", "MHSTDTC"): "SDTMIG 3.4 label",
    ("MH", "MHENDTC"): "SDTMIG 3.4 label",
    ("ADAE", "ADURN"): "ADaM OCCDS label",
    ("ADAE", "ADURU"): "ADaM OCCDS label",
    ("ADAE", "TRTEMFL"): "ADaMIG 1.3 label",
    **{
        (ds, var): "ADaMIG 1.3 label"
        for ds in ("ADLB", "ADVS")
        for var in ("ANL01FL", "SRCDOM", "SRCVAR", "SRCSEQ")
    },
}


def load(standard: str) -> dict[tuple[str, str], dict[str, str]]:
    with open(ASSETS / standard / "variables.tsv", newline="") as f:
        return {(r["dataset"], r["variable"]): r for r in csv.DictReader(f, delimiter="\t")}


def main() -> int:
    ours = {**load("sdtmig-3.4"), **load("adamig-1.3")}
    theirs: dict[tuple[str, str], tuple[str, str]] = {}
    for d in list_sdtm_domains():
        for v in get_sdtm_domain(d).variables:
            theirs[(d, v.name)] = (v.label, v.dtype)
    bds_targets = ["ADLB", "ADVS"]
    for d in list_adam_datasets():
        for target in bds_targets if d == "BDS" else [d]:
            for v in get_adam_dataset(d).variables:
                theirs[(target, v.name)] = (v.label, getattr(v, "dtype", ""))
    problems = 0
    for key in sorted(ours.keys() & theirs.keys()):
        mine, (label, dtype) = ours[key], theirs[key]
        diffs = []
        if mine["label"] != label:
            diffs.append(f"label {mine['label']!r} vs Pointblank {label!r}")
        if dtype in ("Char", "Num") and mine["type"] != dtype:
            diffs.append(f"type {mine['type']} vs Pointblank {dtype}")
        if diffs:
            note = KNOWN.get(key)
            status = f"known: {note}" if note else "UNEXPECTED"
            problems += note is None
            print(f"{key[0]}.{key[1]}: {'; '.join(diffs)} [{status}]")
    shared = len(ours.keys() & theirs.keys())
    print(f"{shared} shared variables compared; {problems} unexpected difference(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
