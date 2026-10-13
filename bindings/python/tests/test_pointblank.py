"""Verification level 3: Pointblank's SDTMIG 3.4 conformance rules (skipped without Pointblank).

Phase 1 exit criterion: `validate_sdtmig` reports zero errors. KNOWN lists rule findings that
are Pointblank defects rather than problems with the data; remove entries as Pointblank fixes
them.
"""

from typing import Any

import pytest

import trialforge as tf

pb = pytest.importorskip("pointblank")

# (rule, column) → (values Pointblank wrongly rejects, why).
KNOWN: dict[tuple[str, str], tuple[set[str], str]] = {
    ("SDTM-089", "AEACN"): (
        {"DRUG INTERRUPTED", "DRUG WITHDRAWN"},
        "Pointblank's bundled ACN codelist lacks these CDISC CT (C66767) submission values",
    ),
}

CASES: list[tuple[str, dict[str, Any]]] = [
    ("parallel_continuous", {}),
    ("parallel_continuous", {"indication": "hypertension", "weeks": 12}),
    ("parallel_continuous", {"indication": "obesity", "dropout_rate": 0.8}),
    ("oncology_solid_tumor", {}),
]


@pytest.mark.parametrize(("archetype", "controls"), CASES)
def test_validate_sdtmig_reports_no_errors(archetype: str, controls: dict[str, Any]) -> None:
    study = tf.study(archetype, n=150, seed=5, **controls)
    findings = pb.validate_sdtmig(dict(study.sdtm), study_id=study.studyid).findings_df()
    unexpected = [
        f
        for f in findings.to_dicts()
        if f["checked_value"] not in KNOWN.get((f["rule_id"], f["checked_column"]), (set(), ""))[0]
    ]
    assert not unexpected, unexpected[:10]
