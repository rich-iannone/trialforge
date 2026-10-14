"""Verification level 3: Pointblank's SDTMIG 3.4 conformance engine, on the written package.

Rules that fail because of known Pointblank defects (not the data) are waived, each with the
reason and its entry in Pointblank's BUGS.md, so a real problem is never hidden behind them.
"""

from __future__ import annotations

import dataclasses
import importlib
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from trialforge._study import Finding, Study

# rule id → (values Pointblank wrongly rejects, or None for any dataset-level failure; why)
KNOWN_POINTBLANK_DEFECTS: dict[str, tuple[frozenset[str] | None, str]] = {
    "SDTM-089": (
        frozenset({"DRUG INTERRUPTED", "DRUG WITHDRAWN"}),
        "Pointblank's ACN codelist lacks these CDISC CT (C66767) values (Pointblank BUGS.md #29)",
    ),
    "SDTM-161": (
        None,
        "requires IEINCLCR in TI, not an SDTMIG TI variable (Pointblank BUGS.md #33)",
    ),
    "SDTM-240": (
        None,
        "requires IEINCLCR in TI, not an SDTMIG TI variable (Pointblank BUGS.md #33)",
    ),
    "SDTM-162": (
        None,
        "requires VISITNUM in AE, which is only Permissible (Pointblank BUGS.md #34)",
    ),
    "SDTM-434": (
        frozenset(
            {f"TREATMENT {k}" for k in range(1, 5)}
            | {f"WASHOUT {k}" for k in range(1, 4)}
            | {"CONTINUATION TREATMENT"}
        ),
        "rejects sponsor extensions of the extensible EPOCH codelist that define.xml declares "
        "(def:ExtendedValue) (Pointblank BUGS.md #37), and CONTINUATION TREATMENT (C123452), "
        "missing from its EPOCH list (BUGS.md #40)",
    ),
    "SDTM-066": (
        frozenset({"EVERY WEEK", "EVERY 8 WEEKS", "EVERY 12 WEEKS", "EVERY 16 WEEKS"}),
        "Pointblank's FREQ codelist lacks these CDISC CT (C71113) values (Pointblank BUGS.md #38)",
    ),
    "SDTM-435": (
        frozenset({"EVERY WEEK", "EVERY 8 WEEKS", "EVERY 12 WEEKS", "EVERY 16 WEEKS"}),
        "Pointblank's FREQ codelist lacks these CDISC CT (C71113) values (Pointblank BUGS.md #38)",
    ),
    "SDTM-335": (
        None,
        "requires VISITNUM in every VS record; it is Expected, so null for diary temperatures "
        "(Pointblank BUGS.md #39)",
    ),
    "SDTM-170": (None, "requires LBORNRLO to be numeric; SDTMIG defines it as Char (BUGS.md #35)"),
    "SDTM-171": (None, "requires LBORNRHI to be numeric; SDTMIG defines it as Char (BUGS.md #35)"),
}


def _import_pointblank() -> Any:
    try:
        return importlib.import_module("pointblank")
    except ImportError as e:  # pragma: no cover - exercised only without Pointblank
        raise ImportError(
            "verify(level=3) needs Pointblank's CDISC engine: pip install 'trialforge[verify]'"
        ) from e


def _classify(report: Any) -> tuple[list[Finding], list[Finding]]:
    """Split a Pointblank report's failed rules into (real failures, waived known defects)."""
    from trialforge._study import Finding

    failed: list[Finding] = []
    waived: list[Finding] = []
    for rule in report.rules():
        if rule.status != "fail":
            continue
        allowed, why = KNOWN_POINTBLANK_DEFECTS.get(rule.rule_id, (frozenset(), ""))
        for row in rule.row_findings or [None]:
            value = None if row is None else row.checked_value
            finding = Finding(
                "L3",
                rule.rule_id,
                rule.dataset or "",
                None if row is None else row.checked_column,
                None if row is None else row.row,
                rule.message or rule.description or "",
            )
            if allowed is None or (value is not None and value in allowed):
                waived.append(
                    dataclasses.replace(finding, message=f"{finding.message} [waived: {why}]")
                )
            else:
                failed.append(finding)
    return failed, waived


def pointblank_findings(study: Study) -> tuple[list[Finding], list[Finding]]:
    """Run Pointblank on the study's SDTM package (with define.xml). Returns (failed, waived)."""
    pb = _import_pointblank()
    with tempfile.TemporaryDirectory() as tmp:
        study.write_submission(tmp)
        report = pb.SubmissionPackage.from_folder(str(Path(tmp) / "sdtm")).validate_conformance()
    return _classify(report)


def pointblank_findings_for(datasets: dict[str, Any]) -> list[Finding]:
    """Run Pointblank on in-memory SDTM datasets; returns the failures that are not waived."""
    pb = _import_pointblank()
    return _classify(pb.validate_sdtmig(datasets))[0]
