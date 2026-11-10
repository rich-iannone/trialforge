"""Verification level 3: Pointblank's SDTMIG 3.4 conformance engine, on the written package.

Rules that fail because of known Pointblank defects (not the data) are waived, each with the
reason and its entry in Pointblank's BUGS.md, so a real problem is never hidden behind them.
There are none at present: Pointblank 1.1.1 fixed every defect trialforge had found
(posit-dev/pointblank#455 and #462).
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
KNOWN_POINTBLANK_DEFECTS: dict[str, tuple[frozenset[str] | None, str]] = {}


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
