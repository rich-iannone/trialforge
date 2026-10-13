"""Verification level 3 (Pointblank's SDTMIG 3.4 rules); skipped without Pointblank."""

from typing import Any

import pytest

import trialforge as tf
from trialforge._conformance import KNOWN_POINTBLANK_DEFECTS

pb = pytest.importorskip("pointblank")

CASES: list[tuple[str, dict[str, Any]]] = [
    ("parallel_continuous", {}),
    ("parallel_continuous", {"indication": "hypertension", "weeks": 12}),
    ("parallel_continuous", {"indication": "obesity", "dropout_rate": 0.8}),
    ("oncology_solid_tumor", {}),
]


@pytest.mark.parametrize(("archetype", "controls"), CASES)
def test_level_3_passes(archetype: str, controls: dict[str, Any]) -> None:
    report = tf.study(archetype, n=120, seed=5, **controls).verify(level=3)
    assert report.passed, report.findings[:10]
    # Every waived finding is a documented Pointblank defect.
    assert all(f.rule in KNOWN_POINTBLANK_DEFECTS for f in report.waived)
    assert all("[waived:" in f.message for f in report.waived)


def test_level_3_runs_define_rules() -> None:
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmp:
        tf.study("parallel_continuous", n=40, seed=2).write_submission(tmp)
        package = pb.SubmissionPackage.from_folder(str(Path(tmp) / "sdtm"))
        assert package.define is not None
        rules = package.validate_conformance().rules()
    define_rules = [r for r in rules if "Define-XML" in (r.description or "")]
    assert define_rules and all(r.status == "pass" for r in define_rules)


def test_level_3_catches_a_real_problem() -> None:
    """Corrupting the data must surface as a level-3 failure, not be waived."""
    import polars as pl

    from trialforge import _conformance

    study = tf.study("parallel_continuous", n=20, seed=3)
    bad_dm = study.sdtm["DM"].with_columns(pl.lit("FEMALE").alias("SEX"))
    failed = _conformance.pointblank_findings_for(dict(study.sdtm, DM=bad_dm))
    assert any(f.dataset == "DM" and f.variable == "SEX" for f in failed)


def test_level_must_be_2_or_3() -> None:
    with pytest.raises(ValueError, match="level"):
        tf.study("parallel_continuous", n=5, seed=1).verify(level=1)  # pyright: ignore[reportArgumentType]
