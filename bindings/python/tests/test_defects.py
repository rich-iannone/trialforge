"""Defect injection: `study(defects=...)`, `Study.defects` and `Study.score()`."""

import json
from pathlib import Path

import polars as pl
import pytest

import trialforge as tf

SPEC = {"ae_end_before_start": 3, "duplicate_seq": 3, "out_of_ct": 5}


def test_defects_are_caught_exactly() -> None:
    s = tf.study("parallel_continuous", n=60, seed=4, defects=SPEC)
    table = s.defects
    assert isinstance(table, pl.DataFrame)
    assert table.height == 11
    assert table.columns == [
        "id", "dataset", "row", "paired_row", "usubjid", "key", "variable", "original",
        "injected", "expected",
    ]  # fmt: skip
    assert table["id"].value_counts().sort("id").rows() == [
        ("ae_end_before_start", 3), ("duplicate_seq", 3), ("out_of_ct", 5),
    ]  # fmt: skip
    report = s.verify()
    assert not report.passed
    score = s.score(report)
    assert isinstance(score, tf.Score)
    assert (score.injected, score.detected, score.recall) == (11, 11, 1.0)
    assert score.precision == 1.0 and score.unmatched == [] and score.missed == []
    assert set(score.by_type) == set(SPEC)
    ct = score.by_type["out_of_ct"]
    assert isinstance(ct, tf.DefectScore)
    assert (ct.injected, ct.detected, ct.false_positives) == (5, 5, 0)


def test_injected_values_are_in_the_data() -> None:
    s = tf.study("parallel_continuous", n=60, seed=4, defects={"out_of_ct": 4})
    for d in s.defects.iter_rows(named=True):
        assert s.sdtm[d["dataset"]][d["variable"]][d["row"]] == d["injected"]
        assert d["original"] != d["injected"]
        assert d["expected"] == f"L0.CT {d['dataset']} {d['row']}"


def test_clean_study_is_unchanged() -> None:
    clean = tf.study("parallel_continuous", n=40, seed=2)
    assert clean.defects is None
    with pytest.raises(ValueError, match="no injected defects"):
        clean.score(clean.verify())
    none = tf.study("parallel_continuous", n=40, seed=2, defects="none")
    assert none.defects is None
    for name, frame in clean.sdtm.items():
        assert frame.equals(none.sdtm[name]), name
    dirty = tf.study("parallel_continuous", n=40, seed=2, defects={"duplicate_seq": 2})
    assert dirty.truth == clean.truth
    assert not dirty.sdtm[dirty.defects["dataset"][0]].equals(
        clean.sdtm[dirty.defects["dataset"][0]]
    )


def test_outputs_and_files(tmp_path: Path) -> None:
    s = tf.study("parallel_continuous", n=40, seed=2, defects={"out_of_ct": 2}, output="pandas")
    assert len(s.defects) == 2
    written = s.write(tmp_path, ["csv"])
    assert tmp_path / "defects.json" in written
    manifest = json.loads((tmp_path / "defects.json").read_text())
    assert manifest["spec"] == {"out_of_ct": 2}
    assert manifest["rows"] == "zero-based"
    assert len(manifest["defects"]) == 2


def test_rates_and_spec_files() -> None:
    s = tf.study("parallel_continuous", n=60, seed=7, defects={"out_of_ct": 0.05})
    assert s.defects is not None and s.defects.height > 0
    from_spec = tf.study_from_spec(
        {
            "study": {"archetype": "parallel_continuous", "seed": 7},
            "population": {"n": 60},
            "data_quality": {"defects": {"out_of_ct": 0.05}},
        }
    )
    assert from_spec.defects is not None and from_spec.defects.equals(s.defects)


@pytest.mark.parametrize("profile", ["light", "realistic", "adversarial"])
def test_profiles_are_caught_exactly(profile: str, tmp_path: Path) -> None:
    s = tf.study("cns_alzheimers", n=80, seed=3, defects=profile)
    assert s.defects is not None and s.defects.height > 0
    score = s.score(s.verify())
    assert (score.recall, score.precision) == (1.0, 1.0)
    s.write(tmp_path, ["csv"])
    manifest = json.loads((tmp_path / "defects.json").read_text())
    assert manifest["profile"] == profile
    # Types the study has no eligible record for are requested but absent.
    assert set(score.by_type) <= set(manifest["spec"])


def test_bad_defects() -> None:
    with pytest.raises(ValueError, match="unknown defect"):
        tf.study("parallel_continuous", n=10, seed=1, defects={"nope": 1})
    with pytest.raises(ValueError, match="light, realistic, adversarial"):
        tf.study("parallel_continuous", n=10, seed=1, defects="heavy")
    with pytest.raises(ValueError, match="rate"):
        tf.study("parallel_continuous", n=10, seed=1, defects={"out_of_ct": 1.5})


def test_level_3_findings_are_scored() -> None:
    pytest.importorskip("pointblank")
    s = tf.study("parallel_continuous", n=60, seed=4, defects=SPEC)
    report = s.verify(level=3)
    score = s.score(report)
    assert score.findings == len(report.findings)
    # L0-L2 alone already detect every defect.
    assert score.recall == 1.0
    # Pointblank on its own: reported, not required (PLAN §8.4).
    l3 = tf.Report(False, [f for f in report.findings if f.level == "L3"])
    alone = s.score(l3)
    assert alone.findings == len(l3.findings)
    assert alone.recall is not None and 0.0 <= alone.recall <= 1.0
