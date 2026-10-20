"""Replay the shared conformance cases in spec/conformance/ (also replayed by the R suite)."""

import hashlib
import json
from pathlib import Path
from typing import Any

import pytest

import trialforge

SPEC = Path(__file__).resolve().parents[3] / "spec" / "conformance"


def _cases(name: str) -> list[dict[str, Any]]:
    return json.loads((SPEC / name).read_text())["cases"]


@pytest.mark.parametrize("case", _cases("rng_streams.json"))
def test_rng_streams(case: dict[str, Any]) -> None:
    assert str(trialforge.stream_seed(case["seed"], case["path"])) == case["expected"]


@pytest.mark.parametrize("case", _cases("study_xpt.json"))
def test_study_xpt(case: dict[str, Any], tmp_path: Path) -> None:
    study = trialforge.study(case["archetype"], n=case["n"], seed=case["seed"], **case["controls"])
    rows = {name: frame.height for name, frame in {**study.sdtm, **study.adam}.items()}
    assert rows == case["rows"]
    written = {
        p.name: hashlib.md5(p.read_bytes()).hexdigest()
        for p in study.write(tmp_path)
        if p.suffix == ".xpt"
    }
    assert written == case["xpt_md5"]
    defines = {
        std: hashlib.md5(study.define_xml(std).encode()).hexdigest() for std in ("sdtm", "adam")
    }
    assert defines == case["define_md5"]


def _expected(defect: dict[str, Any]) -> str:
    return "; ".join(f"{e['rule']} {e['dataset']} {e['row']}" for e in defect["expected"])


@pytest.mark.parametrize("case", _cases("defects.json"))
def test_defects(case: dict[str, Any], tmp_path: Path) -> None:
    study = trialforge.study(
        case["archetype"],
        n=case["n"],
        seed=case["seed"],
        defects=case["defects"],
        **case["controls"],
    )
    written = study.write(tmp_path)
    manifest = case["manifest"]
    assert json.loads((tmp_path / "defects.json").read_text()) == manifest
    md5 = {p.name: hashlib.md5(p.read_bytes()).hexdigest() for p in written}
    assert {name: md5[name] for name in case["xpt_md5"]} == case["xpt_md5"]
    # The table's columns are the manifest's fields, `expected` flattened to text.
    table = study.defects
    assert table is not None
    assert table.rows() == [
        (*(d[c] for c in table.columns[:-1]), _expected(d)) for d in manifest["defects"]
    ]
    score = study.score(study.verify())
    assert {k: getattr(score, k) for k in case["score"]} == case["score"]
    assert (score.recall, score.precision) == (1.0, 1.0)
