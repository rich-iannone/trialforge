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
