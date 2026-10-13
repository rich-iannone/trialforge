"""Replay the shared conformance cases in spec/conformance/ (also replayed by the R suite)."""

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
