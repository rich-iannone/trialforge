"""Verification level 4: greenwood re-analyses ADTTE; skipped without greenwood.

greenwood is validated against R's survival package, so agreement here ties trialforge's truth
manifest (Kaplan-Meier medians, stratified Cox hazard ratios) to an independent implementation.
The coverage of the hazard ratio's confidence interval over replicate studies (the Phase 3 exit
criterion) is a Rust test: `cargo test --release -p trialforge --test coverage -- --ignored`.
"""

import math
from typing import Any

import polars as pl
import pytest

import trialforge as tf

gw = pytest.importorskip("greenwood")

DAYS_PER_MONTH = 30.4375
ENDPOINTS = [("OS", "overall_survival"), ("PFS", "progression_free_survival")]


@pytest.fixture(scope="module")
def trial() -> tf.Study:
    return tf.study("oncology_solid_tumor", n=400, seed=11)


def _surv(t: pl.DataFrame) -> Any:
    return gw.Surv(t["AVAL"].to_numpy(), 1 - t["CNSR"].to_numpy())


@pytest.mark.parametrize(("paramcd", "key"), ENDPOINTS)
def test_cox_hazard_ratio_matches_greenwood(trial: tf.Study, paramcd: str, key: str) -> None:
    t = trial.adam["ADTTE"].filter(pl.col("PARAMCD") == paramcd)
    x = (t["TRT01P"] == "Drug A").cast(pl.Float64).to_numpy().reshape(-1, 1)
    fit = gw.CoxPH(ties="breslow").fit(_surv(t), x, strata=t["STRAT1R"].to_numpy())
    truth = trial.truth[key]["stratified_cox_hr"]
    assert fit.hazard_ratio_[0] == pytest.approx(truth["hr"], abs=1e-4)
    # greenwood reports the interval on the coefficient scale.
    ci = [math.exp(float(c[0])) for c in (fit.conf_low_, fit.conf_high_)]
    assert ci == pytest.approx(truth["ci95"], abs=1e-4)


@pytest.mark.parametrize(("paramcd", "key"), ENDPOINTS)
def test_km_medians_match_greenwood(trial: tf.Study, paramcd: str, key: str) -> None:
    t = trial.adam["ADTTE"].filter(pl.col("PARAMCD") == paramcd)
    for arm, summary in trial.truth[key]["by_arm"].items():
        a = t.filter(pl.col("TRT01P") == arm)
        median = gw.KaplanMeier().fit(_surv(a)).median()["time"][0]
        assert summary["n"] == a.height
        assert summary["events"] == int((a["CNSR"] == 0).sum())
        assert summary["km_median_months"] == pytest.approx(median / DAYS_PER_MONTH, abs=1e-4)


def test_treatment_effect_direction(trial: tf.Study) -> None:
    for _, key in ENDPOINTS:
        assert trial.truth[key]["stratified_cox_hr"]["hr"] < 1
