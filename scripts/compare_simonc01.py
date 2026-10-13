"""Distributional comparison of trialforge's oncology archetype with SIMONC01.

SIMONC01 is the simulated solid-tumour study of the clinical-trials-info book
(`data-raw/sim-onco.R`), designed with the same targets as trialforge's defaults: 300 subjects
1:1, 18 months' enrolment, 30 months' follow-up, CT every 8 weeks, confirmed ORR 35% vs 20%,
median PFS 7 vs 4.5 months and median OS 20 vs 15 months.

The script prints each summary for SIMONC01 next to trialforge's mean and range over replicate
studies of the same size, so differences can be read against sampling variation.

Usage (from the repo root, after `make py-install`):

    .venv/bin/python scripts/compare_simonc01.py [SIMONC01-dir] [--reps N]

SIMONC01-dir defaults to ~/other_projects/clinical-trials-info/data-raw/sim-onco (with sdtm/
and adam/ Parquet files).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl
import trialforge as tf

Frames = dict[str, pl.DataFrame]
DAYS_PER_MONTH = 30.4375


def km_median(time: pl.Series, event: pl.Series) -> float | None:
    """Kaplan-Meier median (the first time survival drops to 0.5 or below)."""
    df = pl.DataFrame({"t": time, "e": event}).group_by("t").agg(pl.col("e").sum(), pl.len()).sort("t")
    at_risk = len(time)
    s = 1.0
    for t, d, n in df.iter_rows():
        if d:
            s *= 1 - d / at_risk
        if s <= 0.5:
            return float(t)
        at_risk -= n
    return None


def tte(adam: Frames, paramcd: str, arm: str) -> pl.DataFrame:
    return adam["ADTTE"].filter(pl.col("PARAMCD") == paramcd, pl.col("TRT01P") == arm)


def metrics(sdtm: Frames, adam: Frames) -> dict[str, float | None]:
    n = adam["ADSL"].height
    tu, tr, rs, adrs, adsl = sdtm["TU"], sdtm["TR"], sdtm["RS"], adam["ADRS"], adam["ADSL"]
    targets = tu.filter(pl.col("TUSTRESC") == "TARGET")
    ovr = rs.filter(pl.col("RSTESTCD") == "OVRLRESP")
    target_tr = tr.filter(pl.col("TRGRPID") == "TARGET", pl.col("TRTESTCD") != "SUMDIAM")
    baseline_sum = tr.filter(pl.col("TRTESTCD") == "SUMDIAM", pl.col("VISIT") == "SCREENING")
    out: dict[str, float | None] = {
        "target lesions per subject": targets.height / n,
        "non-target lesions per subject": tu.filter(pl.col("TUSTRESC") == "NON-TARGET").height / n,
        "% target lesions nodal": 100
        * targets.filter(pl.col("TULOC") == "LYMPH NODE").height
        / targets.height,
        "baseline sum of diameters, median (mm)": baseline_sum["TRSTRESN"].median(),  # type: ignore[arg-type]
        "post-baseline assessments per subject": ovr.height / n,
        "% target measurements not done": 100
        * target_tr.filter(pl.col("TRSTAT") == "NOT DONE").height
        / target_tr.height,
        "% subjects with a new lesion": 100
        * tu.filter(pl.col("TUSTRESC") == "NEW")["USUBJID"].n_unique()
        / n,
    }
    for r in ["CR", "PR", "SD", "PD", "NE"]:
        out[f"% overall responses {r}"] = (
            100 * ovr.filter(pl.col("RSSTRESC") == r).height / max(ovr.height, 1)
        )
    for arm in ["Drug A", "Control"]:
        cbor = adrs.filter(pl.col("PARAMCD") == "CBOR", pl.col("TRT01P") == arm)
        bor = adrs.filter(pl.col("PARAMCD") == "BOR", pl.col("TRT01P") == arm)
        out[f"confirmed ORR % ({arm})"] = 100 * cbor["AVALC"].is_in(["CR", "PR"]).mean()  # type: ignore[operator]
        out[f"unconfirmed ORR % ({arm})"] = 100 * bor["AVALC"].is_in(["CR", "PR"]).mean()  # type: ignore[operator]
        for p in ["OS", "PFS"]:
            t = tte(adam, p, arm)
            m = km_median(t["AVAL"], 1 - t["CNSR"])
            out[f"{p} KM median, months ({arm})"] = None if m is None else m / DAYS_PER_MONTH
            out[f"{p} events % ({arm})"] = 100 * (t["CNSR"] == 0).mean()  # type: ignore[operator]
    out["% died"] = 100 * (adsl["DTHFL"] == "Y").fill_null(False).mean()  # type: ignore[operator]
    out["% withdrew (other than death)"] = (
        100 * adsl.filter(pl.col("EOSSTT") == "DISCONTINUED", pl.col("DCSREAS") != "DEATH").height / n
    )
    return out


def read(dir: Path) -> tuple[Frames, Frames]:
    def load(kind: str) -> Frames:
        return {p.stem.upper(): pl.read_parquet(p) for p in sorted((dir / kind).glob("*.parquet"))}

    return load("sdtm"), load("adam")


def fmt(x: float | None) -> str:
    return "–" if x is None else f"{x:.1f}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    default = Path.home() / "other_projects/clinical-trials-info/data-raw/sim-onco"
    ap.add_argument("simonc01", nargs="?", type=Path, default=default)
    ap.add_argument("--reps", type=int, default=20)
    args = ap.parse_args()

    reference = metrics(*read(args.simonc01))
    runs: list[dict[str, float | None]] = []
    for seed in range(args.reps):
        s = tf.study("oncology_solid_tumor", n=300, seed=seed + 1)
        runs.append(metrics(s.sdtm, s.adam))

    def agg(key: str) -> tuple[float, float, float] | None:
        xs = [r[key] for r in runs if r[key] is not None]
        if not xs:
            return None
        values = [float(x) for x in xs if x is not None]
        return sum(values) / len(values), min(values), max(values)

    print(f"| Summary | SIMONC01 | trialforge mean [range over {args.reps} studies of 300] |")
    print("|---|---:|---|")
    for key, ref in reference.items():
        a = agg(key)
        tfs = "–" if a is None else f"{a[0]:.1f} [{a[1]:.1f}, {a[2]:.1f}]"
        flag = "" if a is None or ref is None or a[1] <= ref <= a[2] else " ←"
        print(f"| {key} | {fmt(ref)} | {tfs}{flag} |")


if __name__ == "__main__":
    main()
