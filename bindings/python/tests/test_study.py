import datetime
import math
from pathlib import Path
from typing import Any

import pandas as pd
import polars as pl
import pyarrow as pa
import pyreadstat
import pytest

import trialforge as tf


@pytest.fixture(scope="module")
def par() -> tf.Study:
    return tf.study("parallel_continuous", n=80, seed=3)


@pytest.fixture(scope="module")
def onc() -> tf.Study:
    return tf.study("oncology_solid_tumor", n=60, seed=3)


def test_archetypes_and_controls() -> None:
    assert set(tf.archetypes()) == {"parallel_continuous", "oncology_solid_tumor"}
    names = {c.name: c.default for c in tf.controls("parallel_continuous")}
    assert names["indication"] == "type_2_diabetes"
    assert names["weeks"] == 24
    assert "xpt" in tf.formats()


def test_parallel_study_shape(par: tf.Study) -> None:
    assert par.studyid == "SIMPAR01"
    assert par.seed == 3
    assert list(par.sdtm) == [
        "TA", "TE", "TV", "TI", "TS", "DM", "SE", "SV", "EX", "CM", "AE", "DS", "MH", "LB", "VS",
        "SUPPAE",
    ]  # fmt: skip
    assert list(par.adam) == ["ADSL", "ADAE", "ADLB", "ADVS"]
    assert par.sdtm["DM"].height == 80
    assert par.controls["indication"] == "type_2_diabetes"
    assert par.truth["endpoint"]["paramcd"] == "HBA1C"
    assert par.verify().passed


def test_oncology_study_shape(onc: tf.Study) -> None:
    assert onc.studyid == "SIMONC02"
    assert onc.cutoff == datetime.date(2025, 8, 30)  # 30 months of 30.4375 days
    assert "DM" in onc.sdtm and "ADSL" in onc.adam
    assert "STRAT1R" in onc.adam["ADSL"].columns
    assert "overall_survival" in onc.truth


def test_controls_are_applied() -> None:
    s = tf.study("parallel_continuous", n=20, seed=1, indication="hypertension", weeks=12)
    assert s.controls["weeks"] == 12
    assert s.truth["endpoint"]["paramcd"] == "SYSBP"
    assert s.sdtm["TV"]["VISIT"].to_list()[-3] == "WEEK 12"


def test_outputs() -> None:
    def adsl(output: tf.Output) -> Any:
        return tf.study("parallel_continuous", n=20, seed=1, output=output).adam["ADSL"]

    assert isinstance(adsl("polars"), pl.DataFrame)
    assert isinstance(adsl("pandas"), pd.DataFrame)
    table = adsl("pyarrow")
    assert isinstance(table, pa.Table)
    assert table.schema.field("TRTSDT").type == pa.date32()
    assert table.schema.field("USUBJID").metadata == {b"label": b"Unique Subject Identifier"}
    assert table.schema.metadata[b"label"] == b"Subject-Level Analysis Dataset"
    with pytest.raises(ValueError, match="output"):
        adsl("spark")  # pyright: ignore[reportArgumentType]


def test_seed_handling() -> None:
    a = tf.study("parallel_continuous", n=10)
    assert 0 <= a.seed <= tf.MAX_SEED
    b = tf.study("parallel_continuous", n=10, seed=a.seed)
    assert a.sdtm["LB"].equals(b.sdtm["LB"])
    tf.study("parallel_continuous", n=1, seed=tf.MAX_SEED)
    with pytest.raises(ValueError, match="2\\^53"):
        tf.study("parallel_continuous", n=1, seed=tf.MAX_SEED + 1)


def test_bad_input() -> None:
    with pytest.raises(ValueError, match="unknown archetype"):
        tf.study("nope", n=10, seed=1)
    with pytest.raises(ValueError, match="at least 1"):
        tf.study("parallel_continuous", n=0, seed=1)
    with pytest.raises(ValueError, match="unknown control"):
        tf.study("parallel_continuous", n=5, seed=1, wekks=3)
    with pytest.raises(ValueError, match="invalid controls"):
        tf.study("parallel_continuous", n=5, seed=1, weeks="twelve")
    with pytest.raises(ValueError, match="unknown format"):
        tf.study("parallel_continuous", n=5, seed=1).write("unused", formats=["sas7bdat"])


def _expected(values: list[Any], dtype: Any) -> list[Any]:
    """In-memory values as they should read back from XPT."""
    if dtype == pl.String:
        # XPT has no character missing value: missing reads back as "".
        return ["" if v is None else v for v in values]
    return values


@pytest.mark.parametrize("which", ["par", "onc"])
def test_xpt_round_trip(which: str, request: pytest.FixtureRequest, tmp_path: Path) -> None:
    study: tf.Study = request.getfixturevalue(which)
    paths = [p for p in study.write(tmp_path) if p.suffix == ".xpt"]
    frames = {**study.sdtm, **study.adam}
    assert sorted(p.stem.upper() for p in paths) == sorted(frames)
    for path in paths:
        back: Any
        back, meta = pyreadstat.read_xport(str(path))
        name = path.stem.upper()
        mem: pl.DataFrame = frames[name]
        assert meta.table_name == name
        assert list(back.columns) == mem.columns
        assert back.shape == mem.shape
        for col in mem.columns:
            expected = _expected(mem[col].to_list(), mem[col].dtype)
            got = back[col].tolist()
            if mem[col].dtype == pl.Date:
                got = [None if v is None or v != v else v for v in got]
            elif mem[col].dtype != pl.String:
                got = [None if isinstance(v, float) and math.isnan(v) else v for v in got]
            assert got == expected, f"{name}.{col}"


def test_parquet_and_csv_round_trip(par: tf.Study, tmp_path: Path) -> None:
    par.write(tmp_path, ["parquet", "csv"])
    for name, mem in {**par.sdtm, **par.adam}.items():
        assert pl.read_parquet(tmp_path / f"{name.lower()}.parquet").equals(mem), name
        back = pl.read_csv(tmp_path / f"{name.lower()}.csv", schema=mem.schema)
        expected = mem.with_columns(pl.col(pl.String).fill_null(""))
        assert back.with_columns(pl.col(pl.String).fill_null("")).equals(expected), name
    assert (tmp_path / "truth.json").exists()
