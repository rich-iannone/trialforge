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
def onc() -> tf.Study:
    return tf.study("oncology_solid_tumor", n=120, seed=3)


def test_archetypes() -> None:
    assert "oncology_solid_tumor" in tf.archetypes()


def test_study_shape(onc: tf.Study) -> None:
    assert onc.studyid == "SIMONC02"
    assert onc.archetype == "oncology_solid_tumor"
    assert onc.seed == 3
    assert onc.cutoff == datetime.date(2025, 8, 29)
    assert list(onc.sdtm) == ["DM", "EX", "DS"]
    assert list(onc.adam) == ["ADSL"]
    assert onc.sdtm["DM"].height == 120
    assert "SIMONC02" in repr(onc)


def test_outputs() -> None:
    def adsl(output: tf.Output) -> Any:
        return tf.study("oncology_solid_tumor", n=20, seed=1, output=output).adam["ADSL"]

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
    a = tf.study("oncology_solid_tumor", n=10)
    assert 0 <= a.seed <= tf.MAX_SEED
    b = tf.study("oncology_solid_tumor", n=10, seed=a.seed)
    assert a.sdtm["DM"].equals(b.sdtm["DM"])
    tf.study("oncology_solid_tumor", n=1, seed=tf.MAX_SEED)
    with pytest.raises(ValueError, match="2\\^53"):
        tf.study("oncology_solid_tumor", n=1, seed=tf.MAX_SEED + 1)


def test_bad_input() -> None:
    with pytest.raises(ValueError, match="unknown archetype"):
        tf.study("nope", n=10, seed=1)
    with pytest.raises(ValueError, match="at least 1"):
        tf.study("oncology_solid_tumor", n=0, seed=1)
    with pytest.raises(ValueError, match="unknown format"):
        tf.study("oncology_solid_tumor", n=5, seed=1).write("unused", formats=["sas7bdat"])


def test_xpt_round_trip(onc: tf.Study, tmp_path: Path) -> None:
    paths = onc.write(tmp_path)
    assert sorted(p.name for p in paths) == ["adsl.xpt", "dm.xpt", "ds.xpt", "ex.xpt"]
    frames = {**onc.sdtm, **onc.adam}
    for path in paths:
        back: Any
        back, meta = pyreadstat.read_xport(str(path))
        name = path.stem.upper()
        mem: pl.DataFrame = frames[name]
        assert meta.table_name == name
        assert list(back.columns) == mem.columns
        assert back.shape == mem.shape
        for col in mem.columns:
            expected = mem[col].to_list()
            got = back[col].tolist()
            if mem[col].dtype == pl.String:
                # XPT has no character missing value: missing reads back as "".
                expected = ["" if v is None else v for v in expected]
            elif mem[col].dtype == pl.Date:
                got = [None if v is None or v != v else v for v in got]
            else:
                got = [None if isinstance(v, float) and math.isnan(v) else v for v in got]
            assert got == expected, f"{name}.{col}"
