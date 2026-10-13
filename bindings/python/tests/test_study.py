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
    assert {"TU", "TR", "RS"} <= set(onc.sdtm)
    assert list(onc.adam) == ["ADSL", "ADRS", "ADTTE", "ADTR"]
    assert set(onc.adam["ADTTE"]["PARAMCD"].unique()) >= {"OS", "PFS"}
    assert onc.truth["targets"]["os_hr"] == 0.75
    assert "stratified_cox_hr" in onc.truth["overall_survival"]
    assert "confirmed_orr" in onc.truth


def test_oncology_controls() -> None:
    s = tf.study("oncology_solid_tumor", n=30, seed=2, orr=[0.5, 0.1], os_median_months=20)
    assert s.controls["orr"] == [0.5, 0.1]
    assert s.truth["targets"]["confirmed_orr"] == {"Drug A": 0.5, "Control": 0.1}
    with pytest.raises(ValueError, match="pfs_median_months"):
        tf.study("oncology_solid_tumor", n=5, seed=1, pfs_median_months=30)


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


def test_write_submission_layout(tmp_path: Path) -> None:
    study = tf.study("parallel_continuous", n=30, seed=4)
    written = study.write_submission(tmp_path)
    assert (tmp_path / "sdtm" / "define.xml") in written
    assert (tmp_path / "adam" / "define.xml") in written
    assert sorted(p.stem.upper() for p in (tmp_path / "sdtm").glob("*.xpt")) == sorted(study.sdtm)
    assert sorted(p.stem.upper() for p in (tmp_path / "adam").glob("*.xpt")) == sorted(study.adam)
    assert study.define_xml("adam").startswith('<?xml version="1.0" encoding="UTF-8"?>')
    with pytest.raises(ValueError, match="sdtm"):
        study.define_xml("send")  # pyright: ignore[reportArgumentType]


@pytest.mark.parametrize("archetype", ["parallel_continuous", "oncology_solid_tumor"])
@pytest.mark.parametrize("standard", ["sdtm", "adam"])
def test_define_xml_is_schema_valid(archetype: str, standard: str) -> None:
    xmlschema = pytest.importorskip("xmlschema")
    odmlib = pytest.importorskip("odmlib")
    xsd = Path(odmlib.__file__).parent / "schemas" / "define" / "2.1" / "define2-1-0.xsd"
    schema = xmlschema.XMLSchema(str(xsd))
    xml = tf.study(archetype, n=40, seed=6).define_xml(standard)  # pyright: ignore[reportArgumentType]
    errors = [str(e.reason) for e in schema.iter_errors(xml)]
    assert not errors, errors[:5]


def test_dataset_json_round_trip(par: tf.Study, tmp_path: Path) -> None:
    import json

    par.write(tmp_path, ["dataset-json"])
    for name, mem in {**par.sdtm, **par.adam}.items():
        doc = json.loads((tmp_path / f"{name.lower()}.json").read_text())
        assert doc["datasetJSONVersion"] == "1.1.0"
        assert doc["itemGroupOID"] == f"IG.{name}"
        assert doc["records"] == mem.height == len(doc["rows"])
        assert [c["name"] for c in doc["columns"]] == mem.columns
        for j, col in enumerate(mem.columns):
            got = [row[j] for row in doc["rows"]]
            expected = mem[col].to_list()
            if mem[col].dtype == pl.Date:
                expected = [None if v is None else v.isoformat() for v in expected]
            assert got == expected, f"{name}.{col}"


def test_study_from_spec_matches_study(tmp_path: Path) -> None:
    yaml = (Path(__file__).resolve().parents[3] / "examples" / "hypertension.yaml").read_text()
    from_file = tf.study_from_spec(
        Path(__file__).resolve().parents[3] / "examples" / "hypertension.yaml"
    )
    from_text = tf.study_from_spec(yaml)
    direct = tf.study(
        "parallel_continuous", n=240, seed=20261013, studyid="HTN001",
        indication="hypertension", weeks=12, change=[-4, -11], dropout_rate=0.12,
    )  # fmt: skip
    assert from_file.sdtm["LB"].equals(direct.sdtm["LB"])
    assert from_text.adam["ADSL"].equals(direct.adam["ADSL"])
    as_dict = tf.study_from_spec(
        {"study": {"archetype": "parallel_continuous", "seed": 1}, "population": {"n": 10}}
    )
    assert as_dict.sdtm["DM"].height == 10
    with pytest.raises(ValueError, match="not supported yet"):
        tf.study_from_spec(
            {"study": {"archetype": "parallel_continuous", "seed": 1}, "dropout": {}}
        )
