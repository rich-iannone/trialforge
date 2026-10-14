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


@pytest.fixture(scope="module")
def pk() -> tf.Study:
    return tf.study("phase1_sad_mad", n=24, seed=3, part="mad", doses=[50, 100, 200], mad_days=4)


def test_archetypes_and_controls() -> None:
    assert set(tf.archetypes()) == {
        "parallel_continuous",
        "oncology_solid_tumor",
        "phase1_sad_mad",
        "crossover_be",
        "dose_ranging_ph2",
        "cns_alzheimers",
        "cv_outcomes",
        "ophthalmology",
        "pediatric_growth",
        "respiratory_exacerbation",
        "vaccine",
        "dose_escalation",
        "heme_lymphoma",
    }
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


def test_phase1_study_shape(pk: tf.Study) -> None:
    assert pk.studyid == "SIMPK01"
    assert {"PC", "PP"} <= set(pk.sdtm)
    assert {"ADPC", "ADPP"} <= set(pk.adam)
    assert pk.adam["ADPC"].schema["ADTM"] == pl.Datetime("ms")
    ex = pk.sdtm["EX"]
    assert ex["EXSTDTC"].str.len_chars().min() == 16  # dosing times
    assert set(pk.sdtm["PP"]["PPTESTCD"]) >= {"CMAX", "AUCTAU", "ARAUC"}
    assert pk.truth["pk"]["model"]["compartments"] == 1
    assert pk.verify().passed


def test_dose_ranging_and_non_inferiority() -> None:
    dr = tf.study("dose_ranging_ph2", n=60, seed=2, doses=[0, 5, 20])
    assert dr.adam["ADSL"]["TRT01P"].n_unique() == 3
    assert dr.truth["dose_response"]["target"]["ed50"] == 7.5
    ni = tf.study("parallel_continuous", n=40, seed=2, comparator="active", ni_margin=0.3)
    assert ni.truth["non_inferiority"]["margin"] == 0.3
    assert "Comparator 10 mg" in set(ni.adam["ADSL"]["TRT01P"])


def test_cns_study_shape() -> None:
    s = tf.study("cns_alzheimers", n=40, seed=3)
    ft = s.sdtm["FT"]
    assert set(ft["FTCAT"]) == {"ADAS-COG"}
    assert ft["FTTESTCD"].n_unique() == 11
    adas = s.adam["ADADAS"]
    assert "ACTOT" in set(adas["PARAMCD"])
    assert "adas_cog" in s.truth
    assert s.verify().passed


def test_cv_outcomes_study_shape() -> None:
    s = tf.study("cv_outcomes", n=500, seed=3, target_events=40)
    assert set(s.sdtm["CE"]["CECAT"]) == {"ADJUDICATED ENDPOINT EVENT"}
    assert set(s.adam["ADTTE"]["PARAMCD"]) == {"MACE", "CVDEATH", "MI", "STROKE", "DEATH"}
    assert s.truth["outcomes"]["realised_first_mace_events"] >= 40
    assert s.verify().passed


def test_ophthalmology_study_shape() -> None:
    s = tf.study("ophthalmology", n=60, seed=3)
    assert {"OE", "SC"} <= set(s.sdtm)
    oe = s.sdtm["OE"]
    assert set(oe["OETESTCD"]) == {"NUMLCOR", "THICK"}
    assert set(oe["OELAT"]) == {"LEFT", "RIGHT"}
    assert set(s.sdtm["EX"]["EXROUTE"]) == {"INTRAVITREAL"}
    adsl = s.adam["ADSL"]
    assert set(adsl["STUDYEYE"]) <= {"LEFT", "RIGHT"}
    bcva = s.adam["ADBCVA"]
    assert set(bcva["PARAMCD"]) == {"SBCVA", "FBCVA", "SBCVALOG", "FBCVALOG"}
    assert s.truth["bcva"]["non_inferiority"]["margin"] == 4.0
    assert s.verify().passed


def test_pediatric_growth_study_shape() -> None:
    s = tf.study("pediatric_growth", n=60, seed=3)
    dm = s.sdtm["DM"]
    assert dm["BRTHDTC"].null_count() == 0
    assert dm["AGE"].min() >= 3
    advs = s.adam["ADVS"]
    assert {"HGTSDS", "HGTPCTL", "WGTASDS", "BMISDS", "HTVEL"} <= set(advs["PARAMCD"])
    assert "AAGECUR" in advs.columns
    assert s.truth["height_velocity"]["non_inferiority"]["margin"] == 2.0
    assert s.verify().passed


def test_respiratory_study_shape() -> None:
    s = tf.study("respiratory_exacerbation", n=80, seed=3)
    assert set(s.sdtm["CE"]["CECAT"]) == {"EXACERBATION"}
    assert set(s.sdtm["CE"]["CESEV"]) <= {"MODERATE", "SEVERE"}
    assert set(s.sdtm["RE"]["RETESTCD"]) == {"FEV1"}
    assert set(s.adam["ADEXAC"]["PARAMCD"]) == {"EXACN", "SEVEXACN", "FUPYR", "EXACRATE"}
    assert set(s.adam["ADTTE"]["PARAMCD"]) == {"TTFEXAC"}
    assert "STRAT1R" in s.adam["ADSL"].columns
    assert "rate_ratio" in s.truth["exacerbations"]
    assert s.verify().passed


def test_vaccine_study_shape() -> None:
    s = tf.study("vaccine", n=60, seed=3)
    assert {"IS", "FACE"} <= set(s.sdtm)
    assert set(s.sdtm["FACE"]["DOMAIN"]) == {"FA"}
    assert set(s.sdtm["CE"]["CECAT"]) == {"REACTOGENICITY"}
    assert {"ADIS", "ADCE", "ADFACE"} <= set(s.adam)
    assert {"VAX01DT", "VAX02DT", "AP02SDT"} <= set(s.adam["ADSL"].columns)
    assert "gmt_ratio" in s.truth["vaccine"]
    assert s.verify().passed


def test_dose_escalation_study_shape() -> None:
    s = tf.study("dose_escalation", n=30, seed=3, method="3+3")
    assert s.sdtm["DM"].height <= 30
    assert "AEDLT" in set(s.sdtm["SUPPAE"]["QNAM"])
    assert "DLTEVLFL" in s.adam["ADSL"].columns
    assert s.truth["escalation"]["method"] == "3+3"
    assert s.verify().passed


def test_lymphoma_study_shape() -> None:
    s = tf.study("heme_lymphoma", n=60, seed=3)
    rs = s.sdtm["RS"]
    assert set(rs["RSCAT"]) == {"LUGANO CLASSIFICATION"}
    assert set(rs["RSSCAT"]) == {"CT", "PET-CT"}
    assert set(s.adam) == {"ADSL", "ADRS", "ADTTE"}
    assert "orr" in s.truth
    assert s.verify().passed


def test_crossover_on_progression() -> None:
    s = tf.study("oncology_solid_tumor", n=200, seed=4, crossover_prob=0.8)
    assert "CONTINUATION TREATMENT" in set(s.sdtm["SE"]["EPOCH"])
    assert {"CROSSFL", "CROSSDT"} <= set(s.adam["ADSL"].columns)
    assert s.truth["crossover"]["crossed_over"] > 0
    assert s.verify().passed


def test_multiple_cuts() -> None:
    final = tf.study("oncology_solid_tumor", n=200, seed=4)
    interim = tf.study("oncology_solid_tumor", n=200, seed=4, interim_months=10)
    assert interim.cutoff < final.cutoff
    ts = interim.sdtm["TS"].filter(pl.col("TSPARMCD") == "DCUTDESC")
    assert ts["TSVAL"].to_list() == ["INTERIM ANALYSIS"]
    ids = set(interim.adam["ADSL"]["USUBJID"])
    assert ids < set(final.adam["ADSL"]["USUBJID"])
    assert interim.truth["interim"]["n_randomised"] == len(ids)
    cv = tf.study("cv_outcomes", n=600, seed=2, target_events=60, interim_events=30)
    assert cv.truth["outcomes"]["realised_first_mace_events"] == 30
    assert cv.verify().passed


def test_crossover_study_shape() -> None:
    s = tf.study("crossover_be", n=24, seed=4, design="2x2x4")
    ta = s.sdtm["TA"]
    assert ta.filter(pl.col("ARMCD") == "TRTR")["EPOCH"].to_list() == [
        "SCREENING",
        "TREATMENT 1",
        "WASHOUT 1",
        "TREATMENT 2",
        "WASHOUT 2",
        "TREATMENT 3",
        "WASHOUT 3",
        "TREATMENT 4",
        "FOLLOW-UP",
    ]
    adsl = s.adam["ADSL"]
    assert {"TRT04P", "AP04SDT", "TRTSEQP"} <= set(adsl.columns)
    assert "APERIOD" in s.adam["ADPP"].columns
    be = s.truth["bioequivalence"]
    assert be["target"]["gmr_auc"] == 0.95 and "AUCIFO" in be
    assert s.verify().passed


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


@pytest.mark.parametrize("which", ["par", "onc", "pk"])
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
            if mem[col].dtype in (pl.Date, pl.Datetime("ms")):
                got = [None if v is None or v != v else v for v in got]
            elif mem[col].dtype != pl.String:
                got = [None if isinstance(v, float) and math.isnan(v) else v for v in got]
            assert got == expected, f"{name}.{col}"


@pytest.mark.parametrize("which", ["par", "pk"])
def test_parquet_and_csv_round_trip(
    which: str, request: pytest.FixtureRequest, tmp_path: Path
) -> None:
    study: tf.Study = request.getfixturevalue(which)
    study.write(tmp_path, ["parquet", "csv"])
    for name, mem in {**study.sdtm, **study.adam}.items():
        parquet = pl.read_parquet(tmp_path / f"{name.lower()}.parquet")
        assert parquet.with_columns(pl.col(pl.Datetime).cast(pl.Datetime("ms"))).equals(mem), name
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


@pytest.mark.parametrize(
    "archetype", ["parallel_continuous", "oncology_solid_tumor", "phase1_sad_mad", "crossover_be"]
)
@pytest.mark.parametrize("standard", ["sdtm", "adam"])
def test_define_xml_is_schema_valid(archetype: str, standard: str) -> None:
    xmlschema = pytest.importorskip("xmlschema")
    odmlib = pytest.importorskip("odmlib")
    xsd = Path(odmlib.__file__).parent / "schemas" / "define" / "2.1" / "define2-1-0.xsd"
    schema = xmlschema.XMLSchema(str(xsd))
    xml = tf.study(archetype, n=40, seed=6).define_xml(standard)  # pyright: ignore[reportArgumentType]
    errors = [str(e.reason) for e in schema.iter_errors(xml)]
    assert not errors, errors[:5]


@pytest.mark.parametrize("which", ["par", "pk"])
def test_dataset_json_round_trip(
    which: str, request: pytest.FixtureRequest, tmp_path: Path
) -> None:
    import json

    study: tf.Study = request.getfixturevalue(which)
    study.write(tmp_path, ["dataset-json"])
    for name, mem in {**study.sdtm, **study.adam}.items():
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
            elif mem[col].dtype == pl.Datetime("ms"):
                expected = [
                    None if v is None else v.isoformat(timespec="seconds") for v in expected
                ]
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
