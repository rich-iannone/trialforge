"""The `Study` object and `study()` constructor."""

from __future__ import annotations

import datetime
import importlib
import json
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from os import PathLike
from pathlib import Path
from typing import Any, Literal

from trialforge import _core

Output = Literal["polars", "pandas", "pyarrow"]
_OUTPUTS: tuple[str, ...] = ("polars", "pandas", "pyarrow")
_EPOCH = datetime.date(1970, 1, 1)


def study(
    archetype: str,
    *,
    n: int = 300,
    seed: int | None = None,
    studyid: str | None = None,
    output: Output = "polars",
    defects: str | Mapping[str, float] | None = None,
    **controls: Any,
) -> Study:
    """Generate a mock study.

    Builds a complete study of one archetype: the trial design datasets, SDTM, ADaM and the
    truth manifest. The archetype fixes the kind of study (its arms, visits, domains and
    derivations), and its controls set what varies from one study of that kind to the next,
    such as the treatment duration, the effect size or the dropout rate. Every study passes
    trialforge's structural, referential and derivational checks before it is returned.

    Parameters
    ----------
    archetype
        The study type, e.g. `"parallel_continuous"`. `archetypes()` lists them all.
    n
        Number of randomized subjects. In `phase1_sad_mad` it is split evenly across the dose
        cohorts, and in `dose_escalation` it is the maximum sample size.
    seed
        Root seed, a whole number in `[0, 2**53]` (see `MAX_SEED`). The same archetype, `n`,
        seed and controls always give the same study, in Python and R. With `None`, a fresh
        seed is drawn and recorded on the result as `Study.seed`, so the study can always be
        reproduced.
    studyid
        `STUDYID`. Each archetype has a default, such as `"SIMPAR01"`.
    output
        DataFrame type for `Study.sdtm`, `Study.adam` and `Study.defects`: `"polars"`
        (`polars.DataFrame`), `"pandas"` (`pandas.DataFrame`) or `"pyarrow"`
        (`pyarrow.Table`, whose fields carry the variable labels as `label` metadata). Only the
        chosen backend needs to be installed.
    defects
        Data-quality defects to inject into copies of the datasets after the clean study is
        generated and verified:

        - `None` or `"none"`: no defects.
        - A profile: `"light"` (about five defects of the common kinds), `"realistic"` (rates
          per eligible record) or `"adversarial"` (three of every type, each in its subtlest
          form).
        - A mapping of defect id to a count (an integer) or to a rate per eligible record (a
          float in `[0, 1]`), e.g. `{"ae_end_before_start": 3, "out_of_ct": 0.02}`.

        See `Study.defects` and `Study.score()`.
    **controls
        The archetype's controls, as keyword arguments, e.g.
        `indication="hypertension", weeks=12`. `controls()` lists each archetype's controls
        with their descriptions and defaults; any control you leave out keeps its default.
        Values are JSON-compatible: numbers, strings, booleans, `None`, and lists for per-arm
        values (in arm order, e.g. `ratio=[1, 2]`).

    Returns
    -------
    Study
        The generated study.

    Raises
    ------
    ValueError
        If the archetype is unknown, a control is unknown, of the wrong type, out of range or
        inconsistent with another control, `n` is less than 1, the seed is above `MAX_SEED`,
        or the output type or the defects requested aren't among the accepted values.

    See Also
    --------
    study_from_spec : Generate a study from a JSON or YAML study spec.
    controls : The controls an archetype accepts.

    Examples
    --------
    A 12-week hypertension study with 200 subjects:

    >>> import trialforge as tf
    >>> trial = tf.study("parallel_continuous", n=200, seed=1, indication="hypertension", weeks=12)
    >>> trial.truth["arms"]
    ['Placebo', 'Drug A 10 mg']
    >>> list(trial.adam)
    ['ADSL', 'ADAE', 'ADLB', 'ADVS']
    >>> len(trial.adam["ADSL"])
    200

    `Study.controls` records every control the study used, defaults included:

    >>> trial.controls["weeks"], trial.controls["dropout_rate"]
    (12, 0.15)

    Controls are checked before anything is generated, so a misspelled one is an error that
    lists the valid names:

    >>> tf.study("parallel_continuous", wekks=12)
    Traceback (most recent call last):
        ...
    ValueError: unknown control "wekks"; available: ae_rate_multiplier, change, ...

    Ask for seeded defects to test a validator:

    >>> dirty = tf.study("parallel_continuous", n=200, seed=1, defects={"out_of_ct": 3})
    >>> len(dirty.defects)
    3
    """
    if output not in _OUTPUTS:
        raise ValueError(f"`output` must be one of {', '.join(map(repr, _OUTPUTS))}")
    payload = json.dumps(controls) if controls else None
    requested = None if defects is None else json.dumps(defects)
    return Study(_core.study(archetype, n, seed, studyid, payload, requested), output)


def study_from_spec(
    spec: str | PathLike[str] | dict[str, Any], *, output: Output = "polars"
) -> Study:
    """Generate a study from a study spec.

    A study spec is a JSON or YAML document that describes a study: its archetype, seed and
    `STUDYID` (in `study`), the number of subjects (`population.n`), the archetype's
    controls (`controls`) and, optionally, defects (`defects`). It is the same request as
    `study()`'s arguments, in a file that can be versioned and shared, and gives the same
    study. The schema is `spec/study-spec.schema.json` in the repository.

    Parameters
    ----------
    spec
        A path to a JSON or YAML file, the JSON or YAML text itself, or a dict with the
        spec's structure. The spec's `outputs` section is for the command-line tool and is
        ignored here; write files with `Study.write()` or `Study.write_submission()`.
    output
        DataFrame type for `Study.sdtm`, `Study.adam` and `Study.defects`: `"polars"`,
        `"pandas"` or `"pyarrow"` (see `study()`).

    Returns
    -------
    Study
        The generated study.

    Raises
    ------
    ValueError
        If the spec doesn't parse, doesn't match the schema, or names an unknown archetype
        or an invalid control.

    Examples
    --------
    >>> import trialforge as tf
    >>> spec = {
    ...     "study": {"archetype": "parallel_continuous", "seed": 1},
    ...     "population": {"n": 30},
    ...     "controls": {"indication": "obesity"},
    ... }
    >>> trial = tf.study_from_spec(spec)
    >>> trial.controls["indication"], len(trial.adam["ADSL"])
    ('obesity', 30)

    The same spec in YAML:

    >>> text = '''
    ... study: {archetype: parallel_continuous, seed: 1}
    ... population: {n: 30}
    ... controls: {indication: obesity}
    ... '''
    >>> tf.study_from_spec(text).truth == trial.truth
    True
    """
    if output not in _OUTPUTS:
        raise ValueError(f"`output` must be one of {', '.join(map(repr, _OUTPUTS))}")
    if isinstance(spec, dict):
        text = json.dumps(spec)
    elif isinstance(spec, PathLike) or (
        "\n" not in spec and not spec.lstrip().startswith("{") and Path(spec).exists()
    ):
        text = Path(spec).read_text(encoding="utf-8")
    else:
        text = spec
    return Study(_core.study_from_spec(text), output)


_EXTRAS = {"pl": "polars", "pd": "pandas", "pa": "pyarrow"}


def _backend(module: str, extra: str) -> Any:
    """Import a DataFrame backend, or explain which extra installs it."""
    try:
        return importlib.import_module(module)
    except ImportError as e:
        others = " or ".join(f"output={o!r}" for o in _OUTPUTS if o != _EXTRAS[extra])
        raise ImportError(
            f"output={_EXTRAS[extra]!r} needs {module}, which is not installed. Install it with "
            f"pip install 'trialforge[{extra}]', or pass {others} to study() if you have that "
            "backend instead. Writing files (Study.write(), Study.write_submission()) needs "
            "no backend."
        ) from e


def archetypes() -> dict[str, str]:
    """The available study archetypes.

    An archetype is a kind of study: a design, a model of how subjects' courses unfold, and
    the SDTM and ADaM datasets a study of that kind has. Pass its name to `study()`, and see
    `controls()` for what each one lets you change.

    Returns
    -------
    dict[str, str]
        Archetype names mapped to one-line descriptions.

    Examples
    --------
    >>> import trialforge as tf
    >>> len(tf.archetypes())
    13
    >>> tf.archetypes()["vaccine"]
    'Placebo-controlled vaccine study: neutralizing antibody titers (IS), 7-day reactogenicity diaries (CE, FACE), ADIS/ADCE/ADFACE'
    """  # noqa: E501
    return dict(_core.archetypes())


@dataclass(frozen=True)
class Control:
    """One archetype control, as listed by `controls()`.

    Attributes
    ----------
    name
        The control's name: the keyword argument of `study()` and the key in a study spec's
        `controls` section.
    description
        What the control sets, its unit, and the values it accepts where they are a fixed
        set.
    default
        The value the archetype uses when the control is left out, as a JSON-compatible
        value: a number, string, boolean, list (one entry per arm, in arm order, for per-arm
        values) or `None` (for an optional setting that is off, or that the archetype derives
        from another control, as the description says).
    """

    name: str
    description: str
    default: Any


def controls(archetype: str) -> list[Control]:
    """The controls an archetype accepts, with descriptions and defaults.

    Controls are the settings that vary from one study of a kind to the next, in the terms a
    protocol or statistical analysis plan would use: the treatment duration, the expected
    effect, the dropout rate, the allocation ratio. Some are clinical targets (a median
    survival, a hazard ratio, a response rate) that trialforge calibrates its model to, and
    the truth manifest reports them next to what the study realized.

    Pass controls to `study()` as keyword arguments, or in a study spec's `controls`
    section. Every control has a default, so `study(archetype)` alone gives a complete
    study, and `Study.controls` records the full set a study used. The user guide lists
    [every archetype's controls](../user-guide/archetypes.html#all-controls-by-archetype).

    Parameters
    ----------
    archetype
        The archetype's name, as listed by `archetypes()`.

    Returns
    -------
    list[Control]
        One `Control` per control, with its name, description and default.

    Raises
    ------
    ValueError
        If the archetype is unknown.

    Examples
    --------
    >>> import trialforge as tf
    >>> [c.name for c in tf.controls("parallel_continuous")]
    ['indication', 'weeks', 'change', 'ratio', 'dose', 'dropout_rate', 'ae_rate_multiplier', 'missed_visit', 'comparator', 'comparator_dose', 'ni_margin']
    >>> tf.controls("parallel_continuous")[1]
    Control(name='weeks', description='Treatment duration in weeks; visits at Weeks 2, 4 and every 4 weeks after', default=24)

    Set a control by passing it to `study()`:

    >>> trial = tf.study("parallel_continuous", n=50, seed=1, weeks=12)
    >>> trial.controls["weeks"]
    12
    """  # noqa: E501
    return [
        Control(c["name"], c["description"], c["default"])
        for c in json.loads(_core.controls(archetype))
    ]


def formats() -> list[str]:
    """The file formats `Study.write()` supports.

    - `"xpt"`: SAS transport (XPORT) version 5, as FDA submissions expect, with variable
      labels and `DATE9.` dates.
    - `"csv"`: RFC 4180 CSV with a header of variable names; missing values are empty and
      dates are ISO 8601.
    - `"dataset-json"`: CDISC Dataset-JSON 1.1 (`.json` files), with variable metadata
      whose OIDs match the study's define.xml.
    - `"parquet"`: Parquet, with the variable labels as Arrow field metadata.

    Each format writes one file per dataset. Define-XML 2.1 is written by
    `Study.write_submission()`, or returned as text by `Study.define_xml()`.

    Returns
    -------
    list[str]
        The format names, to pass to `Study.write()`.

    Examples
    --------
    >>> import trialforge as tf
    >>> tf.formats()
    ['xpt', 'csv', 'dataset-json', 'parquet']
    """
    return _core.formats()


@dataclass(frozen=True)
class Finding:
    """One failed verification check, as listed in a `Report`.

    Attributes
    ----------
    level
        The verification level: `"L0"` (each dataset against the SDTMIG or ADaMIG metadata
        and controlled terminology), `"L1"` (references and consistency between datasets),
        `"L2"` (ADaM derivations, recomputed from SDTM) or `"L3"` (Pointblank's conformance
        rules, from `Study.verify(level=3)`).
    rule
        The rule identifier: for levels 0 to 2, the level and a short name, e.g. `"L0.CT"`
        or `"L1.ORDER"`; for level 3, Pointblank's own rule identifier.
    dataset
        The dataset the finding is in, e.g. `"AE"`.
    variable
        The variable at fault, or `None` for a finding about a whole record or dataset.
    row
        The zero-based row of the record in `Study.sdtm` or `Study.adam`, or `None` for a
        finding about a whole dataset.
    message
        A description of what failed.
    """

    level: str
    rule: str
    dataset: str
    variable: str | None
    row: int | None
    message: str


@dataclass(frozen=True)
class Report:
    """The result of `Study.verify()`.

    Attributes
    ----------
    passed
        `True` when there are no findings. A study without defects always passes levels 0 to
        2.
    findings
        The failed checks, one `Finding` each.
    waived
        Level 3 findings caused by known defects in Pointblank's rule catalog (each message
        names the defect). They don't affect `passed` and aren't scored.

    Examples
    --------
    >>> import trialforge as tf
    >>> tf.study("parallel_continuous", n=50, seed=1).verify()
    Report(passed=True, findings=[], waived=[])
    """

    passed: bool
    findings: list[Finding]
    waived: list[Finding] = field(default_factory=list[Finding])


@dataclass(frozen=True)
class DefectScore:
    """Recall and precision for one defect type, in `Score.by_type`.

    Attributes
    ----------
    injected
        Defects of this type in the study.
    detected
        How many of them at least one finding matched.
    recall
        `detected / injected`, or `None` when none were injected.
    findings
        Findings that matched a defect of this type.
    false_positives
        Findings no defect explains, from a level 0 to 2 rule that catches this type.
    precision
        `findings / (findings + false_positives)`, or `None` when both are 0.
    """

    injected: int
    detected: int
    recall: float | None
    findings: int
    false_positives: int
    precision: float | None


@dataclass(frozen=True)
class Score:
    """How well a report's findings recover a study's injected defects, from `Study.score()`.

    A finding matches a defect when it names the dataset and row the defect expects a finding
    at and, for trialforge's own levels (L0–L2), the expected rule. Findings at level 3 (or
    from another validator) match on dataset and row only, or on the defect's paired row.

    Attributes
    ----------
    injected
        Defects in the study.
    detected
        Defects at least one finding matched.
    recall
        `detected / injected`, or `None` when no defects were injected.
    findings
        Findings scored (waived findings are not).
    matched
        Findings that matched at least one defect.
    precision
        `matched / findings`, or `None` when there are no findings.
    by_type
        A `DefectScore` per defect type.
    missed
        Row indices in `Study.defects` of the defects no finding matched.
    unmatched
        The findings no defect explains.

    Examples
    --------
    >>> import trialforge as tf
    >>> dirty = tf.study("parallel_continuous", n=200, seed=1, defects={"out_of_ct": 3})
    >>> score = dirty.score(dirty.verify())
    >>> score.recall, score.precision
    (1.0, 1.0)
    >>> score.by_type["out_of_ct"]
    DefectScore(injected=3, detected=3, recall=1.0, findings=3, false_positives=0, precision=1.0)
    """

    injected: int
    detected: int
    recall: float | None
    findings: int
    matched: int
    precision: float | None
    by_type: dict[str, DefectScore]
    missed: list[int]
    unmatched: list[Finding]


class Study:
    """A generated mock study: its datasets, the request that made it, and its truth manifest.

    Create one with `study()` or `study_from_spec()`, not with the constructor. A study holds
    the trial design and SDTM datasets (`sdtm`), the ADaM datasets (`adam`), the controls it
    was generated with (`controls`), the truth manifest (`truth`) and, when defects were
    requested, the injected defects (`defects`). Its methods write the datasets to files,
    return define.xml, re-run the checks and score a validator's findings.

    Datasets are converted to the `output` type chosen in `study()` the first time they are
    accessed, then kept.

    Examples
    --------
    >>> import trialforge as tf
    >>> trial = tf.study("parallel_continuous", n=200, seed=1, indication="hypertension", weeks=12)
    >>> trial
    <trialforge.Study SIMPAR01 (parallel_continuous), seed=1, cutoff=2025-05-05: TA, TE, TV, TI, TS, DM, SE, SV, EX, CM, AE, DS, MH, LB, VS, SUPPAE, ADSL, ADAE, ADLB, ADVS>
    >>> trial.studyid, trial.archetype, trial.seed
    ('SIMPAR01', 'parallel_continuous', 1)
    """  # noqa: E501

    def __init__(self, core: _core.Study, output: Output) -> None:
        self._core = core
        self._output: Output = output
        self._cache: dict[str, Any] = {}

    def _cached(self, key: str, make: Callable[[], Any]) -> Any:
        if key not in self._cache:
            self._cache[key] = make()
        return self._cache[key]

    @property
    def studyid(self) -> str:
        """`STUDYID`, the study identifier in every dataset.

        The `studyid` passed to `study()`, or the archetype's default.
        """
        return self._core.studyid

    @property
    def archetype(self) -> str:
        """The archetype the study was generated from, e.g. `"parallel_continuous"`."""
        return self._core.archetype

    @property
    def seed(self) -> int:
        """The root seed.

        The seed passed to `study()`, or the one drawn for it when none was. The same
        archetype, `n`, seed and controls give the same study again:

        >>> import trialforge as tf
        >>> trial = tf.study("vaccine", n=50)
        >>> again = tf.study("vaccine", n=50, seed=trial.seed)
        >>> again.truth == trial.truth
        True
        """
        return self._core.seed

    @property
    def cutoff(self) -> datetime.date:
        """The data cutoff date.

        For an event-driven study (`oncology_solid_tumor`, `heme_lymphoma`, `cv_outcomes`) it
        is the data cutoff of the analysis, or of the interim cut when one was requested; for
        a fixed-duration study it is the last subject's last contact. It is also the creation
        date recorded in the XPT headers and define.xml, so files don't depend on the date
        they were written.

        >>> import trialforge as tf
        >>> tf.study("parallel_continuous", n=200, seed=1, weeks=12).cutoff
        datetime.date(2025, 5, 5)
        """
        return _EPOCH + datetime.timedelta(days=self._core.cutoff_days)

    @property
    def controls(self) -> dict[str, Any]:
        """The full set of controls the study used: the ones you passed, and the defaults of the rest.

        Together with `archetype`, `n`, `seed` and `studyid`, they reproduce the study:
        `study(trial.archetype, n=..., seed=trial.seed, **trial.controls)`.

        >>> import trialforge as tf
        >>> trial = tf.study("parallel_continuous", n=200, seed=1, indication="hypertension", weeks=12)
        >>> trial.controls
        {'ae_rate_multiplier': 1.5, 'change': None, 'comparator': 'placebo', 'comparator_dose': 10.0, 'dose': 10.0, 'dropout_rate': 0.15, 'indication': 'hypertension', 'missed_visit': 0.04, 'ni_margin': None, 'ratio': [1, 1], 'weeks': 12}
        """  # noqa: E501
        return self._cached("controls", lambda: json.loads(self._core.controls_json))

    @property
    def truth(self) -> dict[str, Any]:
        """The truth manifest: what the study was built to show, and what its data realized.

        A dict of JSON-compatible values. Every manifest has the request (`archetype`,
        `studyid`, `n`, `seed`, `controls`) and the arms in order (`arms`); the other
        sections depend on the archetype and hold targets (what the controls asked for) next
        to the realized values (what a standard analysis of this study's ADaM gives), such as
        the mean change from baseline at the end of treatment, Kaplan-Meier medians and
        stratified Cox hazard ratios, or response rates. `write()` and `write_submission()`
        save it as `truth.json`.

        >>> import trialforge as tf
        >>> trial = tf.study("parallel_continuous", n=200, seed=1, indication="hypertension", weeks=12)
        >>> list(trial.truth)
        ['archetype', 'arms', 'controls', 'discontinuation', 'endpoint', 'n', 'seed', 'studyid', 'teae_incidence']
        >>> trial.truth["endpoint"]["target_difference"]
        -8.0
        """  # noqa: E501
        return self._cached("truth", lambda: json.loads(self._core.truth_json))

    @property
    def defects(self) -> Any:
        """The injected defects, as a DataFrame of the `output` type, or `None` without defects.

        One row per defect, with the columns:

        - `id`: the defect type, e.g. `out_of_ct`.
        - `dataset` and `row`: the record it was injected into (rows are zero-based, as in
          `sdtm` and `adam`).
        - `paired_row`: another record the defect involves, which a validator may name
          instead (the record a `duplicate_seq` copied), or null.
        - `usubjid` and `key`: the subject, and the record's key before the edit, e.g.
          `AESEQ=3`.
        - `variable`, `original` and `injected`: the value that was changed, and how.
        - `expected`: the findings trialforge's checks should produce, as
          `RULE DATASET row`, separated by `; `.

        The datasets in `sdtm` and `adam` are the defective copies; `truth` describes the
        clean study.

        >>> import trialforge as tf
        >>> dirty = tf.study("parallel_continuous", n=200, seed=1, defects={"out_of_ct": 3})
        >>> list(dirty.defects.columns)
        ['id', 'dataset', 'row', 'paired_row', 'usubjid', 'key', 'variable', 'original', 'injected', 'expected']
        """  # noqa: E501

        def make() -> Any:
            table = self._core.defects_table()
            return None if table is None else self._convert(table)

        return self._cached("defects", make)

    @property
    def sdtm(self) -> dict[str, Any]:
        """The trial design and SDTM datasets, by name, as DataFrames of the `output` type.

        The trial design datasets (TA, TE, TV, TI, TS) come first, then the subject domains
        (DM, SE, SV, EX, DS and those of the archetype).

        >>> import trialforge as tf
        >>> trial = tf.study("parallel_continuous", n=50, seed=1)
        >>> list(trial.sdtm)
        ['TA', 'TE', 'TV', 'TI', 'TS', 'DM', 'SE', 'SV', 'EX', 'CM', 'AE', 'DS', 'MH', 'LB', 'VS', 'SUPPAE']
        >>> len(trial.sdtm["DM"])
        50
        """  # noqa: E501
        return self._cached(
            "sdtm", lambda: {name: self._frame(name) for name in self._core.sdtm_names}
        )

    @property
    def adam(self) -> dict[str, Any]:
        """The ADaM datasets, by name, as DataFrames of the `output` type.

        ADSL comes first. Every ADaM dataset is derived from the SDTM datasets alone.

        >>> import trialforge as tf
        >>> trial = tf.study("oncology_solid_tumor", n=50, seed=1)
        >>> list(trial.adam)
        ['ADSL', 'ADRS', 'ADTTE', 'ADTR']
        """
        return self._cached(
            "adam", lambda: {name: self._frame(name) for name in self._core.adam_names}
        )

    def write(self, path: str | PathLike[str], formats: Sequence[str] = ("xpt",)) -> list[Path]:
        """Write every dataset to files, plus the truth manifest.

        Writes one file per dataset and format to the directory `path`, named by the
        lower-case dataset name and the format's extension (`dm.xpt`, `dm.csv`, `dm.json`,
        `dm.parquet`), together with `truth.json` and, for a study with defects,
        `defects.json`. Files are deterministic: the same study writes the same bytes.

        Parameters
        ----------
        path
            The output directory. It is created if needed, and existing files of the same
            names are overwritten.
        formats
            One or more of the formats `formats()` lists: `"xpt"`, `"csv"`,
            `"dataset-json"` or `"parquet"`. A single string is one format.

        Returns
        -------
        list[Path]
            The paths written.

        Raises
        ------
        ValueError
            If a format is unknown.
        OSError
            If a file can't be written.

        See Also
        --------
        write_submission : Write SDTM and ADaM folders with define.xml.

        Examples
        --------
        >>> import tempfile
        >>> import trialforge as tf
        >>> trial = tf.study("parallel_continuous", n=50, seed=1)
        >>> paths = trial.write(tempfile.mkdtemp(), formats=["xpt", "csv"])
        >>> [p.name for p in paths[:4]]
        ['ta.xpt', 'ta.csv', 'te.xpt', 'te.csv']
        >>> paths[-1].name
        'truth.json'
        """
        if isinstance(formats, str):
            formats = [formats]
        return [Path(p) for p in self._core.write(Path(path), list(formats))]

    def write_submission(self, path: str | PathLike[str]) -> list[Path]:
        """Write a submission-style package.

        Creates `sdtm/` and `adam/` in the directory `path`, each with one XPT (version 5)
        file per dataset and a Define-XML 2.1 `define.xml` describing them, plus `truth.json`
        (and `defects.json` for a study with defects) in `path` itself.

        Parameters
        ----------
        path
            The output directory. It is created if needed, and existing files of the same
            names are overwritten.

        Returns
        -------
        list[Path]
            The paths written.

        Raises
        ------
        OSError
            If a file can't be written.

        Examples
        --------
        >>> import tempfile
        >>> from pathlib import Path
        >>> import trialforge as tf
        >>> trial = tf.study("parallel_continuous", n=50, seed=1)
        >>> out = Path(tempfile.mkdtemp())
        >>> paths = trial.write_submission(out)
        >>> [p.relative_to(out).as_posix() for p in paths if p.suffix != ".xpt"]
        ['sdtm/define.xml', 'adam/define.xml', 'truth.json']
        """
        return [Path(p) for p in self._core.write_submission(Path(path))]

    def define_xml(self, standard: Literal["sdtm", "adam"] = "sdtm") -> str:
        """The Define-XML 2.1 document for the SDTM or ADaM datasets.

        The document describes every dataset and variable (label, type, length, origin), the
        codelists with the terms the study uses, value-level metadata for findings tests and
        analysis parameters and, for ADaM, how each variable is derived. It is the
        `define.xml` that `write_submission()` writes.

        Parameters
        ----------
        standard
            `"sdtm"` (the trial design and SDTM datasets) or `"adam"`.

        Returns
        -------
        str
            The XML document.

        Raises
        ------
        ValueError
            If `standard` is neither `"sdtm"` nor `"adam"`.

        Examples
        --------
        >>> import trialforge as tf
        >>> xml = tf.study("parallel_continuous", n=50, seed=1).define_xml("adam")
        >>> xml.splitlines()[0]
        '<?xml version="1.0" encoding="UTF-8"?>'
        """
        return self._core.define_xml(standard)

    def verify(self, level: Literal[2, 3] = 2) -> Report:
        """Verify the study.

        Parameters
        ----------
        level
            - `2` (the default) re-runs trialforge's structural (L0), referential (L1) and
              derivational (L2) checks on the study's datasets. A study without defects always
              passes them (generation guarantees it); with seeded defects, the report lists
              what the checks find in the defective datasets, for `score()`.
            - `3` also runs Pointblank's SDTMIG 3.4 conformance rules on the written SDTM
              package, including its define.xml checks. It needs
              `pip install 'trialforge[verify]'`.

        Returns
        -------
        Report
            Whether the study passed, and its findings.

        Raises
        ------
        ValueError
            If `level` is neither 2 nor 3.
        ImportError
            At level 3, if Pointblank isn't installed.

        Examples
        --------
        >>> import trialforge as tf
        >>> tf.study("parallel_continuous", n=50, seed=1).verify().passed
        True
        >>> dirty = tf.study("parallel_continuous", n=200, seed=1, defects={"out_of_ct": 3})
        >>> report = dirty.verify()
        >>> report.passed, len(report.findings)
        (False, 3)
        >>> report.findings[0].rule
        'L0.CT'
        """
        if level not in (2, 3):
            raise ValueError("`level` must be 2 or 3")
        passed, raw = self._core.verify()
        findings = [Finding(*f) for f in raw]
        waived: list[Finding] = []
        if level == 3:
            from trialforge._conformance import pointblank_findings

            failed, waived = pointblank_findings(self)
            findings.extend(failed)
            passed = passed and not failed
        return Report(passed, findings, waived)

    def score(self, report: Report) -> Score:
        """Score a report's findings against the study's injected defects.

        Matches each finding to the defects that explain it, for recall (the share of defects
        found) and precision (the share of findings that are real), overall and per defect
        type. The findings can come from `verify()`, at either level, or from another
        validator, as a `Report` of `Finding`s with zero-based rows.

        Parameters
        ----------
        report
            The findings to score. Waived findings are not scored.

        Returns
        -------
        Score
            Recall and precision, the defects missed and the findings unmatched.

        Raises
        ------
        ValueError
            If the study has no defects.

        Examples
        --------
        >>> import trialforge as tf
        >>> dirty = tf.study("parallel_continuous", n=200, seed=1, defects="light")
        >>> score = dirty.score(dirty.verify())
        >>> score.recall, score.precision
        (1.0, 1.0)
        """
        findings = report.findings
        raw = json.loads(self._core.score([(f.level, f.rule, f.dataset, f.row) for f in findings]))
        return Score(
            injected=raw["injected"],
            detected=raw["detected"],
            recall=raw["recall"],
            findings=raw["findings"],
            matched=raw["matched"],
            precision=raw["precision"],
            by_type={k: DefectScore(**v) for k, v in raw["by_type"].items()},
            missed=raw["missed"],
            unmatched=[findings[i] for i in raw["unmatched"]],
        )

    def _frame(self, name: str) -> Any:
        return self._convert(self._core.table(name))

    def _convert(self, table: _core.ArrowTable) -> Any:
        if self._output == "polars":
            pl = _backend("polars", "pl")
            return pl.DataFrame(table)
        pa = _backend("pyarrow", "pa" if self._output == "pyarrow" else "pd")
        arrow = pa.table(table)
        if self._output == "pyarrow":
            return arrow
        _backend("pandas", "pd")
        return arrow.to_pandas()

    def __repr__(self) -> str:
        names = ", ".join([*self._core.sdtm_names, *self._core.adam_names])
        return (
            f"<trialforge.Study {self.studyid} ({self.archetype}), seed={self.seed}, "
            f"cutoff={self.cutoff}: {names}>"
        )
