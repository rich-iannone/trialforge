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
    Generate a 12-week hypertension study with 200 subjects. The `indication` and `weeks`
    controls choose the preset and the treatment duration, and fixing the seed makes the
    study reproducible:

    ```{python}
    import trialforge as tf

    trial = tf.study("parallel_continuous", n=200, seed=1, indication="hypertension", weeks=12)
    trial
    ```

    The study's representation gives its `STUDYID`, seed and data cutoff, then lists its
    datasets: the trial design datasets (TA to TS), the SDTM domains and the ADaM datasets.
    Each dataset is a DataFrame of the `output` type, Polars by default. Here are the first
    subjects in ADSL:

    ```{python}
    trial.adam["ADSL"].select("USUBJID", "TRT01P", "AGE", "SEX", "TRTSDT", "TRTEDT").head()
    ```

    `Study.controls` records every control the study used: the two passed here, and the
    defaults of the rest.

    ```{python}
    trial.controls
    ```

    The truth manifest reports what the study was built to show next to what its data gives.
    For this study, the endpoint is the change from baseline in systolic blood pressure at the
    end of treatment:

    ```{python}
    trial.truth["endpoint"]
    ```

    The target difference from placebo comes from the preset's changes in each arm. The
    realized difference, computed from the study's ADVS, differs from it by sampling
    variation, as a real study's estimate would.

    Controls are checked before anything is generated, so a misspelled control is an error
    that lists the valid names:

    ```{python}
    try:
        tf.study("parallel_continuous", wekks=12)
    except ValueError as e:
        print(e)
    ```

    To test a validator, ask for seeded defects. The datasets then carry the defects, and
    `Study.defects` lists each one:

    ```{python}
    dirty = tf.study("parallel_continuous", n=200, seed=1, defects={"out_of_ct": 3})
    dirty.defects
    ```

    Each row gives the record a defect is in, the value it changed and the finding
    trialforge's checks should report for it (in `expected`).
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
    A spec can be a dict with the same structure as the JSON or YAML file. This one asks for
    an obesity study with 30 subjects:

    ```{python}
    import trialforge as tf

    spec = {
        "study": {"archetype": "parallel_continuous", "seed": 1},
        "population": {"n": 30},
        "controls": {"indication": "obesity"},
    }
    trial = tf.study_from_spec(spec)
    trial
    ```

    It gives the same study as the equivalent call to `study()`:

    ```{python}
    same = tf.study("parallel_continuous", n=30, seed=1, indication="obesity")
    same.truth == trial.truth
    ```

    The spec is more often a file. Its text can be passed directly too, here in YAML:

    ```{python}
    text = '''
    study: {archetype: parallel_continuous, seed: 1}
    population: {n: 30}
    controls: {indication: obesity}
    '''
    tf.study_from_spec(text).truth == trial.truth
    ```

    A file's path works the same way, so a spec kept under version control reproduces its
    study wherever it is run.
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
    List the archetypes with their descriptions:

    ```{python}
    import trialforge as tf

    tf.archetypes()
    ```

    Each description names the study's design and the domains particular to it. Any of the
    names can be passed to `study()`, and `controls()` gives what it accepts:

    ```{python}
    trial = tf.study("vaccine", n=50, seed=1)
    trial.truth["arms"]
    ```
    """
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
    List the controls of the `parallel_continuous` archetype:

    ```{python}
    import trialforge as tf

    tf.controls("parallel_continuous")
    ```

    Each `Control` gives a name, what the control sets and its default. Their names alone are
    a quick reminder of what can be changed:

    ```{python}
    [c.name for c in tf.controls("parallel_continuous")]
    ```

    Set a control by passing it to `study()` as a keyword argument. This study treats for 12
    weeks instead of the default 24, which also shortens the visit schedule:

    ```{python}
    trial = tf.study("parallel_continuous", n=50, seed=1, weeks=12)
    trial.sdtm["TV"].select("VISITNUM", "VISIT", "VISITDY")
    ```
    """
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
    List the supported formats:

    ```{python}
    import trialforge as tf

    tf.formats()
    ```

    Pass one or more of them to `Study.write()`. This writes each dataset as both XPT and
    Dataset-JSON, to a temporary directory:

    ```{python}
    import tempfile

    trial = tf.study("parallel_continuous", n=50, seed=1)
    paths = trial.write(tempfile.mkdtemp(), formats=["xpt", "dataset-json"])
    [p.name for p in paths if p.stem == "dm"]
    ```
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
    A study generated without defects passes every check, so its report has no findings:

    ```{python}
    import trialforge as tf

    tf.study("parallel_continuous", n=50, seed=1).verify()
    ```

    With seeded defects, the report lists what the checks find in the defective datasets:

    ```{python}
    dirty = tf.study("parallel_continuous", n=200, seed=1, defects={"out_of_ct": 3})
    report = dirty.verify()
    report.passed, len(report.findings)
    ```

    Each finding names the rule that failed and the record it failed on:

    ```{python}
    report.findings[0]
    ```
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
    Generate a study with the `light` profile of defects, and score trialforge's own checks
    against them:

    ```{python}
    import trialforge as tf

    dirty = tf.study("parallel_continuous", n=200, seed=1, defects="light")
    score = dirty.score(dirty.verify())
    score.recall, score.precision
    ```

    The checks find every defect and report nothing else, so recall and precision are both
    1. A validator that missed some defects would have a lower recall, with the defects it
    missed in `missed`; one that flagged clean records would have a lower precision, with
    those findings in `unmatched`. `by_type` breaks the score down by defect type:

    ```{python}
    score.by_type
    ```
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
    Generate a study, then look at what it holds:

    ```{python}
    import trialforge as tf

    trial = tf.study("parallel_continuous", n=200, seed=1, indication="hypertension", weeks=12)
    trial
    ```

    The representation lists the datasets in `sdtm` (from TA to SUPPAE) and in `adam` (from
    ADSL). Each is a DataFrame, here the vital signs domain:

    ```{python}
    trial.sdtm["VS"].head()
    ```

    The request that made the study is kept with it, so it can be generated again:

    ```{python}
    trial.studyid, trial.archetype, trial.seed
    ```
    """

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

        The seed passed to `study()`, or the one drawn for it when none was.

        Examples
        --------
        Generate a study without a seed. A fresh one is drawn, and recorded:

        ```{python}
        import trialforge as tf

        trial = tf.study("vaccine", n=50)
        trial.seed
        ```

        Passing that seed back, with the same archetype, `n` and controls, gives the same
        study again:

        ```{python}
        again = tf.study("vaccine", n=50, seed=trial.seed)
        again.truth == trial.truth
        ```
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

        Examples
        --------
        A 12-week study's cutoff is the date of its last subject's last visit:

        ```{python}
        import trialforge as tf

        trial = tf.study("parallel_continuous", n=200, seed=1, weeks=12)
        trial.cutoff
        ```

        No date in the study's data is later. In ADSL, the last date of contact with a
        subject is at the cutoff at the latest:

        ```{python}
        trial.adam["ADSL"]["LSTALVDT"].max()
        ```
        """
        return _EPOCH + datetime.timedelta(days=self._core.cutoff_days)

    @property
    def controls(self) -> dict[str, Any]:
        """Every control the study used: the ones you passed, and the defaults of the rest.

        Together with `archetype`, `n`, `seed` and `studyid`, they reproduce the study:
        `study(trial.archetype, n=..., seed=trial.seed, **trial.controls)`.

        Examples
        --------
        Generate a study with two controls set:

        ```{python}
        import trialforge as tf

        trial = tf.study("parallel_continuous", n=200, seed=1, indication="hypertension", weeks=12)
        trial.controls
        ```

        Besides `indication` and `weeks`, the dict has every other control at its default.
        Passing it back to `study()` reproduces the study:

        ```{python}
        again = tf.study(trial.archetype, n=200, seed=trial.seed, **trial.controls)
        again.truth == trial.truth
        ```
        """
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

        Examples
        --------
        Generate a hypertension study and list its manifest's sections:

        ```{python}
        import trialforge as tf

        trial = tf.study("parallel_continuous", n=200, seed=1, indication="hypertension", weeks=12)
        list(trial.truth)
        ```

        Next to the request, this archetype reports its endpoint, the incidence of
        treatment-emergent adverse events and discontinuation. The endpoint section holds the
        target mean change in systolic blood pressure in each arm and the difference between
        them, with the realized values from the study's ADVS:

        ```{python}
        trial.truth["endpoint"]
        ```

        An analysis of this study's ADaM should reproduce the realized values; the target is
        what a very large study would show.
        """
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

        Examples
        --------
        Inject three values outside their controlled terminology:

        ```{python}
        import trialforge as tf

        dirty = tf.study("parallel_continuous", n=200, seed=1, defects={"out_of_ct": 3})
        dirty.defects
        ```

        Each defect changed a specimen type in LB (`LBSPEC`) to a value its codelist doesn't
        have. The defective records are in `sdtm`, at the rows listed:

        ```{python}
        rows = dirty.defects["row"].to_list()
        dirty.sdtm["LB"][rows].select("USUBJID", "LBTESTCD", "LBSPEC")
        ```

        A study generated without defects has none to list:

        ```{python}
        tf.study("parallel_continuous", n=200, seed=1).defects is None
        ```
        """

        def make() -> Any:
            table = self._core.defects_table()
            return None if table is None else self._convert(table)

        return self._cached("defects", make)

    @property
    def sdtm(self) -> dict[str, Any]:
        """The trial design and SDTM datasets, by name, as DataFrames of the `output` type.

        The trial design datasets (TA, TE, TV, TI, TS) come first, then the subject domains
        (DM, SE, SV, EX, DS and those of the archetype).

        Examples
        --------
        List a study's SDTM datasets:

        ```{python}
        import trialforge as tf

        trial = tf.study("parallel_continuous", n=50, seed=1)
        list(trial.sdtm)
        ```

        Besides the domains every study has, this archetype has the safety domains (CM, AE,
        MH, LB, VS) and the supplemental qualifiers of AE (SUPPAE). DM has one row per
        subject:

        ```{python}
        trial.sdtm["DM"].head()
        ```
        """
        return self._cached(
            "sdtm", lambda: {name: self._frame(name) for name in self._core.sdtm_names}
        )

    @property
    def adam(self) -> dict[str, Any]:
        """The ADaM datasets, by name, as DataFrames of the `output` type.

        ADSL comes first. Every ADaM dataset is derived from the SDTM datasets alone.

        Examples
        --------
        List the ADaM datasets of an oncology study:

        ```{python}
        import trialforge as tf

        trial = tf.study("oncology_solid_tumor", n=50, seed=1)
        list(trial.adam)
        ```

        After ADSL come the response (ADRS), time-to-event (ADTTE) and tumor results (ADTR)
        datasets. ADTTE has one row per subject and parameter: progression-free and overall
        survival for every subject, and the duration of response for responders.

        ```{python}
        trial.adam["ADTTE"].select("USUBJID", "PARAMCD", "AVAL", "CNSR", "EVNTDESC").head()
        ```
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
        Write a study's datasets as XPT and CSV files, to a temporary directory:

        ```{python}
        import tempfile

        import trialforge as tf

        trial = tf.study("parallel_continuous", n=50, seed=1)
        paths = trial.write(tempfile.mkdtemp(), formats=["xpt", "csv"])
        [p.name for p in paths]
        ```

        There are two files per dataset, one per format, in the order of `sdtm` then `adam`,
        and `truth.json` last. The files can be read back with any reader of their format,
        such as Polars for the CSV files:

        ```{python}
        import polars as pl

        dm = next(p for p in paths if p.name == "dm.csv")
        pl.read_csv(dm).head()
        ```
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
        Write a submission package to a temporary directory:

        ```{python}
        import tempfile
        from pathlib import Path

        import trialforge as tf

        trial = tf.study("parallel_continuous", n=50, seed=1)
        out = Path(tempfile.mkdtemp())
        paths = trial.write_submission(out)
        [p.relative_to(out).as_posix() for p in paths]
        ```

        Each folder has its datasets as XPT files, followed by the `define.xml` that describes
        them, and `truth.json` is at the top. A conformance checker that reads a submission's
        SDTM or ADaM folder can be pointed at `sdtm/` or `adam/`.
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
        Get the ADaM define.xml of a study, and look at its start:

        ```{python}
        import trialforge as tf

        xml = tf.study("parallel_continuous", n=50, seed=1).define_xml("adam")
        print(xml[:600])
        ```

        The document is plain text, so it can be parsed with any XML library. Here, the
        datasets it describes, from its `ItemGroupDef` elements:

        ```{python}
        import xml.etree.ElementTree as ET

        odm = "{http://www.cdisc.org/ns/odm/v1.3}"
        root = ET.fromstring(xml)
        [group.get("Name") for group in root.iter(f"{odm}ItemGroupDef")]
        ```
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
        A study generated without defects passes every check:

        ```{python}
        import trialforge as tf

        tf.study("parallel_continuous", n=50, seed=1).verify()
        ```

        Inject three defects, values outside their controlled terminology, and verify again:

        ```{python}
        dirty = tf.study("parallel_continuous", n=200, seed=1, defects={"out_of_ct": 3})
        report = dirty.verify()
        report.passed, len(report.findings)
        ```

        The study now fails, with one finding per defect. Each names the level and rule that
        caught it, and the dataset, variable and row it is in:

        ```{python}
        for finding in report.findings:
            print(finding)
        ```

        Pass the report to `score()` to compare the findings with the defects injected.
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
        Generate a study with the `light` profile of defects, and score trialforge's own
        checks against them:

        ```{python}
        import trialforge as tf

        dirty = tf.study("parallel_continuous", n=200, seed=1, defects="light")
        score = dirty.score(dirty.verify())
        score.recall, score.precision
        ```

        The checks find every defect, and every finding is a defect. A validator's findings
        are scored the same way. This one only checks AE, so it finds just the defect injected
        there:

        ```{python}
        from trialforge import Report

        found = [f for f in dirty.verify().findings if f.dataset == "AE"]
        partial = dirty.score(Report(passed=False, findings=found))
        partial.recall, partial.precision
        ```

        `missed` gives the rows in `defects` of the defects it didn't find:

        ```{python}
        dirty.defects[partial.missed].select("id", "dataset", "row", "variable")
        ```
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
