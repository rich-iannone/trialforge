"""The `Study` object and `study()` constructor."""

from __future__ import annotations

import datetime
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from functools import cached_property
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

    Parameters
    ----------
    archetype
        The study type, e.g. `"parallel_continuous"`. See `archetypes()`.
    n
        Number of randomized subjects.
    seed
        Root seed, in `[0, 2**53]`. With `None`, a fresh seed is drawn and recorded on the
        result as `Study.seed`, so the study can always be reproduced.
    studyid
        `STUDYID`. Each archetype has a default.
    output
        DataFrame type for `Study.sdtm` and `Study.adam`: `"polars"`, `"pandas"` or
        `"pyarrow"`.
    defects
        Data-quality defects to inject after the clean study is generated and verified:
        `"none"`; a profile, `"light"` (about five defects of the common kinds),
        `"realistic"` (rates per eligible record) or `"adversarial"` (three of every type,
        each in its subtlest form); or a mapping of defect id to a count (an integer) or a
        rate per eligible record (a float in `[0, 1]`), e.g.
        `{"ae_end_before_start": 3, "out_of_ct": 0.02}`. See `Study.defects` and
        `Study.score()`.
    **controls
        Archetype controls, e.g. `indication="hypertension", weeks=12`. See `controls()`.
        Unknown names are an error.
    """
    if output not in _OUTPUTS:
        raise ValueError(f"`output` must be one of {', '.join(map(repr, _OUTPUTS))}")
    payload = json.dumps(controls) if controls else None
    requested = None if defects is None else json.dumps(defects)
    return Study(_core.study(archetype, n, seed, studyid, payload, requested), output)


def study_from_spec(
    spec: str | PathLike[str] | dict[str, Any], *, output: Output = "polars"
) -> Study:
    """Generate a study from a study spec (see `spec/study-spec.schema.json`).

    `spec` is a path to a JSON or YAML file, the JSON/YAML text itself, or a dict. The spec's
    `outputs` section is for the command-line tool; write files with `Study.write()` or
    `Study.write_submission()`.
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


def archetypes() -> dict[str, str]:
    """The available study archetypes, mapped to one-line descriptions."""
    return dict(_core.archetypes())


@dataclass(frozen=True)
class Control:
    """One archetype control."""

    name: str
    description: str
    default: Any


def controls(archetype: str) -> list[Control]:
    """The controls an archetype accepts, with descriptions and defaults."""
    return [
        Control(c["name"], c["description"], c["default"])
        for c in json.loads(_core.controls(archetype))
    ]


def formats() -> list[str]:
    """Output formats `Study.write()` supports."""
    return _core.formats()


@dataclass(frozen=True)
class Finding:
    """One failed verification check."""

    level: str
    rule: str
    dataset: str
    variable: str | None
    row: int | None
    message: str


@dataclass(frozen=True)
class Report:
    """The result of `Study.verify()`.

    `waived` lists level-3 findings caused by known defects in Pointblank's rule catalog (each
    message names the defect); they do not affect `passed`.
    """

    passed: bool
    findings: list[Finding]
    waived: list[Finding] = field(default_factory=list[Finding])


@dataclass(frozen=True)
class DefectScore:
    """Recall and precision for one defect type (see `Score`)."""

    injected: int
    detected: int
    recall: float | None
    findings: int
    false_positives: int
    precision: float | None


@dataclass(frozen=True)
class Score:
    """How well a report's findings recover a study's injected defects (`Study.score()`).

    A finding matches a defect when it names the dataset and row the defect expects a finding
    at and, for trialforge's own levels (L0–L2), the expected rule; Pointblank (L3) findings
    match on dataset and row only (or a duplicate's paired row). `recall` is
    `detected / injected` and `precision` `matched / findings` (`None` when undefined).
    `missed` holds row indices of `Study.defects` that no finding matched; `unmatched` the
    findings no defect explains.
    In `by_type`, `false_positives` counts unmatched L0–L2 findings from a rule that catches
    that type.
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
    """A generated mock study: SDTM and ADaM datasets plus how they were made."""

    def __init__(self, core: _core.Study, output: Output) -> None:
        self._core = core
        self._output: Output = output

    @property
    def studyid(self) -> str:
        """`STUDYID`."""
        return self._core.studyid

    @property
    def archetype(self) -> str:
        """The archetype the study was generated from."""
        return self._core.archetype

    @property
    def seed(self) -> int:
        """The root seed; pass it to `study()` to reproduce this study."""
        return self._core.seed

    @property
    def cutoff(self) -> datetime.date:
        """The data cutoff (event-driven studies) or last contact date (fixed-duration)."""
        return _EPOCH + datetime.timedelta(days=self._core.cutoff_days)

    @cached_property
    def controls(self) -> dict[str, Any]:
        """The full control set used (defaults merged with the ones you passed)."""
        return json.loads(self._core.controls_json)

    @cached_property
    def truth(self) -> dict[str, Any]:
        """The truth manifest: targets and the values realised in the generated ADaM."""
        return json.loads(self._core.truth_json)

    @cached_property
    def defects(self) -> Any:
        """The injected defects as a DataFrame of the chosen `output` type, or `None`.

        One row per defect: `id`, `dataset`, `row` (zero-based), `paired_row` (another record
        the defect involves, which a validator may name instead: the record a `duplicate_seq`
        copied), `usubjid`, `key` (the record's key before the edit, e.g. `AESEQ=3`),
        `variable`, `original`, `injected` and `expected` (the findings it should produce, as
        `RULE DATASET row`, `; `-separated).
        """
        table = self._core.defects_table()
        return None if table is None else self._convert(table)

    @cached_property
    def sdtm(self) -> dict[str, Any]:
        """SDTM datasets by name, as DataFrames of the chosen `output` type."""
        return {name: self._frame(name) for name in self._core.sdtm_names}

    @cached_property
    def adam(self) -> dict[str, Any]:
        """ADaM datasets by name, as DataFrames of the chosen `output` type."""
        return {name: self._frame(name) for name in self._core.adam_names}

    def write(self, path: str | PathLike[str], formats: Sequence[str] = ("xpt",)) -> list[Path]:
        """Write every dataset to the directory `path` (created if needed), plus `truth.json`.

        A study with defects also gets `defects.json`. Files are named by lower-case dataset
        name, e.g. `dm.xpt`. See `formats()` for the supported formats. Returns the paths
        written.
        """
        if isinstance(formats, str):
            formats = [formats]
        return [Path(p) for p in self._core.write(Path(path), list(formats))]

    def write_submission(self, path: str | PathLike[str]) -> list[Path]:
        """Write a submission-style package to the directory `path`.

        Creates `sdtm/` and `adam/`, each with one XPT file per dataset and a Define-XML 2.1
        `define.xml`, plus `truth.json` (and `defects.json`). Returns the paths written.
        """
        return [Path(p) for p in self._core.write_submission(Path(path))]

    def define_xml(self, standard: Literal["sdtm", "adam"] = "sdtm") -> str:
        """The Define-XML 2.1 document for the SDTM or ADaM datasets."""
        return self._core.define_xml(standard)

    def verify(self, level: Literal[2, 3] = 2) -> Report:
        """Verify the study.

        - Level 2 (default) re-runs trialforge's structural, referential and derivational
          checks. Generation already guarantees these pass; this is for after you modify data.
        - Level 3 also runs Pointblank's SDTMIG 3.4 conformance rules on the written SDTM
          package, including its define.xml checks (needs `pip install 'trialforge[verify]'`).
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
        """Score `report`'s findings (from `verify()`, at either level) against the defects.

        Waived findings are not scored. Raises `ValueError` if the study has no defects.
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
            import polars as pl

            return pl.DataFrame(table)
        import pyarrow as pa

        arrow = pa.table(table)
        return arrow if self._output == "pyarrow" else arrow.to_pandas()

    def __repr__(self) -> str:
        names = ", ".join([*self._core.sdtm_names, *self._core.adam_names])
        return (
            f"<trialforge.Study {self.studyid} ({self.archetype}), seed={self.seed}, "
            f"cutoff={self.cutoff}: {names}>"
        )
