"""The `Study` object and `study()` constructor."""

from __future__ import annotations

import datetime
import json
from collections.abc import Sequence
from dataclasses import dataclass
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
    **controls
        Archetype controls, e.g. `indication="hypertension", weeks=12`. See `controls()`.
        Unknown names are an error.
    """
    if output not in _OUTPUTS:
        raise ValueError(f"`output` must be one of {', '.join(map(repr, _OUTPUTS))}")
    payload = json.dumps(controls) if controls else None
    return Study(_core.study(archetype, n, seed, studyid, payload), output)


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
    """The result of `Study.verify()`."""

    passed: bool
    findings: list[Finding]


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
    def sdtm(self) -> dict[str, Any]:
        """SDTM datasets by name, as DataFrames of the chosen `output` type."""
        return {name: self._frame(name) for name in self._core.sdtm_names}

    @cached_property
    def adam(self) -> dict[str, Any]:
        """ADaM datasets by name, as DataFrames of the chosen `output` type."""
        return {name: self._frame(name) for name in self._core.adam_names}

    def write(self, path: str | PathLike[str], formats: Sequence[str] = ("xpt",)) -> list[Path]:
        """Write every dataset to the directory `path` (created if needed), plus `truth.json`.

        Files are named by lower-case dataset name, e.g. `dm.xpt`. See `formats()` for the
        supported formats. Returns the paths written.
        """
        if isinstance(formats, str):
            formats = [formats]
        return [Path(p) for p in self._core.write(Path(path), list(formats))]

    def verify(self) -> Report:
        """Re-run verification levels 0–2 (structural, referential, derivational).

        Generation already guarantees these pass; this is for after you modify the data.
        """
        passed, findings = self._core.verify()
        return Report(passed, [Finding(*f) for f in findings])

    def _frame(self, name: str) -> Any:
        table = self._core.table(name)
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
