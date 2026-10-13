"""The `Study` object and `study()` constructor."""

from __future__ import annotations

import datetime
from collections.abc import Sequence
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
) -> Study:
    """Generate a mock study.

    Parameters
    ----------
    archetype
        The study type, e.g. `"oncology_solid_tumor"`. See `archetypes()`.
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
    """
    if output not in _OUTPUTS:
        raise ValueError(f"`output` must be one of {', '.join(map(repr, _OUTPUTS))}")
    return Study(_core.study(archetype, n, seed, studyid), output)


def archetypes() -> dict[str, str]:
    """The available study archetypes, mapped to one-line descriptions."""
    return dict(_core.archetypes())


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
        """The data cutoff date."""
        return _EPOCH + datetime.timedelta(days=self._core.cutoff_days)

    @cached_property
    def sdtm(self) -> dict[str, Any]:
        """SDTM datasets by name, as DataFrames of the chosen `output` type."""
        return {name: self._frame(name) for name in self._core.sdtm_names}

    @cached_property
    def adam(self) -> dict[str, Any]:
        """ADaM datasets by name, as DataFrames of the chosen `output` type."""
        return {name: self._frame(name) for name in self._core.adam_names}

    def write(self, path: str | PathLike[str], formats: Sequence[str] = ("xpt",)) -> list[Path]:
        """Write every dataset to the directory `path` (created if needed).

        Files are named by lower-case dataset name, e.g. `dm.xpt`. Returns the paths written.
        """
        if isinstance(formats, str):
            formats = [formats]
        return [Path(p) for p in self._core.write(Path(path), list(formats))]

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
