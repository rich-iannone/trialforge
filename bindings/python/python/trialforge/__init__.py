"""trialforge: synthetic CDISC clinical-trial data from declarative study specs.

Pre-alpha. `study("parallel_continuous", ...)` and `study("oncology_solid_tumor", ...)`
generate a verified SDTM + ADaM study; see `archetypes()` and `controls()`.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from trialforge._core import MAX_SEED, core_version, stream_seed
from trialforge._study import (
    Control,
    Finding,
    Output,
    Report,
    Study,
    archetypes,
    controls,
    formats,
    study,
    study_from_spec,
)

try:
    __version__ = version("trialforge")
except PackageNotFoundError:  # pragma: no cover - running from a source tree without install
    __version__ = core_version()

__all__ = [
    "MAX_SEED",
    "Control",
    "Finding",
    "Output",
    "Report",
    "Study",
    "__version__",
    "archetypes",
    "controls",
    "core_version",
    "formats",
    "stream_seed",
    "study",
    "study_from_spec",
]
