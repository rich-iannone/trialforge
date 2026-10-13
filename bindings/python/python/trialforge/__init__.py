"""trialforge: synthetic CDISC clinical-trial data from declarative study specs.

This is a pre-alpha preview: `study("oncology_solid_tumor", ...)` generates DM, EX, DS and
ADSL.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from trialforge._core import MAX_SEED, core_version, stream_seed
from trialforge._study import Output, Study, archetypes, study

try:
    __version__ = version("trialforge")
except PackageNotFoundError:  # pragma: no cover - running from a source tree without install
    __version__ = core_version()

__all__ = [
    "MAX_SEED",
    "Output",
    "Study",
    "__version__",
    "archetypes",
    "core_version",
    "stream_seed",
    "study",
]
