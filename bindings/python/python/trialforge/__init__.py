"""trialforge: synthetic CDISC clinical-trial data from declarative study specs.

This is a pre-alpha scaffold.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from trialforge._core import core_version, stream_seed

try:
    __version__ = version("trialforge")
except PackageNotFoundError:  # pragma: no cover - running from a source tree without install
    __version__ = core_version()

__all__ = ["__version__", "core_version", "stream_seed"]
