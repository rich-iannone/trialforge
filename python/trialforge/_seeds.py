"""Seeds and the core version, documented here for the reference (the work is in Rust)."""

from __future__ import annotations

from collections.abc import Sequence

from trialforge import _core

MAX_SEED: int = _core.MAX_SEED
"""The largest root seed, `2**53`.

Seeds are whole numbers in `[0, MAX_SEED]`: the integers R's doubles hold exactly, so a seed means
the same in Python and R.
"""


def stream_seed(seed: int, path: Sequence[int | str]) -> int:
    """The seed of the random stream at `path` under a root seed.

    Every random draw in a study comes from a stream addressed by a path under the study's root
    seed (the subject, then the process, then the analyte or term), so that each independent
    process has its own stream. Python and R give the same seed for the same path.

    Parameters
    ----------
    seed
        The root seed.
    path
        The stream's path. Integers are index components and strings are label components,
        e.g. `["subject", 17, "ae"]`.

    Returns
    -------
    int
        The stream's 64-bit seed.
    """
    return _core.stream_seed(seed, path)


def core_version() -> str:
    """The version of the Rust core the package was built with."""
    return _core.core_version()
