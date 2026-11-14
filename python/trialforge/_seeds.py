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

    Raises
    ------
    TypeError
        If a path component is neither an integer nor a string.

    Examples
    --------
    The seed of the adverse-event stream of subject 17, in a study with root seed 1:

    ```{python}
    import trialforge as tf

    tf.stream_seed(1, ["subject", 17, "ae"])
    ```

    The R package gives the same number for the same path (as a string, since R can't hold a
    64-bit integer). Each component of the path changes the seed, so subject 18's stream is
    independent of subject 17's:

    ```{python}
    tf.stream_seed(1, ["subject", 18, "ae"])
    ```
    """
    return _core.stream_seed(seed, path)


def core_version() -> str:
    """The version of the Rust core the package was built with.

    A study depends on the core's version as well as on its request: see the reproducibility
    guide for what a version change can alter. The package's `__version__` is the same
    number.

    Returns
    -------
    str
        The version, e.g. `"0.1.0"`.

    Examples
    --------
    ```{python}
    import trialforge as tf

    tf.core_version()
    ```

    Record it with a study's request, alongside the archetype, `n`, seed and controls, to
    know which version to install to generate the study again.
    """
    return _core.core_version()
