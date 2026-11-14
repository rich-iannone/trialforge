"""The executable examples in the public docstrings run.

The reference renders each docstring's `{python}` cells, in order, in one session; here they run
the same way, so an example that would fail the docs build fails a test first.
"""

import inspect
import re
from collections.abc import Iterator

import pytest

import trialforge as tf

_CELL = re.compile(r"^```\{python\}\n(.*?)^```$", re.MULTILINE | re.DOTALL)


def _docstrings() -> Iterator[tuple[str, str]]:
    for name in tf.__all__:
        obj = getattr(tf, name)
        if inspect.isfunction(obj) or inspect.isclass(obj):
            yield name, inspect.getdoc(obj) or ""
        if inspect.isclass(obj):
            for member, value in vars(obj).items():
                public = not member.startswith("_")
                if public and (isinstance(value, property) or inspect.isfunction(value)):
                    yield f"{name}.{member}", inspect.getdoc(value) or ""


EXAMPLES = {name: cells for name, doc in _docstrings() if (cells := _CELL.findall(doc))}


def test_public_functions_have_examples() -> None:
    functions = {name for name in tf.__all__ if inspect.isfunction(getattr(tf, name))}
    assert functions <= set(EXAMPLES)
    assert {"Study.write", "Study.verify", "Study.score"} <= set(EXAMPLES)


@pytest.mark.parametrize("name", sorted(EXAMPLES))
def test_docstring_examples_run(name: str) -> None:
    namespace: dict[str, object] = {}
    for cell in EXAMPLES[name]:
        exec(compile(cell, f"<{name} example>", "exec"), namespace)
