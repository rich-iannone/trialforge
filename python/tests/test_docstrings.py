"""The examples in the public docstrings run and give the output they show."""

import doctest
from types import ModuleType

import pytest

from trialforge import _seeds, _study


@pytest.mark.parametrize("module", [_study, _seeds], ids=lambda m: m.__name__)
def test_docstring_examples(module: ModuleType) -> None:
    flags = doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE
    result = doctest.testmod(module, optionflags=flags, raise_on_error=False, verbose=False)
    assert result.attempted > 0
    assert result.failed == 0
