import pytest

import trialforge


def test_versions_agree() -> None:
    assert trialforge.__version__ == trialforge.core_version()


def test_stream_seed_is_deterministic() -> None:
    path = ["subject", 17, "ae"]
    assert trialforge.stream_seed(1, path) == trialforge.stream_seed(1, path)
    assert trialforge.stream_seed(1, path) != trialforge.stream_seed(2, path)


def test_stream_seed_rejects_bad_parts() -> None:
    with pytest.raises(TypeError):
        trialforge.stream_seed(1, [1.5])  # type: ignore[list-item]
