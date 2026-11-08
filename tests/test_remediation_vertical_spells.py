import pandas as pd
import pytest

from traminex import HSPELL_to_STS


@pytest.mark.parametrize("overwrite", [True, False])
@pytest.mark.parametrize("priority", ["later", "earlier", "longer", "shorter"])
def test_vertical_overlap_respects_requested_ownership(overwrite, priority):
    frame = pd.DataFrame({"id": [1, 1], "begin": [1, 2], "end": [4, 3], "state": ["A", "B"]})
    result = HSPELL_to_STS(
        frame,
        id_col="id",
        start_col="begin",
        end_col="end",
        state_col="state",
        overwrite=overwrite,
        overlap_priority=priority,
    )
    labels = [result._state_for_code(int(code)) for code in result.data[0]]
    assert labels == list("ABBA" if overwrite and priority in {"later", "shorter"} else "AAAA")


@pytest.mark.parametrize(
    "priority,expected",
    [("later", "BBBA"), ("earlier", "BAAA"), ("longer", "BBBA"), ("shorter", "BAAA")],
)
def test_vertical_equal_duration_ties_follow_declared_rows(priority, expected):
    frame = pd.DataFrame({"id": [1, 1], "begin": [2, 1], "end": [4, 3], "state": ["A", "B"]})
    result = HSPELL_to_STS(
        frame,
        id_col="id",
        start_col="begin",
        end_col="end",
        state_col="state",
        overlap_priority=priority,
    )
    assert [result._state_for_code(int(code)) for code in result.data[0]] == list(expected)
