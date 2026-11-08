import pandas as pd
import pytest

import traminex as tx


def _decode_states(x):
    return [[str(x.alphabet[int(v)]) for v in row] for row in x.data]


def test_hspell_vertical_basic():
    df = pd.DataFrame(
        {
            "id": [1, 1, 1, 2],
            "start": [1, 3, 5, 2],
            "end": [2, 4, 6, 3],
            "state": ["A", "B", "A", "B"],
        }
    )
    X = tx.HSPELL_to_STS(df, id_col="id", start_col="start", end_col="end", state_col="state")
    assert X.nseq == 2
    assert X.ids == ["1", "2"]
    # Check grid covers 1..6
    assert X.seqlen == 6


def test_hspell_horizontal_overwrite_and_fill():
    # Two spells per person in horizontal format
    df = pd.DataFrame(
        {
            "id": ["p1", "p2"],
            "begin1": [1, 1],
            "end1": [3, 2],
            "status1": ["A", "A"],
            "begin2": [2, 3],
            "end2": [4, 4],
            "status2": ["B", "B"],
        }
    )
    X = tx.HSPELL_to_STS(
        df,
        id_col="id",
        begin=["begin1", "begin2"],
        end=["end1", "end2"],
        status=["status1", "status2"],
        overwrite=True,
        fillblanks="_",
    )
    assert X.nseq == 2
    assert X.ids == ["p1", "p2"]
    assert _decode_states(X) == [["A", "B", "B", "B"], ["A", "A", "B", "B"]]


def test_hspell_endobs_numeric_id_horizontal():
    df = pd.DataFrame(
        {
            "id": [1, 2],
            "begin1": [1, 1],
            "end1": [2, 2],
            "status1": ["A", "A"],
            "begin2": [3, 3],
            "end2": [4, 4],
            "status2": ["B", "B"],
            "endObs": [2, 2],
        }
    )
    X = tx.HSPELL_to_STS(
        df,
        id_col="id",
        begin=["begin1", "begin2"],
        end=["end1", "end2"],
        status=["status1", "status2"],
        fillblanks="_",
        endObs="endObs",
    )
    # Live TraMineRextras 0.6.9 blanks every cell from time endObs - 1 on
    # (seqindex >= endObs - tmin), after fillblanks, so endObs = 2 leaves nothing.
    # traminex blanks from endObs on (TX-REP-20261017-005): time 1 is observed,
    # times 2..4 become NA, i.e. void once seqdef drops the trailing NAs.
    assert _decode_states(X) == [["A", "%", "%", "%"], ["A", "%", "%", "%"]]


def test_hspell_pvar_realign_preserves_rendered_row_ids():
    df = pd.DataFrame(
        {
            "id": ["case-z", "case-a"],
            "begin1": [1, 1],
            "end1": [2, 2],
            "status1": ["A", "B"],
            "origin": [1, 1],
        }
    )

    X = tx.HSPELL_to_STS(
        df,
        id_col="id",
        begin=["begin1"],
        end=["end1"],
        status=["status1"],
        pvar="origin",
        missing_color="orange",
    )

    # One row per input row in input order, as R keeps it (this used to pin the
    # ids sorted as text; 2026-09-22 audit S2-01).
    assert X.ids == ["case-z", "case-a"]
    assert X.missing_color == "orange"
    assert dict(zip(X.ids, _decode_states(X), strict=True)) == {
        "case-a": ["B", "B"],
        "case-z": ["A", "A"],
    }


def test_hspell_endobs_numeric_id_vertical():
    df = pd.DataFrame(
        {
            "id": [1, 1, 2, 2],
            "start": [1, 3, 1, 3],
            "end": [2, 4, 2, 4],
            "state": ["A", "B", "A", "B"],
            "endObs": [2, 2, 2, 2],
        }
    )
    X = tx.HSPELL_to_STS(
        df,
        id_col="id",
        start_col="start",
        end_col="end",
        state_col="state",
        fillblanks="_",
        endObs="endObs",
    )
    # Vertical spells give the horizontal result (TX-AUD-20261011-019): the
    # periods from endObs = 2 on are blanked after the filler (TX-REP-20261017-005),
    # so only time 1 stays and the trailing NAs become void.
    assert _decode_states(X) == [["A", "%", "%", "%"], ["A", "%", "%", "%"]]


@pytest.mark.parametrize("id_type", [int, float])
@pytest.mark.parametrize("mode", ["bounds", "endobs", "origin"])
def test_vertical_numeric_ids_use_same_identity_for_all_person_columns(id_type, mode):
    df = pd.DataFrame(
        {
            "id": [id_type(1), id_type(2)],
            "start": [1, 1],
            "end": [4, 4],
            "state": ["A", "B"],
            "lower": [2, 1],
            "upper": [3, 4],
            "observed": [2, 4],
            "origin": [1, 2],
        }
    )
    original = df.copy(deep=True)
    options = {
        "bounds": {"tmin_col": "lower", "tmax_col": "upper"},
        "endobs": {"endObs": "observed"},
        "origin": {"pvar": "origin"},
    }[mode]
    result = tx.HSPELL_to_STS(
        df,
        id_col="id",
        start_col="start",
        end_col="end",
        state_col="state",
        fillblanks="_",
        **options,
    )
    # The horizontal route's rules (TX-AUD-20261011-019): per-person bounds blank
    # outside [lower, upper] after the filler, as endObs does (TX-AUD-20261019-037),
    # leaving missing before and void after; endObs blanks from endObs on after
    # it (TX-REP-20261017-005), leaving void; pvar subtracts the origin from the
    # times, so the grid is -1..3 (person 1 covers 0..3, person 2 covers -1..2).
    expected = {
        "bounds": [["*", "A", "A", "%"], ["B", "B", "B", "B"]],
        "endobs": [["A", "%", "%", "%"], ["B", "B", "B", "%"]],
        "origin": [["_", "A", "A", "A", "A"], ["B", "B", "B", "B", "_"]],
    }[mode]
    assert result.ids == ["1", "2"]
    assert _decode_states(result) == expected
    pd.testing.assert_frame_equal(df, original)
