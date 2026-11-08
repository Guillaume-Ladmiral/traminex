"""Resource-boundary regressions for event-subsequence matching."""

import json
import subprocess
import sys
import textwrap

import numpy as np
import pandas as pd
import pytest

import traminex as tx


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!\]\ count\ method\ \d+\ is\ only\ for\ internal\ use\.:UserWarning"
)
@pytest.mark.parametrize(
    "count_method",
    ["COBJ", "CDIST_O", "CWIN", "CMINWIN", "CDIST", 6],
)
def test_long_exact_event_subsequence_does_not_depend_on_python_recursion_limit(count_method):
    event_count = 1_100
    labels = [f"E{index:04d}" for index in range(event_count)]
    event_sequence = tx.seqecreate(
        pd.DataFrame(
            {
                "id": ["one"] * event_count,
                "timestamp": np.arange(1, event_count + 1, dtype=float),
                "event": labels,
            }
        )
    )
    pattern = "(" + ")-(".join(labels) + ")"
    constraint = tx.seqeconstraint(
        count_method=count_method,
        window_size=2_000 if count_method == "CWIN" else -1,
    )

    result = tx.seqefsub(event_sequence, str_subseq=[pattern], constraint=constraint)

    assert result.data.shape[0] == 1
    assert np.isfinite(result.data.loc[0, "Support"])
    assert np.isfinite(result.data.loc[0, "Count"])


def _repeated_event_sequence(event_count):
    return tx.seqecreate(
        pd.DataFrame(
            {
                "id": ["one"] * event_count,
                "timestamp": np.arange(event_count, dtype=float),
                "event": ["A"] * event_count,
            }
        )
    )


def test_seqefsub_prefix_tree_handles_1100_event_chain_at_default_recursion_limit():
    event_count = 1_100

    result = tx.seqefsub(
        _repeated_event_sequence(event_count),
        min_support=1,
        constraint=tx.seqeconstraint(max_gap=1, age_min=0, age_max=0),
        max_k=event_count,
    )

    assert len(result) == event_count
    assert result.subseq[0] == "(A)"
    assert result.subseq[-1].count("(A)") == event_count
    np.testing.assert_array_equal(result.data["Support"], np.ones(event_count))
    np.testing.assert_array_equal(result.data["Count"], np.ones(event_count))


@pytest.mark.parametrize("recursion_limit", [80, 120])
def test_seqefsub_prefix_tree_result_is_independent_of_lower_recursion_limit(recursion_limit):
    program = textwrap.dedent(f"""
        import json
        import sys
        import numpy as np
        import pandas as pd
        import traminex as tx

        event_count = 150
        events = tx.seqecreate(pd.DataFrame({{
            "id": ["one"] * event_count,
            "timestamp": np.arange(event_count, dtype=float),
            "event": ["A"] * event_count,
        }}))
        sys.setrecursionlimit({recursion_limit})
        result = tx.seqefsub(
            events,
            min_support=1,
            constraint=tx.seqeconstraint(max_gap=1, age_min=0, age_max=0),
            max_k=event_count,
        )
        print(json.dumps({{
            "length": len(result),
            "first": result.subseq[0],
            "last_size": result.subseq[-1].count("(A)"),
            "support_sum": float(result.data["Support"].sum()),
            "count_sum": float(result.data["Count"].sum()),
        }}, sort_keys=True))
        """)

    completed = subprocess.run(
        [sys.executable, "-c", program],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout) == {
        "count_sum": 150.0,
        "first": "(A)",
        "last_size": 150,
        "length": 150,
        "support_sum": 150.0,
    }
