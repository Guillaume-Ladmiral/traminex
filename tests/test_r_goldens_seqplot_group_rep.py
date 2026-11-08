"""seqplot(group=) prepares types "r" and "rf" on the whole set, as R does (TX-AUD-20261010-044).

Before it splits the sequences by group, R's seqplot computes the distances
from ``method`` and ``dmax = max(diss)`` for type "r", and the default MDS sort
variable for type "rf" (seqplot.R:96-156); each panel then receives its subset
(seqplot.R:331-339, :368). It stops for type "r" when a group has fewer than
two cases (seqplot.R:97-101). traminex computed those values on each panel's
own subset, so the grouped plots showed other medoids and representatives
than R's. The panel inputs are read from the renderer's calls to seqrf and
seqrep. A one-dimensional MDS axis has no sign of its own
(TX-REP-20261017-006), so the rf medoids are compared with R's grouped plot
for the sign of traminex's axis. Fixture:
``scripts/r_goldens/make_goldens_seqplot_group_rep.R``.
"""

from __future__ import annotations

import importlib
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

import traminex as tx  # noqa: E402

representatives = importlib.import_module("traminex.visualization.plots.representatives")

HERE = Path(__file__).resolve().parent / "r_goldens"
GOLDEN = pd.read_csv(HERE / "seqplot_group_rep.csv", dtype=str, keep_default_na=False)
INPUT = pd.read_csv(HERE / "seqplot_group_rep_input.csv", dtype=str, keep_default_na=False)
DISS = np.loadtxt(HERE / "seqplot_group_rep_diss.csv", delimiter=",")
ALPHABET = [str(state) for state in range(8)]
ROWS = np.array([sequence.split("-") for sequence in INPUT["sequence"]], dtype=object)
SEQS = tx.seqdef(ROWS, alphabet=ALPHABET)
WEIGHTED = tx.seqdef(ROWS, alphabet=ALPHABET, weights=INPUT["weight"].astype(float).to_numpy())
GROUP = INPUT["sex"].tolist()
LEVELS = ["man", "woman"]
R_OPTIONS = {
    "r_lcs_dmax": {"diss": DISS, "coverage": 0.5, "pradius": 0.2},
    # R's sm is traminex's subm.
    "r_trate_method": {"method": "OM", "subm": "TRATE", "indel": 1},
}


def _values(case: str, level: str, quantity: str) -> list[str]:
    rows = GOLDEN[
        (GOLDEN["case"] == case) & (GOLDEN["level"] == level) & (GOLDEN["quantity"] == quantity)
    ]
    return rows.sort_values("position", key=lambda c: c.astype(int))["value"].tolist()


def _numbers(case: str, level: str, quantity: str) -> np.ndarray:
    values = [value.strip() for value in _values(case, level, quantity)]
    return np.array([np.nan if value == "NA" else float(value) for value in values])


def _index(level: str) -> list[int]:
    return [i for i, value in enumerate(GROUP) if value == level]


def _spy(monkeypatch, name: str) -> list[tuple[dict, dict]]:
    original = getattr(representatives, name)
    calls: list[tuple[dict, dict]] = []

    def spy(*args, **kwargs):
        result = original(*args, **kwargs)
        calls.append((kwargs, result))
        return result

    monkeypatch.setattr(representatives, name, spy)
    return calls


@pytest.mark.parametrize("case, seqs", [("rf_k4", SEQS), ("rf_weighted_k4", WEIGHTED)])
def test_rf_panels_sort_on_the_whole_set_axis(case, seqs, monkeypatch) -> None:
    calls = _spy(monkeypatch, "seqrf")
    fig, _ = tx.seqplot(seqs, type="rf", group=GROUP, diss=DISS, k=4)
    plt.close(fig)
    # R's wcmdscale gives a zero-weight case a missing coordinate (NA, read as
    # NaN); R's order() puts it last, as NumPy's sort puts NaN. Copies of one
    # sequence share a coordinate up to rounding (about 5e-14 here), so the
    # orders and medoids are compared as sequences, not row numbers.
    axis = _numbers(case, "all", "global_axis")
    assert len(calls) == len(LEVELS)
    signs = set()
    for (kwargs, result), level in zip(calls, LEVELS):
        assert not isinstance(kwargs["sortv"], str), "the panel computed its own MDS"
        texts = INPUT["sequence"].to_numpy()[_index(level)]
        order = texts[np.argsort(np.asarray(kwargs["sortv"], dtype=float), kind="stable")]
        own = texts[np.argsort(axis[_index(level)], kind="stable")].tolist()
        negated = texts[np.argsort(-axis[_index(level)], kind="stable")].tolist()
        assert order.tolist() in (own, negated)
        suffix = "" if order.tolist() == own else "_negated"
        signs.add(suffix)
        medoids = [int(m) - 1 for m in _values(case, level, "medoids" + suffix)]
        assert texts[np.asarray(result["medoids"])].tolist() == texts[medoids].tolist()
        np.testing.assert_allclose(
            np.asarray(result["heights"], dtype=float),
            _numbers(case, level, "heights" + suffix),
            rtol=1e-12,
        )
    assert len(signs) == 1


@pytest.mark.parametrize("case", sorted(R_OPTIONS))
def test_r_panels_take_whole_set_distances_and_dmax(case, monkeypatch) -> None:
    calls = _spy(monkeypatch, "seqrep")
    fig, _ = tx.seqplot(SEQS, type="r", group=GROUP, **R_OPTIONS[case])
    plt.close(fig)
    assert len(calls) == len(LEVELS)
    for (kwargs, result), level in zip(calls, LEVELS):
        assert kwargs.get("diss") is not None, "the panel computed its own distances"
        assert float(np.sum(kwargs["diss"])) == pytest.approx(
            _numbers(case, level, "diss_sum")[0], rel=1e-12
        )
        assert kwargs.get("dmax") == pytest.approx(_numbers(case, level, "dmax")[0], rel=1e-12)
        drawn = ["-".join(ALPHABET[code] for code in row) for row in result["representatives"]]
        assert drawn == _values(case, level, "representatives")
        np.testing.assert_allclose(
            result["statistics"]["na(%)"].to_numpy(dtype=float),
            _numbers(case, level, "na_pct"),
            rtol=1e-12,
        )


def test_the_fixture_tells_the_whole_set_rule_from_the_per_panel_rule() -> None:
    rf_apart = any(
        _values("rf_k4", level, "subset_medoids")
        not in (_values("rf_k4", level, "medoids"), _values("rf_k4", level, "medoids_negated"))
        for level in LEVELS
    )
    r_apart = _values("r_lcs_dmax", "man", "representatives") != _values(
        "r_lcs_dmax", "man", "subset_representatives"
    )
    method_apart = _values("r_trate_method", "man", "diss_sum") != _values(
        "r_trate_method", "man", "subset_diss_sum"
    )
    assert rf_apart and r_apart and method_apart


@pytest.mark.parametrize("function", [tx.seqplot, tx.seqrplot])
def test_r_refuses_a_group_with_one_case(function) -> None:
    message = _values("r_one_case", "all", "error")[0]
    options = {"type": "r"} if function is tx.seqplot else {}
    with pytest.raises(ValueError, match=re.escape(message)):
        function(SEQS, group=["g"] + ["h"] * 59, diss=DISS, **options)


@pytest.mark.parametrize("function", [tx.seqplot, tx.seqrplot])
def test_r_refuses_one_sequence_without_group(function) -> None:
    # Without group, R's check counts the whole set as one group.
    message = _values("r_one_sequence", "all", "error")[0]
    options = {"type": "r"} if function is tx.seqplot else {}
    one = tx.seqdef(ROWS[:1], alphabet=ALPHABET)
    with pytest.raises(ValueError, match=re.escape(message)):
        function(one, diss=DISS[:1, :1], **options)
