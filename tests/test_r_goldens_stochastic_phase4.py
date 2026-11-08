"""Focused evidence for public stochastic contracts re-audited in Phase 4."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import get_args, get_type_hints

import numpy as np

from traminex import dissrep, seqdef, seqgen
from traminex.cluster import seqrep

GOLDEN = Path(__file__).resolve().parent / "r_goldens" / "stochastic_phase4.csv"


def _golden_rows(case: str, field: str) -> list[str]:
    with GOLDEN.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    selected = [row for row in rows if row["case"] == case and row["field"] == field]
    return [row["value"] for row in sorted(selected, key=lambda row: int(row["position"]))]


def _distance_fixture() -> np.ndarray:
    points = np.array([0.0, 0.7, 1.9, 3.6, 5.8, 8.5])
    return np.abs(points[:, None] - points[None, :])


def test_seqgen_is_locally_reproducible_but_not_the_r_global_stream() -> None:
    kwargs = {"p": [0.55, 0.30, 0.15], "seed": 24680}
    first = seqgen(6, 7, ["A", "B", "C"], **kwargs)
    second = seqgen(6, 7, ["A", "B", "C"], **kwargs)

    np.testing.assert_array_equal(first, second)
    assert first.shape == (6, 7)
    assert first.dtype == object
    assert set(first.ravel()) <= {"A", "B", "C"}
    assert first.ravel().tolist() != _golden_rows("seqgen", "draw")
    assert _golden_rows("seqgen", "repeat_identical") == ["TRUE"]
    assert _golden_rows("seqgen", "global_rng_advanced") == ["TRUE"]


def test_dissrep_random_uses_repeatable_python_global_draws_with_r_invariants() -> None:
    distance = _distance_fixture()
    original = distance.copy()

    np.random.seed(97531)
    first = dissrep(distance, criterion="random", nrep=3, pradius=0)
    np.random.seed(97531)
    second = dissrep(distance, criterion="random", nrep=3, pradius=0)

    np.testing.assert_array_equal(first["representatives"], second["representatives"])
    np.testing.assert_array_equal(distance, original)
    assert first["representatives"].shape == (3,)
    assert len(set(first["representatives"].tolist())) == 3
    assert first["representatives"].tolist() != [
        int(value) for value in _golden_rows("dissrep", "representative_zero_based")
    ]
    assert _golden_rows("dissrep", "repeat_identical") == ["TRUE"]
    assert _golden_rows("dissrep", "global_rng_advanced") == ["TRUE"]


def test_dissrep_random_state_is_local_reproducible_and_accepts_generator() -> None:
    distance = _distance_fixture()

    np.random.seed(713)
    expected_global = np.random.random(4)
    np.random.seed(713)
    first = dissrep(distance, criterion="random", nrep=3, pradius=0, random_state=911)
    second = dissrep(distance, criterion="random", nrep=3, pradius=0, random_state=911)

    np.testing.assert_array_equal(first["representatives"], second["representatives"])
    np.testing.assert_array_equal(np.random.random(4), expected_global)

    generator_a = np.random.default_rng(912)
    generator_b = np.random.default_rng(912)
    generated_a = dissrep(distance, criterion="random", nrep=3, pradius=0, random_state=generator_a)
    generated_b = dissrep(distance, criterion="random", nrep=3, pradius=0, random_state=generator_b)
    np.testing.assert_array_equal(generated_a["representatives"], generated_b["representatives"])


def test_seqrep_random_is_typed_and_preserves_public_result_schema() -> None:
    criterion = get_type_hints(seqrep)["criterion"]
    assert {"freq", "random"} <= set(get_args(criterion))

    x = seqdef(
        [
            ["A", "A", "B"],
            ["A", "B", "B"],
            ["B", "B", "C"],
            ["B", "C", "C"],
            ["C", "C", "A"],
            ["C", "A", "A"],
        ],
        alphabet=["A", "B", "C"],
    )
    distance = _distance_fixture()
    original_data = x.data.copy()

    np.random.seed(86420)
    first = seqrep(x, diss=distance, criterion="random", nrep=3, pradius=0)
    np.random.seed(86420)
    second = seqrep(x, diss=distance, criterion="random", nrep=3, pradius=0)

    np.testing.assert_array_equal(first["rep_indices"], second["rep_indices"])
    np.testing.assert_array_equal(first["scores"], second["scores"])
    np.testing.assert_array_equal(first["representatives"], x.data[first["rep_indices"]])
    np.testing.assert_array_equal(x.data, original_data)
    assert first["scores"].shape == (x.nseq,)
    assert first["metadata"]["criterion"] == "random"
