"""Public metadata-bearing result contract for ``seqrf``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from r_goldens_utils import load_golden

from traminex import STSList, seqrf

ROOT = Path(__file__).resolve().parents[1]


def _metadata_fixture() -> STSList:
    codes = load_golden("seqrf_input_codes").astype(int)
    return STSList(
        data=codes,
        alphabet=["A", "B", "C", "D", "*", "%"],
        labels=["Alpha", "Beta", "Gamma", "Delta", "Missing", "Void"],
        cpal=["#112233", "#445566", "#778899", "#AABBCC", "#DDEEFF", "#000000"],
        xtlab=[f"Wave {index}" for index in range(1, codes.shape[1] + 1)],
        xtstep=2,
        tick_last=True,
        nr="*",
        void="%",
        nr_code=4,
        void_code=5,
        weights=load_golden("seqrf_input_weights").ravel(),
        start=7,
        ids=[f"case-{index:02d}" for index in range(1, codes.shape[0] + 1)],
        missing_color="#CC00CC",
        issues={"flag": [1, 11, 13, 15], "other": [2, 14]},
    )


def _result(source: STSList) -> dict:
    return seqrf(
        source,
        load_golden("seqrf_input_diss"),
        k=3,
        sortv=load_golden("seqrf_input_sortv").ravel(),
        weighted=True,
    )


def test_seqrf_adds_one_metadata_complete_medoid_sequence_field() -> None:
    source = _metadata_fixture()
    result = _result(source)

    assert set(result) == {
        "seqtoplot",
        "rf",
        "medoids",
        "heights",
        "medoid_sequences",
    }
    medoids = result["medoids"]
    representatives = result["medoid_sequences"]
    assert isinstance(representatives, STSList)
    np.testing.assert_array_equal(representatives.data, source.data[medoids])
    np.testing.assert_array_equal(result["seqtoplot"], source.data[medoids])
    np.testing.assert_allclose(representatives.weights, result["rf"]["heights"])
    assert representatives.ids == [source.ids[index] for index in medoids]
    assert representatives.alphabet == source.alphabet
    assert representatives.labels == source.labels
    assert representatives.cpal == source.cpal
    assert representatives.xtlab == source.xtlab
    assert (
        representatives.xtstep,
        representatives.tick_last,
        representatives.start,
    ) == (2, True, 7)
    assert (
        representatives.nr,
        representatives.void,
        representatives.nr_code,
        representatives.void_code,
    ) == ("*", "%", 4, 5)
    assert representatives.missing_color == "#CC00CC"
    assert representatives.issues == {"flag": [1, 3], "other": [2]}


def test_seqrf_medoid_sequences_are_mutation_isolated() -> None:
    source = _metadata_fixture()
    result = _result(source)
    medoids = result["medoids"].copy()
    legacy_seqtoplot = result["seqtoplot"].copy()
    representatives = result["medoid_sequences"]

    representatives.data[0, 0] = 99
    representatives.alphabet[0] = "changed"
    representatives.labels[0] = "changed"
    representatives.cpal[0] = "changed"
    representatives.xtlab[0] = "changed"
    representatives.weights[0] = 99
    representatives.ids[0] = "changed"
    representatives.issues["flag"].append(99)

    np.testing.assert_array_equal(source.data[medoids], legacy_seqtoplot)
    np.testing.assert_array_equal(result["seqtoplot"], legacy_seqtoplot)
    assert source.alphabet[0] == "A"
    assert source.labels[0] == "Alpha"
    assert source.cpal[0] == "#112233"
    assert source.xtlab[0] == "Wave 1"
    assert source.ids[medoids[0]] != "changed"
    assert source.issues == {"flag": [1, 11, 13, 15], "other": [2, 14]}
    assert result["rf"]["heights"][0] != 99


def test_visualization_contract_documents_public_seqrf_field_and_ownership() -> None:
    contract = (ROOT / "docs/VISUALIZATION_SERIES_CONTRACT.md").read_text(encoding="utf-8")

    assert "`medoid_sequences`" in contract
    assert "metadata-complete" in contract
    assert "renderer-owned drawing summaries" in contract
