"""Public container diagnostics for state- and event-sequence consumers."""

from __future__ import annotations

import numpy as np
import pytest

import traminex as tx

STATE_SEQUENCE_CONSUMERS = (
    "STS_to_SPELL",
    "STS_to_SRS",
    "STS_to_TSE",
    "chronogram",
    "cluster_chronograms",
    "cluster_chronograms_with_representatives",
    "cpal",
    "duration_boxplot",
    "entropy",
    "frequent_ngrams",
    "overlay_representatives_on_index_plot",
    "seqHtplot",
    "seqST",
    "seqalignplot",
    "seqasnum",
    "seqauto",
    "seqcta",
    "seqdHplot",
    "seqdss",
    "seqdur",
    "seqetm",
    "seqfplot",
    "seqgranularity",
    "seqhasmiss",
    "seqici",
    "seqient",
    "seqindic",
    "seqindic_dyn",
    "seqinsecurity",
    "seqintegr",
    "seqistatd",
    "seqivolatility",
    "seqmaintokens",
    "seqmeant",
    "seqmodst",
    "seqmsplot",
    "seqmtplot",
    "seqnum",
    "seqplot",
    "seqplot_rf",
    "seqplot_tentrop",
    "seqprecarity",
    "seqrep",
    "seqrfplot",
    "seqrplot",
    "seqstatd",
    "seqstatf",
    "seqsubsn",
    "seqsurvplot",
    "seqtransn",
    "seqtransn_count",
    "sequence_index_plot",
    "sequence_index_plot_grouped",
    "state_distribution_plot",
    "state_legend",
    "state_legend_figure",
    "stlab",
    "toPersonPeriod",
    "transition_heatmap",
    "turbulence_exact",
    "turbulence_proxy",
)


@pytest.mark.parametrize("name", STATE_SEQUENCE_CONSUMERS)
def test_state_sequence_consumers_raise_explicit_container_type_errors(name):
    function = getattr(tx, name)

    with pytest.raises(TypeError, match="must be an STSList state sequence object"):
        function(np.asarray([[0, 1], [1, 0]], dtype=int))


def test_state_sequence_guard_preserves_keyword_argument_calls():
    sequences = tx.seqdef([["A", "B"], ["B", "A"]])

    result = tx.seqstatd(x=sequences)

    assert result["Frequencies"].shape == (2, 2)


def test_state_sequence_guard_preserves_missing_argument_diagnostic():
    with pytest.raises(TypeError, match="missing 1 required positional argument: 'x'"):
        tx.seqstatd()


def test_seqetm_deprecated_alias_uses_the_explicit_container_diagnostic():
    with pytest.warns(DeprecationWarning, match="seq"):
        with pytest.raises(TypeError, match="x must be an STSList state sequence object"):
            tx.seqetm(seq=np.asarray([[0, 1], [1, 0]], dtype=int))


@pytest.mark.parametrize(
    ("name", "expected_type"),
    [
        ("plot_eseq", "Eseq"),
        ("plot_seqelist", "SeqElist"),
        # seqpcplot also takes a state sequence object, as R does (TX-AUD-20261011-027).
        ("seqpcplot", r"Eseq, SeqElist\) or state-sequence \(STSList"),
        ("seqedplot", "SeqElist"),
        ("plot_subseqelist", "SubSeqElist"),
        ("plot_subseqelistchisq", "SubSeqElistChisq"),
    ],
)
def test_event_sequence_consumers_raise_explicit_container_type_errors(name, expected_type):
    function = getattr(tx, name)

    with pytest.raises(TypeError, match=expected_type):
        function(np.asarray([[0, 1], [1, 0]], dtype=int))


def test_seqplotmd_refuses_channels_that_are_not_sequence_objects():
    # 2026-09-22 audit S6-15: a list of plain arrays raised AttributeError.
    matrix = np.asarray([[0, 1], [1, 0]], dtype=int)

    with pytest.raises(TypeError, match="every channel must be an STSList"):
        tx.seqplotMD([matrix, matrix])
