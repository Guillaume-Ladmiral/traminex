"""``cluster_summary`` and ``cluster_representatives`` keep the cluster labels.

Audit 2026-09-22, finding S4-14. Both traminex-only helpers passed every
cluster label through ``int()``: text labels raised ``ValueError`` and
fractional labels were truncated, so the distinct clusters 1.5 and 1.7 were
both reported as cluster ``1``. The helpers document one row per cluster
labelled by that cluster, so the label must be kept as given (integer labels
still come back as Python ``int``, as before).
"""

from __future__ import annotations

import numpy as np
import pytest

from traminex import seqdef
from traminex.cluster import cluster_representatives, cluster_summary

SEQUENCES = seqdef(np.array([list("AAB"), list("ABB"), list("BBB"), list("BAA")], dtype=object))
D = np.array([[0, 1, 2, 2], [1, 0, 1, 2], [2, 1, 0, 2], [2, 2, 2, 0]], dtype=float)
LABELS = {
    "int": (np.array([1, 1, 2, 2]), [1, 2], int),
    "str": (np.array(["high", "high", "low", "low"]), ["high", "low"], str),
    "float": (np.array([1.5, 1.5, 1.7, 1.7]), [1.5, 1.7], float),
}


@pytest.mark.parametrize("helper", [cluster_summary, cluster_representatives])
@pytest.mark.parametrize("kind", sorted(LABELS))
def test_cluster_labels_are_kept(helper, kind):
    labels, expected, python_type = LABELS[kind]
    result = helper(SEQUENCES, labels, D)
    assert result["cluster"].tolist() == expected
    assert all(type(value) is python_type for value in result["cluster"])
