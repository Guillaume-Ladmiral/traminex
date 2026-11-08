"""Index bases of the index-bearing public callables (TX-AUD-20261019-007, -066).

``pairwise(ref_sets=)`` takes zero-based row positions and ``seqdist(refseq=)`` R's one-based
indices; docs/DISTANCES.md once passed the same ``(a, b)`` to both.
"""

from __future__ import annotations

import inspect
import re
import warnings

import numpy as np
import pytest

import traminex as tx
from traminex.api_manifest import PUBLIC_API
from traminex.distances import pairwise


def test_pairwise_ref_sets_example_equals_seqdist() -> None:
    # The documented example: the same block from the zero-based and the one-based route.
    x = tx.seqdef(
        [list("AABBC"), list("ABBCC"), list("CCBAA"), list("BBCCA"), list("AAAAB")],
        alphabet=list("ABC"),
    )
    a, b = [1, 2], [3, 4]
    a0, b0 = [i - 1 for i in a], [i - 1 for i in b]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from_pairwise = pairwise(
            np.asarray(x.data), method="OM", subm="TRATE", indel=1.0, ref_sets=(a0, b0)
        )
        from_seqdist = tx.seqdist(x, method="OM", subm="TRATE", indel=1.0, refseq=[a, b])
    np.testing.assert_allclose(from_pairwise, np.asarray(from_seqdist), rtol=0, atol=1e-12)


# The public callables that take or return row positions, ranks or codes, classified by reading
# each hit of a name-and-docstring census (TX-AUD-20261019-066); "zero" and "one" name the base
# each docstring must state. Ids, labels, grid sizes and indicator names ("complexity index")
# are not positions and are not listed.
_ZERO, _ONE = "zero", "one"
INDEX_BEARING: dict[str, set[str]] = {
    "STSList": {_ONE},
    "seqdef": {_ZERO},
    "TSE_to_STS": {_ONE},
    "FCE_to_TSE": {_ONE},
    "STS_to_SPELL": {_ONE},
    "seqtab": {_ZERO, _ONE},
    "seqtabstocc": {_ONE},
    "sorti": {_ZERO},
    "sortv": {_ZERO},  # batch-8 review M-5: missing from the first census
    "seqasnum": {_ZERO},
    "seqstatl": {_ONE},
    "seqfind": {_ZERO},
    "seqfpos": {_ONE},
    "seqfposend": {_ONE},
    "seqconc": {_ONE},
    "seqpm": {_ONE},
    "seqsamm": {_ONE},
    "seqmaintokens": {_ZERO},
    "seqalign": {_ZERO},
    "seqalignplot": {_ZERO},
    "seqstart": {_ZERO, _ONE},
    "toPersonPeriod": {_ZERO, _ONE},
    "sequence_index_plot": {_ZERO, _ONE},
    "seqpcplot": {_ONE},
    "seqrep": {_ZERO},
    "seqrep_grp": {_ZERO},
    "medoid_index": {_ZERO},
    "cluster_agglomerative": {_ZERO},
    "cluster_representatives": {_ZERO},
    "cluster_summary": {_ZERO},
    "seqrf": {_ZERO},
    "dissrf": {_ZERO},
    "seqdist": {_ZERO, _ONE},
    "disscenter": {_ZERO},
    "weighted_pairwise_sum": {_ZERO},
    "weighted_inertia": {_ZERO},
    "dissrep": {_ZERO},
    "disstree": {_ZERO},
    "seqtree": {_ZERO},
    "disstree_assign": {_ONE},
}
_BASE_WORDS = {_ZERO: re.compile(r"zero-based|0-based"), _ONE: re.compile(r"one-based|1-based")}


@pytest.mark.parametrize("name", sorted(INDEX_BEARING))
def test_index_bearing_callables_state_their_base(name: str) -> None:
    assert name in PUBLIC_API
    doc = inspect.getdoc(getattr(tx, name)) or ""
    missing = [base for base in INDEX_BEARING[name] if not _BASE_WORDS[base].search(doc)]
    assert not missing, f"{name} does not state its {missing} base"
