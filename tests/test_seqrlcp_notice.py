"""seqRLCP prints seqdist's with.missing notice only when it applies (TX-AUD-20261019-025).

seqRLCP (a traminex helper; no R function) called seqdist with ``with_missing=True`` on every call, so data without
missing cells always printed R's " [!!] seqdist: 'with.missing' set as FALSE ..." notice; R's helpers seqLLCS and
seqLLCP print none. It now passes ``with_missing=True`` only when the data hold missing cells.
"""

from __future__ import annotations

import warnings

import numpy as np

import traminex as tx


def test_missing_free_data_print_no_notice() -> None:
    x = tx.seqdef(np.array([list("AABBC"), list("ABBCC"), list("CCBAA")], dtype=object))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = tx.seqRLCP(x)
    np.testing.assert_array_equal(result, tx.seqdist(x, method="RLCP"))


def test_data_with_missing_cells_still_compare_them() -> None:
    x = tx.seqdef(
        np.array(
            [["A", None, "B", "C"], ["A", "B", "B", "C"], ["C", "C", None, "A"]], dtype=object
        ),
        alphabet=["A", "B", "C"],
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = tx.seqRLCP(x)
    np.testing.assert_array_equal(result, tx.seqdist(x, method="RLCP", with_missing=True))
