from __future__ import annotations

import numpy as np
import pandas as pd

from scripts.disstree_demo import load_labels, load_matrix


def test_load_matrix_accepts_documented_named_header_or_no_header(tmp_path):
    expected = np.array([[0.0, 1.0], [1.0, 0.0]])
    with_header = tmp_path / "with-header.csv"
    without_header = tmp_path / "without-header.csv"
    pd.DataFrame(expected, columns=["case_a", "case_b"]).to_csv(with_header, index=False)
    pd.DataFrame(expected).to_csv(without_header, index=False, header=False)

    np.testing.assert_array_equal(load_matrix(str(with_header)), expected)
    np.testing.assert_array_equal(load_matrix(str(without_header)), expected)


def test_load_labels_preserves_headerless_first_value_and_named_label_header(tmp_path):
    expected = np.array(["A", "A", "B"])
    with_header = tmp_path / "with-header.csv"
    without_header = tmp_path / "without-header.csv"
    pd.DataFrame({"label": expected}).to_csv(with_header, index=False)
    pd.Series(expected).to_csv(without_header, index=False, header=False)

    np.testing.assert_array_equal(load_labels(str(with_header)), expected)
    np.testing.assert_array_equal(load_labels(str(without_header)), expected)
