import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDENS_DIR

from traminex import seqdef
from traminex.sequences import STSList


def test_seqdef_unknown_state_reports_actionable_error():
    mat = np.array([["A", "B"], ["B", "A"]], dtype=object)
    with pytest.raises(ValueError, match="not present in provided alphabet"):
        seqdef(mat, alphabet=["A"])


def test_seqdef_weights_length_must_match_nseq():
    mat = np.array([["A"], ["B"]], dtype=object)
    with pytest.raises(ValueError, match="weights must have length 2"):
        seqdef(mat, weights=[1.0])


def test_seqdef_weights_must_be_finite():
    mat = np.array([["A"], ["B"]], dtype=object)
    with pytest.raises(ValueError, match="finite numeric values"):
        seqdef(mat, weights=[1.0, np.nan])


def test_seqdef_weights_must_be_non_negative():
    mat = np.array([["A"], ["B"]], dtype=object)
    with pytest.raises(ValueError, match="non-negative"):
        seqdef(mat, weights=[1.0, -0.1])


def test_seqdef_rejects_unknown_or_lowercase_input_formats():
    mat = np.array([["A"]], dtype=object)
    for informat in ("wide", "sts", "BAD"):
        with pytest.raises(ValueError, match="informat"):
            seqdef(mat, informat=informat)


def test_seqdef_preserves_dataframe_ids_and_builds_auto_ids():
    frame = pd.DataFrame([["A"], ["B"]], index=["case-b", "case-a"])
    preserved = seqdef(frame)
    assert preserved.ids == ["case-b", "case-a"]

    automatic = seqdef(frame.to_numpy(), id="auto")
    assert automatic.ids == ["[1]", "[2]"]

    selected = seqdef(pd.DataFrame({"first": ["A"], "second": ["B"]}), var=[0])
    assert selected.alphabet == ["A"]
    selected_last = seqdef(
        pd.DataFrame({"first": ["A"], "second": ["B"]}),
        var=[-1],
    )
    assert selected_last.alphabet == ["B"]


def test_seqdef_inference_is_r_order_and_never_stringifies_null_cells():
    mat = np.array([["B", None, None], ["A", "B", None]], dtype=object)
    x = seqdef(mat, right=None)

    # R's order A, B (default_alphabet.csv case right_na; owner decision D2), then nr.
    assert x.alphabet == ["A", "B", "*"]
    assert "None" not in x.alphabet
    assert "nan" not in x.alphabet
    assert x.nr_code == 2
    assert x.void_code is None
    np.testing.assert_array_equal(x.data, np.array([[1, 2, 2], [0, 1, 2]]))


def test_seqdef_retains_preencoded_special_markers_outside_explicit_alphabet():
    x = seqdef(
        np.array([["A", "*", "%"]], dtype=object),
        alphabet=["A"],
        nr="*",
        void="%",
    )

    assert x.alphabet == ["A", "*", "%"]
    assert x.nr_code == 1
    assert x.void_code == 2
    np.testing.assert_array_equal(x.data, np.array([[0, 1, 2]]))


def test_seqdef_extends_real_labels_and_colors_for_internal_markers():
    x = seqdef(
        np.array([[None, "A", None]], dtype=object),
        alphabet=["A"],
        labels=["Alpha"],
        cpal=["#010203"],
        missing_color="#A0A0A0",
    )

    assert x.alphabet == ["A", "*", "%"]
    assert x.labels == ["Alpha", "*", "%"]
    assert x.cpal == ["#010203", "#A0A0A0", "white"]


def test_seqdef_accepts_explicit_labels_and_colors_for_inferred_markers():
    x = seqdef(
        np.array([["A", None, "%"]], dtype=object),
        alphabet=["A"],
        labels=["Alpha", "Nonresponse", "Structural void"],
        cpal=["#010203", "#A0A0A0", "#FFFFFF"],
    )

    assert x.alphabet == ["A", "*", "%"]
    assert x.labels == ["Alpha", "Nonresponse", "Structural void"]
    assert x.cpal == ["#010203", "#A0A0A0", "#FFFFFF"]


def test_seqdef_states_mapping_is_total_and_injective():
    mat = np.array([["A", "B"]], dtype=object)
    with pytest.raises(ValueError, match="same length"):
        seqdef(mat, alphabet=["A", "B"], states=["X"])
    with pytest.raises(ValueError, match="unique"):
        seqdef(mat, alphabet=["A", "B"], states=["X", "X"])


def test_seqdef_only_autodecompresses_recognized_one_column_separators():
    extended = seqdef(np.array([["AB"], ["CD"]], dtype=object))
    assert extended.data.shape == (2, 1)
    assert extended.alphabet == ["AB", "CD"]

    compressed = seqdef(np.array([["A-B"], ["C-D"]], dtype=object))
    assert compressed.data.shape == (2, 2)
    assert compressed.alphabet == ["A", "B", "C", "D"]

    sps = pd.DataFrame(
        [["(ignored,1)", "(A,1)", "(B,1)"]],
        columns=["ignored", "first", "second"],
    )
    negative = seqdef(sps, informat="SPS", var=[-2, -1])
    positive = seqdef(sps, informat="SPS", var=[1, 2])
    np.testing.assert_array_equal(negative.data, positive.data)
    assert negative.alphabet == positive.alphabet == ["A", "B"]


def test_seqdef_safely_pads_ragged_rows_before_missing_preprocessing():
    x = seqdef([["A", "@"], ["B"]], missing="@")

    assert x.alphabet == ["A", "B", "%"]
    assert x.void_code == 2
    np.testing.assert_array_equal(x.data, np.array([[0, 2], [1, 2]]))


def test_seqdef_only_refuses_neutral_when_that_region_is_present():
    unused = seqdef(np.array([["A", None]], dtype=object), left="NEUTRAL")
    assert unused.alphabet == ["A", "%"]
    np.testing.assert_array_equal(unused.data, np.array([[0, 1]]))

    with pytest.raises(ValueError, match="NEUTRAL"):
        seqdef(np.array([[None, "A"]], dtype=object), left="NEUTRAL")


def test_seqdef_stslist_input_uses_codes_to_disambiguate_marker_text():
    loose = STSList(
        data=np.array([[0, 1, 2]]),
        alphabet=["A", "*", "%"],
        nr="*",
        void="%",
        nr_code=None,
        void_code=None,
        issues={"source": [1]},
    )
    real_markers = seqdef(loose)
    # Real states in R's sort order (seqformat_input_cases.csv, sts_marker_states).
    assert real_markers.alphabet == ["*", "%", "A"]
    assert real_markers.nr_code is None
    assert real_markers.void_code is None
    assert real_markers.issues == {"source": [1]}

    encoded = loose.copy()
    encoded.nr_code = 1
    encoded.void_code = 2
    special_markers = seqdef(encoded)
    assert special_markers.alphabet == ["A", "*", "%"]
    assert special_markers.nr_code == 1
    assert special_markers.void_code == 2


def test_seqdef_explicit_mapping_keeps_code_aware_marker_text_real():
    loose = STSList(
        data=np.array([[0, 1]]),
        alphabet=["A", "*"],
        nr="*",
        void="%",
        nr_code=None,
        void_code=None,
    )

    mapped = seqdef(
        loose,
        alphabet=["A", "*"],
        states=["Alpha", "Literal star"],
    )

    assert mapped.alphabet == ["Alpha", "Literal star"]
    assert mapped.nr_code is None
    np.testing.assert_array_equal(mapped.data, np.array([[0, 1]]))


def test_seqdef_renames_code_aware_markers_that_collide_with_real_states():
    mixed = STSList(
        data=np.array([[0, 1, 2]]),
        alphabet=["*", "A", "encoded missing"],
        nr="*",
        void="%",
        nr_code=2,
        void_code=None,
    )

    x = seqdef(mixed, alphabet=["*", "A"])

    assert x.alphabet == ["*", "A", "<NR>"]
    assert x.nr == "<NR>"
    assert x.nr_code == 2
    np.testing.assert_array_equal(x.data, np.array([[0, 1, 2]]))

    unobserved_real_marker = STSList(
        data=np.array([[1, 2]]),
        alphabet=["*", "A", "encoded missing"],
        nr="*",
        void="%",
        nr_code=2,
        void_code=None,
    )
    declared = seqdef(unobserved_real_marker, alphabet=["*", "A"])
    assert declared.alphabet == ["*", "A", "<NR>"]
    assert declared.nr_code == 2
    np.testing.assert_array_equal(declared.data, np.array([[1, 2]]))


def test_seqdef_explicit_missing_reprocesses_code_aware_marker_cells():
    source = seqdef(np.array([["A", None]], dtype=object), right=None)

    reprocessed = seqdef(
        source,
        alphabet=["A"],
        missing="*",
        right="DEL",
    )

    assert reprocessed.alphabet == ["A", "%"]
    assert reprocessed.nr_code is None
    assert reprocessed.void_code == 1
    np.testing.assert_array_equal(reprocessed.data, np.array([[0, 1]]))

    loose = STSList(
        data=np.array([[0, 1]]),
        alphabet=["A", "encoded missing"],
        nr="canonical missing",
        nr_code=1,
        void="%",
        void_code=None,
    )
    canonical = seqdef(
        loose,
        alphabet=["A"],
        missing="canonical missing",
        right="DEL",
    )
    assert canonical.alphabet == ["A", "%"]
    np.testing.assert_array_equal(canonical.data, np.array([[0, 1]]))


def test_seqdef_preserves_declared_marker_labels_and_colors():
    x = seqdef(
        np.array([["A", "*", "%"]], dtype=object),
        alphabet=["A", "*", "%"],
        labels=["Alpha", "Missing", "Void"],
        cpal=["red", "gray", "blue"],
    )

    assert x.labels == ["Alpha", "Missing", "Void"]
    assert x.cpal == ["red", "gray", "blue"]

    with_policy = seqdef(
        np.array([[None, "A"]], dtype=object),
        alphabet=["A", "*"],
        left="L",
        right=None,
        labels=["Alpha", "Missing", "Left"],
        cpal=["red", "gray", "green"],
    )
    assert with_policy.alphabet == ["A", "L", "*"]
    assert with_policy.labels == ["Alpha", "Left", "Missing"]
    assert with_policy.cpal == ["red", "green", "gray"]


def test_seqdef_rejects_marker_policy_collisions_and_unknown_sts_options():
    data = np.array([["A", None]], dtype=object)

    with pytest.raises(ValueError, match="distinct character markers"):
        seqdef(data, nr="%", void="%")
    with pytest.raises(ValueError, match="differ from the 'nr' marker"):
        seqdef(data, nr="*", left="*")
    with pytest.raises(TypeError, match="unsupported STS seqdef arguments"):
        seqdef(data, unexpected=True)
    with pytest.raises(ValueError, match="'stsep' must be a character"):
        seqdef(np.array([["A-B"]], dtype=object), stsep=1)


def test_seqdef_rejects_invalid_labels_colors_column_names_and_ids():
    data = np.array([["A", "B"], ["B", "A"]], dtype=object)

    with pytest.raises(ValueError, match="number of labels"):
        seqdef(data, labels=["only one"])
    with pytest.raises(ValueError, match="number of colors"):
        seqdef(data, cpal=["red"])
    with pytest.raises(ValueError, match="cnames must have length 2"):
        seqdef(data, cnames=["one"])
    with pytest.raises(ValueError, match="invalid 'row.names' length"):
        seqdef(data, id=["one"])
    with pytest.raises(ValueError, match="missing values in 'row.names'"):
        seqdef(data, id=["one", None])
    with pytest.raises(ValueError, match="duplicate 'row.names'"):
        seqdef(data, id=["same", "same"])


def test_seqdef_accepts_scalar_selectors_and_colon_compressed_input():
    selected = seqdef(pd.DataFrame({"first": ["A"], "second": ["B"]}), var="second")
    assert selected.alphabet == ["B"]

    compressed = seqdef(np.array([["A:B:B"]], dtype=object))
    np.testing.assert_array_equal(compressed.data, np.array([[0, 1, 1]]))
    assert compressed.alphabet == ["A", "B"]


def test_seqdef_applies_r_numeric_dataframe_promotion():
    x = seqdef(pd.DataFrame({"logical": [True, False], "integer": [1, 2]}))

    # R: alphabet 0 1 2 (default_alphabet.csv case logical_promotion).
    assert x.alphabet == ["0", "1", "2"]
    np.testing.assert_array_equal(x.data, np.array([[1, 1], [0, 2]]))


def test_seqdef_renders_numeric_categorical_levels_like_r_factors():
    """R's seqdef on ``data.frame(state = factor(c(1, 2, 1)))``: the golden case numeric_factor
    of ``scripts/r_goldens/make_goldens_default_alphabet.R`` (TX-AUD-20261019-088)."""
    golden = pd.read_csv(R_GOLDENS_DIR / "default_alphabet.csv", dtype=str, keep_default_na=False)
    case = golden[golden["case"] == "numeric_factor"].sort_values(
        "position", key=lambda s: s.astype(int)
    )
    r_alphabet = case[case["kind"] == "alphabet"]["value"].tolist()
    r_codes = [int(v) - 1 for v in case[case["kind"] == "codes"]["value"]]  # R's codes are 1-based

    x = seqdef(pd.DataFrame({"state": pd.Categorical([1.0, 2.0, 1.0])}))

    assert x.alphabet == r_alphabet == ["1", "2"]
    np.testing.assert_array_equal(x.data, np.array([[code] for code in r_codes]))


def test_seqdef_treats_string_metadata_arguments_as_scalars():
    x = seqdef(
        np.array([["AB"]], dtype=object),
        alphabet="AB",
        states="State AB",
        labels="Long label",
        cpal="red",
        cnames="T1",
    )

    assert x.alphabet == ["State AB"]
    assert x.labels == ["Long label"]
    assert x.cpal == ["red"]
    assert x.xtlab == ["T1"]


def test_seqdef_rechecks_missing_after_state_mapping_and_policy_replacement():
    mapped = seqdef(
        np.array([["A"]], dtype=object),
        alphabet=["A"],
        states=["M"],
        missing="M",
    )
    assert mapped.alphabet == ["M", "%"]
    np.testing.assert_array_equal(mapped.data, np.array([[1]]))

    replaced = seqdef(
        np.array([["M", "A"]], dtype=object),
        missing="M",
        left="M",
        right=None,
    )
    assert replaced.alphabet == ["A", "M", "*"]
    np.testing.assert_array_equal(replaced.data, np.array([[2, 0]]))


def test_seqdef_resets_missing_token_after_sps_and_spell_normalization():
    sps = seqdef(
        pd.DataFrame([["(A,1)", "(A,1)"]]),
        informat="SPS",
        alphabet=["A"],
        states=["@"],
        missing="@",
    )
    assert sps.alphabet == ["@"]
    np.testing.assert_array_equal(sps.data, np.array([[0, 0]]))

    spell = seqdef(
        pd.DataFrame({"id": ["x"], "begin": [1], "end": [1], "status": ["A"]}),
        informat="SPELL",
        alphabet=["A"],
        states=["@"],
        missing="@",
        limit=1,
    )
    assert spell.alphabet == ["@"]
    np.testing.assert_array_equal(spell.data, np.array([[0]]))


def test_seqdef_preserves_literal_marker_text_from_sps_and_spell_codes():
    sps = seqdef(
        pd.DataFrame([["(*,1)", "(A,1)"]]),
        informat="SPS",
    )
    assert sps.alphabet == ["*", "A"]
    assert sps.nr_code is None

    spell = seqdef(
        pd.DataFrame({"id": ["x"], "begin": [1], "end": [1], "status": ["*"]}),
        informat="SPELL",
        limit=1,
    )
    assert spell.alphabet == ["*"]
    assert spell.nr_code is None


def test_seqdef_preserves_literal_sps_markers_alongside_actual_missing_cells():
    nonresponse = seqdef(
        pd.DataFrame([["(*,1)", "(@,1)", "(A,1)"]]),
        informat="SPS",
        missing="@",
        right=None,
    )
    assert nonresponse.alphabet == ["*", "A", "<NR>"]
    assert nonresponse.nr == "<NR>"
    assert nonresponse.nr_code == 2
    np.testing.assert_array_equal(nonresponse.data, np.array([[0, 2, 1]]))

    trailing_void = seqdef(
        pd.DataFrame([["(%,1)", "(A,1)", "(@,1)"]]),
        informat="SPS",
        missing="@",
    )
    assert trailing_void.alphabet == ["%", "A", "<VOID>"]
    assert trailing_void.void == "<VOID>"
    assert trailing_void.void_code == 2
    np.testing.assert_array_equal(trailing_void.data, np.array([[0, 1, 2]]))
