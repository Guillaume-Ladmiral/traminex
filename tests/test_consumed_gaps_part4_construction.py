"""Consumed-component gaps, part 4: sequence construction and conversion.

``seqdef``'s input and column selection, ``seqformat``'s rarely used options, R's string and
number rendering, ICU collation failures and the lazy exports. See the supervision repository's
``docs/TRAMINEX_COMPONENT_GAPS.md`` (part 4, re-measured on 2026-11-08 at ``2ed89d77``).
R's values come from ``scripts/r_goldens/make_goldens_gap_part4.R`` (R 4.6.2, TraMineR 2.2-14).
"""

from __future__ import annotations

import csv
import importlib
import types
import warnings

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDENS_DIR

import traminex
import traminex as tx
from traminex import _categorical_order, _rformat
from traminex.sequences import STSList

seqformat_module = importlib.import_module("traminex.seqformat")


def _read_golden(name: str) -> dict[tuple[str, str], str]:
    with (R_GOLDENS_DIR / name).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(line for line in handle if not line.startswith("#")))
    return {(case, key): value for case, key, value in rows[1:]}


#: R's values, from scripts/r_goldens/make_goldens_gap_part4.R.
R = _read_golden("gap_part4.csv")

FRAME = pd.DataFrame({"t1": ["A", "B", "A"], "t2": ["B", "B", "A"], "t3": ["C", "A", "B"]})
WIDE = pd.DataFrame({"id": [1, 2], "t1": ["A", "B"], "t2": ["B", "B"], "t3": ["C", "A"]})


def _same(left: STSList, right: STSList) -> bool:
    return np.array_equal(left.data, right.data) and list(left.alphabet) == list(right.alphabet)


# --- seqdef: what a table of states is -----------------------------------------------------------


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(["AB", "BA"], id="flat list of strings"),
        pytest.param(np.zeros((2, 2, 2)), id="three-dimensional array"),
        pytest.param(
            [
                np.array([["A", "B"], ["A", "B"]], dtype=object),
                np.array([["A", "B", "C"], ["A", "B", "C"]], dtype=object),
            ],
            id="rows given as matrices",
        ),
    ],
)
def test_seqdef_refuses_what_is_not_a_two_dimensional_table(data):
    with pytest.raises(ValueError, match="seqdef expects a 2D table of states"):
        tx.seqdef(data)


# --- seqdef: var= column selection ---------------------------------------------------------------


def test_seqdef_reads_a_single_integer_selector_as_a_one_element_list():
    assert _same(tx.seqdef(FRAME, var=1), tx.seqdef(FRAME, var=[1]))


@pytest.mark.parametrize("var", [["t2", "t3"], [1, -1]])
def test_selecting_from_an_existing_object_equals_selecting_from_its_table(var):
    named = tx.seqdef(FRAME, cnames=["t1", "t2", "t3"])
    assert _same(tx.seqdef(named, var=var), tx.seqdef(FRAME, var=var))


def test_selecting_from_an_array_equals_selecting_from_a_frame():
    array = FRAME.to_numpy(dtype=object)
    assert _same(tx.seqdef(array, var=[0, 2]), tx.seqdef(FRAME, var=[0, 2]))


@pytest.mark.parametrize(
    ("data", "var"),
    [
        pytest.param("named", ["zz"], id="object, unknown name"),
        pytest.param("named", [9], id="object, position out of range"),
        pytest.param("plain", ["t1"], id="object without time labels, by name"),
        pytest.param("frame", ["zz"], id="frame, unknown name"),
        pytest.param("frame", [9], id="frame, position out of range"),
        pytest.param("array", ["t1"], id="array, by name"),
        pytest.param("array", [True], id="array, boolean"),
        pytest.param("array", [3], id="array, position out of range"),
    ],
)
def test_seqdef_refuses_undefined_columns(data, var):
    sources = {
        "named": tx.seqdef(FRAME, cnames=["t1", "t2", "t3"]),
        "plain": tx.seqdef(FRAME),
        "frame": FRAME,
        "array": FRAME.to_numpy(dtype=object),
    }
    with pytest.raises(ValueError, match="undefined columns selected"):
        tx.seqdef(sources[data], var=var)


# --- seqdef: ids, informat, markers of an existing object -----------------------------------------


def test_one_sequence_takes_a_text_id():
    assert list(tx.seqdef(pd.DataFrame({"a": ["A"], "b": ["B"]}), id="s1").ids) == ["s1"]


def test_seqdef_refuses_spell_formats_for_an_existing_object_or_a_nested_list():
    with pytest.raises(ValueError, match="SPS/SPELL input must be a data frame or matrix"):
        tx.seqdef(tx.seqdef(FRAME), informat="SPS")
    with pytest.raises(ValueError, match="must be one-dimensional and contain only strings"):
        tx.seqdef([["(A,2)", "(B,1)"], ["(B,3)", None]], informat="SPS")


VOIDED = pd.DataFrame({"a": ["A", "B", "A"], "b": ["B", None, "A"], "c": [None, None, "B"]})


def test_void_cells_named_as_missing_come_back_as_void():
    original = tx.seqdef(VOIDED)
    reread = tx.seqdef(original, missing="%")
    assert _same(reread, original)
    assert reread.void_code == original.void_code


def test_an_alphabet_naming_the_void_marker_keeps_it_a_marker():
    original = tx.seqdef(VOIDED)
    declared = tx.seqdef(original, alphabet=["A", "B", "%"])
    assert _same(declared, original)
    assert declared.void == "%" and declared.void_code == original.void_code


def test_a_marker_outside_the_alphabet_gets_a_free_spelling():
    obj = STSList(
        data=np.array([[0, 1, 2, 4]]), alphabet=["*", "<NR>", "A", "?"], nr="*", nr_code=4
    )
    canonical = obj._canonicalized()
    assert canonical.nr == "<NR:1>"
    assert list(canonical.alphabet) == ["*", "<NR>", "A", "?", "<NR:1>"]


# --- seqformat ----------------------------------------------------------------------------------


def test_sps_output_writes_trailing_missing_as_the_right_state():
    df = pd.DataFrame({"t1": ["A", "B"], "t2": ["B", "B"], "t3": ["C", None]})
    assert " | ".join(tx.seqformat(df, from_="STS", to="SPS", right="X")) == R["sps_right_x", "sps"]


def test_complex_and_logical_cells_are_written_as_r_writes_them():
    complex_cells = pd.DataFrame({"a": [1 + 1j, 2 + 0j], "b": [1 + 1j, 1 + 1j]})
    assert " | ".join(tx.seqformat(complex_cells, from_="STS", to="SPS")) == R["sps_complex", "sps"]
    logical = pd.DataFrame({"a": [True, False], "b": [True, True]})
    assert " | ".join(tx.seqformat(logical, from_="STS", to="SPS")) == R["sps_logical", "sps"]


COMPLEX = [
    1 + 1j,
    -1 - 2.5j,
    0.1 + 0j,
    1e10 + 1j,
    0 - 1j,
    complex(1, -0.0),
    complex(np.inf, -np.inf),
    complex(np.nan, 1),
    complex(1 / 3, -2 / 3),
]


def test_complex_states_render_as_r_as_character():
    expected = R["complex_text", "text"].split(" | ")
    assert [seqformat_module._r_state_text(value) for value in COMPLEX] == expected


def test_complex_states_stay_distinct_states():
    # R's seqdef keeps only the real parts, so its alphabet repeats a state (an R defect).
    cells = pd.DataFrame({"a": [2 + 0j, 1 + 1j, 1 - 1j], "b": [1 + 1j, 1 + 1j, 2 + 0j]})
    r_alphabet = R["seqdef_complex", "alphabet"].split()
    assert len(set(r_alphabet)) < len(r_alphabet)
    assert sorted(tx.seqdef(cells).alphabet) == ["1+1i", "1-1i", "2+0i"]


def test_an_empty_separator_splits_into_characters_and_a_bad_one_is_refused():
    split = tx.seqdef(pd.DataFrame({"s": ["ABC", "BBA"]}), stsep="")
    assert [split.alphabet[code] for code in split.data[0]] == R["strsplit_empty", "parts"].split()
    assert split.data.tolist() == [[0, 1, 2], [1, 1, 0]]
    with pytest.raises(ValueError, match="unterminated subpattern"):
        tx.seqformat(pd.DataFrame({"s": ["A(B", "B(A"]}), from_="STS", to="STS", stsep="(")


def _tevent():
    return tx.seqetm(tx.seqdef(WIDE[["t1", "t2", "t3"]]), method="transition")


def test_tse_ids_from_none_a_column_name_or_a_one_name_list():
    with pytest.warns(
        UserWarning, match="replacing original IDs in the output by the sequence indexes"
    ):
        by_index = tx.seqformat(
            WIDE, from_="STS", to="TSE", var=[2, 3, 4], tevent=_tevent(), id=None
        )
    assert sorted(set(by_index["id"])) == [1, 2]
    by_name = tx.seqformat(WIDE, from_="STS", to="TSE", var=[2, 3, 4], tevent=_tevent(), id="id")
    by_list = tx.seqformat(WIDE, from_="STS", to="TSE", var=[2, 3, 4], tevent=_tevent(), id=["id"])
    assert by_name.equals(by_list)


@pytest.mark.parametrize(
    ("id_value", "message"),
    [
        ([1.5], "'id' invalid type or unknown value"),
        (1.5, "'id' invalid type or unknown value"),
        (["a", "b", "c"], "'id' must contain 2 IDs"),
    ],
)
def test_tse_refuses_unusable_ids(id_value, message):
    with pytest.raises(ValueError, match=message):
        tx.seqformat(WIDE, from_="STS", to="TSE", var=[2, 3, 4], tevent=_tevent(), id=id_value)


def test_a_compressed_column_gives_the_same_events_as_the_wide_table():
    tevent = _tevent()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        compressed = tx.seqformat(
            pd.DataFrame({"seq": ["A-B-C", "B-B-A"]}),
            from_="STS",
            to="TSE",
            stsep="-",
            tevent=tevent,
            id=None,
        )
        wide = tx.seqformat(WIDE[["t1", "t2", "t3"]], from_="STS", to="TSE", tevent=tevent, id=None)
    assert compressed.equals(wide)


def test_tevent_cells_that_are_missing_or_numbers():
    numeric = _tevent().astype(object)
    numeric.loc["A", "B"] = 7
    missing = numeric.copy()
    missing.loc["A", "B"] = None
    with_number = tx.seqformat(WIDE, from_="STS", to="TSE", var=[2, 3, 4], tevent=numeric, id="id")
    without = tx.seqformat(WIDE, from_="STS", to="TSE", var=[2, 3, 4], tevent=missing, id="id")
    assert with_number["event"].tolist() == ["A", "7", "B>C", "B", "B>A"]
    assert without["event"].tolist() == ["A", "B>C", "B", "B>A"]


def test_an_object_without_positions_has_no_events():
    empty = STSList(data=np.zeros((2, 0), dtype=int), alphabet=["A", "B"])
    tevent = pd.DataFrame([["A", "A>B"], ["B>A", "B"]], index=["A", "B"], columns=["A", "B"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        events = tx.seqformat(empty, from_="STS", to="TSE", tevent=tevent, id=None)
    assert len(events) == 0


def test_missing_codes_given_as_a_list_are_accepted():
    obj = tx.seqdef(FRAME)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        assert _same(tx.seqformat(obj, from_="STS", to="STS", missing=["*", "?"]), obj)


@pytest.mark.parametrize("var", [1.5, [1, "a"]])
def test_seqformat_refuses_mixed_or_fractional_selectors(var):
    with pytest.raises(ValueError, match="'var' must be a vector of positive integers or strings"):
        tx.seqformat(WIDE, from_="STS", to="SPS", var=var)


def test_seqformat_selects_an_objects_columns_by_label_and_refuses_unknown_ones():
    named = tx.seqdef(WIDE[["t1", "t2", "t3"]], cnames=["t1", "t2", "t3"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        assert tx.seqformat(named, from_="STS", to="SPS", var=["t2", "t3"]) == [
            "(B,1)-(C,1)",
            "(B,1)-(A,1)",
        ]
    for source, var in (
        (named, ["zz"]),
        (tx.seqdef(WIDE[["t1", "t2", "t3"]]), ["t2"]),
        (named, [9]),
    ):
        with pytest.raises(ValueError, match="undefined columns selected"):
            tx.seqformat(source, from_="STS", to="SPS", var=var)


def test_internal_helpers_guard_their_inputs():
    class Unconvertible:
        def __array__(self, *args, **kwargs):
            raise ValueError("no array form")

    assert seqformat_module._cell_is_na(Unconvertible()) is False
    with pytest.raises(ValueError, match="normalized STS input must be two-dimensional"):
        seqformat_module._materialize_sts(
            np.array(["A", "B"], dtype=object), ids=None, columns=["a", "b"]
        )


# --- R's text and number reading ------------------------------------------------------------------


SPELLINGS = ["nan", "NaN", " nan ", "nanx", "0x1G", "0x10", "0x1p5000", "-0x1p5000"]


def test_r_string_to_double_reads_as_r_as_numeric():
    for text, r_value in zip(SPELLINGS, R["as_numeric", "values"].split(), strict=True):
        value = _categorical_order.r_string_to_double(text)
        if r_value == "NA":
            assert value is None, text
        elif r_value == "NaN":
            assert value is not None and np.isnan(value), text
        else:
            assert value == float(r_value), text


@pytest.mark.parametrize("first", ["1", "1.0"])
def test_numerically_equal_states_are_ordered_by_their_text(first):
    other = "1.0" if first == "1" else "1"
    frame = pd.DataFrame({"a": [first, other, "2"], "b": ["2", first, other]})
    assert tx.seqdef(frame).alphabet == ["1", "1.0", "2"]


def test_r_as_character_reads_pandas_missing_values_as_na():
    assert _rformat.r_as_character(pd.NA) is None
    assert _rformat.r_as_character(pd.NaT) is None


# --- ICU collation failures -----------------------------------------------------------------------


def _fake_icu(open_result, keys):
    calls = iter(keys)

    def open_collator(name, status):
        return open_result

    def get_sort_key(collator, source, length, buffer, capacity):
        return next(calls)

    def close_collator(collator):
        return None

    return lambda: (object(), open_collator, get_sort_key, close_collator)


@pytest.mark.parametrize(
    ("open_result", "keys", "message"),
    [
        (0, [], "failed to create the ICU root collator"),
        (1, [0], "ICU returned an invalid character sort key"),
        (1, [4, 3], "ICU returned an inconsistent character sort key"),
    ],
)
def test_icu_failures_are_reported(monkeypatch, open_result, keys, message):
    monkeypatch.setattr(_categorical_order, "_icu_root_collation_api", _fake_icu(open_result, keys))
    with pytest.raises(RuntimeError, match=message):
        _categorical_order.r_character_order(["zeta-gap4", "alpha-gap4"])


# --- lazy exports ---------------------------------------------------------------------------------


def _lazy_name():
    return next(iter(traminex._LAZY_ATTR_TO_MODULE))


def test_a_lazy_export_reraises_an_import_error_that_is_not_an_optional_extra(monkeypatch):
    name = _lazy_name()
    monkeypatch.delitem(traminex.__dict__, name, raising=False)

    def broken(module, package=None):
        raise ImportError("broken install", name="definitely_not_an_extra")

    monkeypatch.setattr(traminex.importlib, "import_module", broken)
    with pytest.raises(ImportError, match="broken install"):
        getattr(traminex, name)


def test_a_lazy_export_missing_from_its_module_is_named(monkeypatch):
    name = _lazy_name()
    monkeypatch.delitem(traminex.__dict__, name, raising=False)
    monkeypatch.setattr(
        traminex.importlib,
        "import_module",
        lambda module, package=None: types.ModuleType("traminex.fake"),
    )
    with pytest.raises(AttributeError, match=f"does not define expected export '{name}'"):
        getattr(traminex, name)
