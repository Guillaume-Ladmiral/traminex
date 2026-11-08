import json
import warnings

import numpy as np
import pytest

import traminex as tx
from traminex.config import load_calibrated_params


def test_load_calibrated_params_does_not_search_cwd_implicitly(tmp_path, monkeypatch):
    cfg = {"trate": {"mode": "row", "normalize": "clip_0_2", "smoothing": 1.0}}
    p = tmp_path / "calibrated_params.json"
    p.write_text(json.dumps(cfg), encoding="utf-8")

    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("TRAMINEX_CALIBRATED_PARAMS", raising=False)

    # Deterministic behavior: no implicit cwd lookup.
    assert load_calibrated_params() == {}


def test_load_calibrated_params_explicit_path_and_env(tmp_path, monkeypatch):
    cfg = {"dhd": {"variant": "sum", "smooth_window": None}}
    p = tmp_path / "params.json"
    p.write_text(json.dumps(cfg), encoding="utf-8")

    assert load_calibrated_params(paths=str(p)) == cfg
    assert load_calibrated_params(paths=[str(p)]) == cfg

    monkeypatch.setenv("TRAMINEX_CALIBRATED_PARAMS", str(p))
    assert load_calibrated_params() == cfg


def test_load_calibrated_params_invalid_json_raises(tmp_path):
    p = tmp_path / "bad.json"
    p.write_text("{invalid json", encoding="utf-8")

    with pytest.raises(ValueError, match="Invalid JSON"):
        load_calibrated_params(paths=str(p))


def test_load_calibrated_params_requires_top_level_object(tmp_path):
    p = tmp_path / "list.json"
    p.write_text(json.dumps([1, 2, 3]), encoding="utf-8")

    with pytest.raises(ValueError, match="JSON object"):
        load_calibrated_params(paths=str(p))


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(r"ignore:\s*Non\-symmetric\ substitution\ cost\ matrix!:UserWarning")
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ at\ least,\ one\ substitution\ cost\ doesn't\ respect\ the\ triangle\ inequality\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
def test_seqdist_trate_calibration_requires_explicit_path(tmp_path, monkeypatch):
    arr = np.array(
        [
            ["A", "A", "B", "B"],
            ["A", "B", "B", "A"],
            ["B", "B", "A", "A"],
        ],
        dtype=object,
    )
    x = tx.seqdef(arr)

    cfg = {"trate": {"mode": "row", "normalize": "clip_0_2", "smoothing": 1.0}}
    p = tmp_path / "calibrated_params.json"
    p.write_text(json.dumps(cfg), encoding="utf-8")

    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("TRAMINEX_CALIBRATED_PARAMS", raising=False)

    d_default = tx.seqdist(x, method="OM", subm="TRATE")
    # No implicit lookup from the cwd file, and no silent fallback to the
    # defaults either (TX-AUD-20261010-032).
    with pytest.raises(ValueError, match="no calibration file loaded"):
        tx.seqdist(x, method="OM", subm="TRATE", trate_profile="calibrated")
    d_cal_path = tx.seqdist(
        x,
        method="OM",
        subm="TRATE",
        trate_profile="calibrated",
        calibrated_params_path=str(p),
    )

    # Explicit path should change behavior.
    assert not np.allclose(d_default, d_cal_path)


def test_seqdist_dhd_calibration_requires_explicit_path(tmp_path, monkeypatch):
    arr = np.array(
        [
            ["A", "A", "B", "B"],
            ["A", "B", "B", "A"],
            ["B", "B", "A", "A"],
        ],
        dtype=object,
    )
    x = tx.seqdef(arr)

    cfg = {"dhd": {"variant": "sum", "smooth_window": None}}
    p = tmp_path / "calibrated_params.json"
    p.write_text(json.dumps(cfg), encoding="utf-8")

    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("TRAMINEX_CALIBRATED_PARAMS", raising=False)

    # Calibrated without a file used to switch silently to the legacy DHD with
    # variant "mean" (TX-AUD-20261010-032); it is now refused.
    d_sum = tx.seqdist(x, method="DHD", variant="sum")
    with pytest.raises(ValueError, match="no calibration file loaded"):
        tx.seqdist(x, method="DHD", dhd_profile="calibrated")
    d_cal_path = tx.seqdist(
        x, method="DHD", dhd_profile="calibrated", calibrated_params_path=str(p)
    )

    assert np.allclose(d_sum, d_cal_path)


def test_invalid_calibrated_profile_values_raise():
    arr = np.array(
        [
            ["A", "B", "A"],
            ["B", "A", "B"],
        ],
        dtype=object,
    )
    x = tx.seqdef(arr)

    with pytest.raises(ValueError):
        tx.seqdist(x, method="OM", subm="TRATE", trate_profile="unknown")
    with pytest.raises(ValueError):
        tx.seqdist(x, method="DHD", dhd_profile="unknown")


# Audit 2026-09-22 S2-17: a single pathlib.Path raised "'PosixPath' object is
# not iterable"; a negative or NaN trate.smoothing in the file surfaced as the
# unrelated "OM indel cost must be strictly positive" (or "substitution costs
# must be finite"); OMConfig accepted booleans as costs (subm=True became 1).


def test_load_calibrated_params_accepts_a_single_path_object(tmp_path):
    cfg = {"dhd": {"variant": "sum"}}
    p = tmp_path / "params.json"
    p.write_text(json.dumps(cfg), encoding="utf-8")

    assert load_calibrated_params(paths=p) == cfg
    assert load_calibrated_params(paths=[p]) == cfg


@pytest.mark.parametrize("smoothing", ["-1", "NaN", "Infinity", '"x"', "true"])
def test_invalid_calibrated_trate_smoothing_is_refused_by_name(tmp_path, smoothing):
    p = tmp_path / "params.json"
    p.write_text('{"trate": {"smoothing": ' + smoothing + "}}", encoding="utf-8")
    x = tx.seqdef(np.array([["A", "B", "A"], ["B", "A", "B"]], dtype=object))

    with pytest.raises(ValueError, match="trate.smoothing .* finite non-negative number"):
        load_calibrated_params(paths=str(p))
    with pytest.raises(ValueError, match="trate.smoothing .* finite non-negative number"):
        tx.seqdist(
            x,
            method="OM",
            subm="TRATE",
            trate_profile="calibrated",
            calibrated_params_path=str(p),
        )


@pytest.mark.parametrize("field", ["indel", "subm"])
def test_om_config_refuses_boolean_costs(field):
    from traminex.config_dataclasses import OMConfig

    with pytest.raises(TypeError, match=f"{field} must be a number"):
        OMConfig(**{field: True})


# TX-AUD-20261010-032: a calibrated profile is an explicit request. Without a
# file it silently used the defaults (TRATE) or the legacy DHD with variant
# "mean"; a file without the profile's section or with a misspelt key was
# ignored the same way; and the two profiles let explicit arguments and file
# values win in opposite orders.

ROWS = np.array(
    [["A", "A", "B", "B", "C"], ["A", "B", "B", "C", "C"], ["C", "B", "A", "A", "B"]],
    dtype=object,
)


def _write(tmp_path, payload, name="calibrated.json"):
    path = tmp_path / name
    path.write_text(json.dumps(payload), encoding="utf-8")
    return str(path)


@pytest.mark.parametrize("paths", [None, "missing.json"])
def test_calibrated_trate_without_a_file_raises(tmp_path, monkeypatch, paths):
    monkeypatch.delenv("TRAMINEX_CALIBRATED_PARAMS", raising=False)
    options = {} if paths is None else {"calibrated_params_path": str(tmp_path / paths)}
    with pytest.raises(ValueError, match="calibrated trate profile requested but no calibration"):
        tx.seqdist(
            tx.seqdef(ROWS), method="OM", subm="TRATE", trate_profile="calibrated", **options
        )


def test_calibrated_dhd_without_a_file_raises(tmp_path, monkeypatch):
    monkeypatch.delenv("TRAMINEX_CALIBRATED_PARAMS", raising=False)
    with pytest.raises(ValueError, match="calibrated dhd profile requested but no calibration"):
        tx.seqdist(
            tx.seqdef(ROWS),
            method="DHD",
            dhd_profile="calibrated",
            calibrated_params_path=str(tmp_path / "missing.json"),
        )


@pytest.mark.parametrize(
    ("method", "profile", "payload"),
    [
        ("OM", {"subm": "TRATE", "trate_profile": "calibrated"}, {"dhd": {"variant": "sum"}}),
        ("DHD", {"dhd_profile": "calibrated"}, {"trate": {"smoothing": 1.0}}),
    ],
)
def test_calibrated_profile_without_its_section_raises(tmp_path, method, profile, payload):
    path = _write(tmp_path, payload)
    with pytest.raises(ValueError, match="has no '(trate|dhd)' section"):
        tx.seqdist(tx.seqdef(ROWS), method=method, calibrated_params_path=path, **profile)


def test_calibrated_profile_unknown_key_raises(tmp_path):
    path = _write(tmp_path, {"trate": {"smothing": 1.0}})
    with pytest.raises(ValueError, match=r"unknown key\(s\) \['smothing'\]"):
        tx.seqdist(
            tx.seqdef(ROWS),
            method="OM",
            subm="TRATE",
            trate_profile="calibrated",
            calibrated_params_path=path,
        )


def test_explicit_argument_overrides_file_value(tmp_path):
    x = tx.seqdef(ROWS)
    smoothed = _write(tmp_path, {"trate": {"smoothing": 1.0, "mode": "row"}}, "smoothed.json")
    unsmoothed = _write(tmp_path, {"trate": {"smoothing": 0.0, "mode": "row"}}, "plain.json")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        explicit = tx.seqdist(
            x,
            method="OM",
            subm="TRATE",
            trate_profile="calibrated",
            calibrated_params_path=smoothed,
            trate_smoothing=0.0,
        )
        reference = tx.seqdist(
            x,
            method="OM",
            subm="TRATE",
            trate_profile="calibrated",
            calibrated_params_path=unsmoothed,
        )
        smoothed_only = tx.seqdist(
            x,
            method="OM",
            subm="TRATE",
            trate_profile="calibrated",
            calibrated_params_path=smoothed,
        )
    np.testing.assert_array_equal(explicit, reference)
    assert not np.array_equal(explicit, smoothed_only)

    dhd_file = _write(tmp_path, {"dhd": {"variant": "sum"}}, "dhd.json")
    np.testing.assert_array_equal(
        tx.seqdist(
            x,
            method="DHD",
            dhd_profile="calibrated",
            calibrated_params_path=dhd_file,
            variant="mean",
        ),
        tx.seqdist(x, method="DHD", variant="mean"),
    )


def test_valid_file_reproduces_the_values_it_names(tmp_path):
    x = tx.seqdef(ROWS)
    path = _write(tmp_path, {"trate": {"smoothing": 1.0, "mode": "row", "normalize": "clip_0_2"}})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from_file = tx.seqdist(
            x, method="OM", subm="TRATE", trate_profile="calibrated", calibrated_params_path=path
        )
        named = tx.seqdist(
            x,
            method="OM",
            subm="TRATE",
            trate_smoothing=1.0,
            trate_mode="row",
            trate_normalize="clip_0_2",
        )
    np.testing.assert_array_equal(from_file, named)
