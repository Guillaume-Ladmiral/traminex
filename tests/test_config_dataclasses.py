import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL
from test_r_goldens_distances import GOLDENS, SEQS

import traminex as tx
from traminex.config_dataclasses import (
    DHDConfig,
    DURConfig,
    LCSLCPConfig,
    OMConfig,
    TRATEConfig,
    default_configs,
)
from traminex.distances.costs import trate_substitution_matrix


def test_trate_config_accepts_runtime_supported_values():
    cfg = TRATEConfig(mode="row", normalize="mean_offdiag")
    kwargs = cfg.as_dict()
    assert kwargs["mode"] == "row"
    assert kwargs["normalize"] == "mean_offdiag"

    cfg2 = TRATEConfig(mode="col", normalize="clip_0_2")
    kwargs2 = cfg2.as_dict()
    assert kwargs2["mode"] == "col"
    assert kwargs2["normalize"] == "clip_0_2"


def test_trate_config_rejects_obsolete_values():
    with pytest.raises(ValueError, match="smoothing must be non-negative"):
        TRATEConfig(smoothing=-1.0)
    with pytest.raises(ValueError):
        TRATEConfig(mode="asymmetric")
    with pytest.raises(ValueError):
        TRATEConfig(mode="mean")
    with pytest.raises(ValueError):
        TRATEConfig(normalize="max")
    with pytest.raises(ValueError):
        TRATEConfig(normalize="range")


def test_trate_config_works_with_trate_substitution_matrix():
    seqmat = np.array(
        [
            [0, 1, 1, 0],
            [1, 1, 0, 0],
            [0, 0, 1, 1],
        ],
        dtype=int,
    )
    cfg = TRATEConfig(mode="symmetric", normalize="none")
    sm = trate_substitution_matrix(seqmat, **cfg.as_dict())
    assert sm.shape == (2, 2)
    assert np.allclose(np.diag(sm), 0.0)


def test_dhd_config_serializes_only_active_options_and_validates_values():
    # No variant by default: seqdist's default DHD is TraMineR's (TX-AUD-20261010-005).
    assert DHDConfig().as_dict() == {"dhd_profile": "default"}
    assert DHDConfig(variant="mean").as_dict() == {"variant": "mean", "dhd_profile": "default"}
    assert DHDConfig(variant="sum", smooth_window=3, profile="calibrated").as_dict() == {
        "variant": "sum",
        "dhd_profile": "calibrated",
        "smooth_window": 3,
    }

    with pytest.raises(ValueError, match="variant must be one of"):
        DHDConfig(variant="legacy")
    with pytest.raises(ValueError, match="smooth_window must be positive"):
        DHDConfig(smooth_window=0)
    for invalid_window in (True, 1.5):
        with pytest.raises(TypeError, match="smooth_window must be an integer"):
            DHDConfig(smooth_window=invalid_window)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="profile must be one of"):
        DHDConfig(profile="implicit")


@pytest.mark.parametrize("config", [DHDConfig(), default_configs()["DHD"]])
def test_default_dhd_config_selects_traminer_dhd(config):
    seqs = tx.seqdef(SEQS, alphabet=["A", "B", "C"])
    np.testing.assert_allclose(
        tx.seqdist(seqs, method="DHD", **config.as_dict()),
        GOLDENS["DHD_default"],
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_calibrated_dhd_config_takes_its_variant_from_the_calibration_file(tmp_path):
    # TX-REP-20261012-002: DHDConfig(profile="calibrated") used to emit
    # variant="mean", which overrode the variant the calibration file names.
    calibration = tmp_path / "calibrated.json"
    calibration.write_text('{"dhd": {"variant": "sum"}}', encoding="utf-8")
    seqs = tx.seqdef(SEQS, alphabet=["A", "B", "C"])
    calibrated = tx.seqdist(
        seqs,
        method="DHD",
        calibrated_params_path=str(calibration),
        **DHDConfig(profile="calibrated").as_dict(),
    )
    np.testing.assert_allclose(calibrated, tx.seqdist(seqs, method="DHD", variant="sum"))


def test_om_config_serializes_optional_and_trate_parameters():
    indels = np.array([1.0, 2.0])
    substitutions = np.zeros((2, 2, 2))
    weights = np.array([0.25, 0.75])
    trate = TRATEConfig(smoothing=0.5, mode="row", normalize="clip_0_2", weights=weights)

    kwargs = OMConfig(
        indel=1.5,
        subm="TRATE",
        norm="max",
        indel_by_pos=indels,
        subm_by_pos=substitutions,
        trate_config=trate,
    ).as_dict()

    assert kwargs["indel"] == 1.5
    assert kwargs["subm"] == "TRATE"
    assert kwargs["norm"] == "max"
    assert kwargs["indel_by_pos"] is indels
    assert kwargs["subm_by_pos"] is substitutions
    assert kwargs["trate_smoothing"] == 0.5
    assert kwargs["trate_mode"] == "row"
    assert kwargs["trate_normalize"] == "clip_0_2"
    assert kwargs["trate_weights"] is weights

    without_weights = OMConfig(subm="trate", trate_config=TRATEConfig()).as_dict()
    assert "trate_weights" not in without_weights

    # A TRATE config with another substitution mode would be ignored, so it is refused, as
    # seqdist refuses TRATE-specific keywords without subm="TRATE" (TX-AUD-20261019-026).
    with pytest.raises(TypeError, match="trate_config requires subm='TRATE'"):
        OMConfig(trate_config=trate)


@pytest.mark.parametrize(
    "factory, message",
    [
        (lambda: OMConfig(indel=-1), "indel must be strictly positive"),
        (lambda: OMConfig(indel=0), "indel must be strictly positive"),
        (lambda: OMConfig(indel=np.nan), "indel must be strictly positive"),
        (lambda: OMConfig(subm=-1), "subm must be non-negative"),
        # seqdist accepts "auto" and "max"; only an unknown spelling is refused
        # (TX-AUD-20261011-043).
        (lambda: OMConfig(norm="bogus"), "norm must be one of"),
        (lambda: LCSLCPConfig(norm="bogus"), "norm must be one of"),
        (lambda: DURConfig(metric="euclidean"), "metric must be one of"),
    ],
)
def test_distance_configs_reject_invalid_values(factory, message):
    with pytest.raises(ValueError, match=message):
        factory()


def test_simple_configs_and_default_registry_are_complete():
    assert LCSLCPConfig().as_dict() == {}
    assert LCSLCPConfig(norm="length").as_dict() == {"norm": "length"}
    assert DURConfig(metric="chi2").as_dict() == {"metric": "chi2"}

    configs = default_configs()
    assert list(configs) == ["OM", "LCS", "LCP", "RLCP", "DHD", "DUR", "TRATE"]
    assert isinstance(configs["OM"], OMConfig)
    assert isinstance(configs["LCS"], LCSLCPConfig)
    assert isinstance(configs["DHD"], DHDConfig)
    assert isinstance(configs["DUR"], DURConfig)
    assert isinstance(configs["TRATE"], TRATEConfig)


def test_trate_smoothing_default_is_r_faithful_zero():
    """TraMineR's seqsubm/seqcost TRATE route applies no smoothing; the
    cycle-8 remediation aligned every entry point on 0.0 and the public
    dataclass route must match (cycle-9 finding N-13)."""
    assert TRATEConfig().smoothing == 0.0
    assert default_configs()["TRATE"].smoothing == 0.0


# --- The dataclasses configure seqdist, so they take its values (TX-AUD-20261011-043). ---

S3 = tx.seqdef(
    np.array([["A", "B", "A"], ["B", "B", "A"], ["A", "A", "B"]], dtype=object), alphabet=["A", "B"]
)


def test_om_config_takes_seqdist_indel_values() -> None:
    for indel in ("auto", [1.0, 2.0], 1.5):
        config = OMConfig(indel=indel)
        np.testing.assert_array_equal(
            tx.seqdist(S3, method="OM", **config.as_dict()),
            tx.seqdist(S3, method="OM", indel=indel, subm=2.0),
        )


def test_om_config_defaults_to_the_automatic_indel_with_unchanged_distances() -> None:
    assert OMConfig().indel == "auto"
    np.testing.assert_array_equal(
        tx.seqdist(S3, method="OM", **OMConfig().as_dict()),
        tx.seqdist(S3, method="OM", indel=1.0, subm=2.0),
    )


def test_om_config_trate_gives_seqdist_trate() -> None:
    np.testing.assert_array_equal(
        tx.seqdist(S3, method="OM", **OMConfig(subm="TRATE").as_dict()),
        tx.seqdist(S3, method="OM", subm="TRATE"),
    )


@pytest.mark.parametrize("norm", ["none", "auto", "maxlength", "gmean", "maxdist", "YujianBo"])
def test_configs_take_seqdist_norms(norm: str) -> None:
    np.testing.assert_array_equal(
        tx.seqdist(S3, method="OM", **OMConfig(norm=norm).as_dict()),
        tx.seqdist(S3, method="OM", indel="auto", subm=2.0, norm=norm),
    )
    np.testing.assert_array_equal(
        tx.seqdist(S3, method="LCS", **LCSLCPConfig(norm=norm).as_dict()),
        tx.seqdist(S3, method="LCS", norm=norm),
    )


@pytest.mark.parametrize(
    "factory, error, field",
    [
        (lambda: TRATEConfig(smoothing="x"), TypeError, "smoothing"),
        (lambda: OMConfig(indel="bogus"), ValueError, "indel"),
        (lambda: OMConfig(indel=[1.0, -1.0]), ValueError, "indel"),
        (lambda: OMConfig(indel=[[1.0]]), ValueError, "indel"),
        (lambda: OMConfig(subm=-np.ones((2, 2))), ValueError, "subm"),
        (lambda: OMConfig(subm=np.ones((2, 3))), ValueError, "subm"),
        (lambda: OMConfig(subm="bogus"), ValueError, "subm"),
        # batch-6 review M-5: float(None) and NumPy's conversion errors leaked.
        (lambda: OMConfig(subm=None), TypeError, "subm"),
        (lambda: OMConfig(subm=object()), TypeError, "subm"),
        (lambda: OMConfig(subm=[["0", "x"], ["x", "0"]]), ValueError, "subm"),
    ],
)
def test_configs_name_the_field_they_refuse(factory, error, field) -> None:
    with pytest.raises(error, match=field):
        factory()
