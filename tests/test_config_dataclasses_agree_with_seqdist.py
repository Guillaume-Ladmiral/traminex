"""The configuration dataclasses accept what seqdist accepts and refuse options it would ignore (TX-AUD-20261019-026).

Ruling R6-l: the dataclasses accept seqdist's values and refuse what seqdist refuses. Three departures:
``DHDConfig(variant="MEAN")`` raised where seqdist accepts ``variant="MEAN"``; a ``calibrated_params_path`` was
silently ignored unless a profile was ``"calibrated"``; ``OMConfig(subm="INDELS", trate_config=...)`` dropped the
TRATE settings. An option that would be ignored is now refused, as seqdist refuses TRATE-specific keywords
without ``subm="TRATE"``.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import traminex as tx
from traminex.config_dataclasses import DHDConfig, OMConfig, TRATEConfig
from traminex.distances import pairwise

X = tx.seqdef(np.array([list("AABBC"), list("ABBCC"), list("CCBAA"), list("BBCCA")], dtype=object))


def test_dhd_config_accepts_the_variant_spellings_seqdist_accepts() -> None:
    config = DHDConfig(variant="MEAN")
    assert config.as_dict()["variant"] == "mean"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.testing.assert_array_equal(
            tx.seqdist(X, method="DHD", **config.as_dict()),
            tx.seqdist(X, method="DHD", variant="MEAN"),
        )


@pytest.mark.parametrize(
    "call",
    [
        lambda: tx.seqdist(X, method="OM", subm="TRATE", calibrated_params_path="cfg.json"),
        lambda: tx.seqdist(X, method="DHD", calibrated_params_path="cfg.json"),
        lambda: pairwise(
            np.asarray(X.data), method="OM", subm="TRATE", calibrated_params_path="cfg.json"
        ),
    ],
    ids=["seqdist-om-trate", "seqdist-dhd", "pairwise-om-trate"],
)
def test_a_calibration_path_without_a_calibrated_profile_is_refused(call) -> None:
    with pytest.raises(TypeError, match="calibrated"):
        call()


def test_trate_config_with_another_cost_method_is_refused() -> None:
    with pytest.raises(TypeError, match="trate_config requires subm='TRATE'"):
        OMConfig(subm="INDELS", trate_config=TRATEConfig(smoothing=1.0))
    assert (
        "trate_smoothing"
        in OMConfig(subm="TRATE", trate_config=TRATEConfig(smoothing=1e-6)).as_dict()
    )
