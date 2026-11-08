import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex.config_dataclasses import DHDConfig, DURConfig


def test_seqdef_from_long_rejects_unknown_kwargs():
    long = pd.DataFrame(
        {
            "id": ["p1", "p1"],
            "time": [1, 2],
            "state": ["A", "B"],
        }
    )

    with pytest.raises(TypeError, match="Unsupported kwargs"):
        tx.seqdef_from_long(long, fill="_")


def test_tse_to_sts_rejects_spell_specific_kwargs():
    long = pd.DataFrame(
        {
            "id": ["p1", "p1"],
            "time": [1, 2],
            "state": ["A", "B"],
        }
    )

    with pytest.raises(TypeError, match="Unsupported kwargs"):
        tx.TSE_to_STS(long, id_col="id", time_col="time", state_col="state", fill="_")


def test_dissdomassoc_accepts_weights_and_requires_two_domains():
    # RP4 T14: dissdomassoc is now the R-faithful inter-domain association
    # (weighted Pearson/Spearman + Cronbach) and (re-)accepts a weights kwarg;
    # it needs at least two domains.
    D1 = np.array([[0.0, 1.0, 2.0], [1.0, 0.0, 1.5], [2.0, 1.5, 0.0]], dtype=float)
    D2 = np.array([[0.0, 0.5, 1.0], [0.5, 0.0, 2.0], [1.0, 2.0, 0.0]], dtype=float)
    out = tx.dissdomassoc([D1, D2], weights=np.array([1.0, 2.0, 1.0], dtype=float))
    assert "Pearson" in out

    with pytest.raises(ValueError):
        tx.dissdomassoc([D1])


@pytest.mark.parametrize("variant", ["traminer", "tmean"])
def test_dhd_config_rejects_removed_variants(variant):
    with pytest.raises(ValueError, match="variant"):
        DHDConfig(variant=variant)


def test_dur_config_removed_eps_parameter():
    cfg = DURConfig(metric="chi2")
    assert cfg.as_dict() == {"metric": "chi2"}

    with pytest.raises(TypeError):
        DURConfig(metric="chi2", eps=1e-6)
