"""Optional smoke test: seqdef builds sequences from TraMineR's actcal data.

It reads ``DATA/TraMineR/actcal.RData``, a local copy not in the distribution, and
checks the shape of the result only: it compares nothing with R's output, so it is
not an R parity check, and no gate runs it (TX-AUD-20261019-077).
"""

import os
import pathlib

import pytest

pytestmark = pytest.mark.skipif(
    os.environ.get("TRAMINEX_PARITY") != "1",
    reason="Set TRAMINEX_PARITY=1 to run the actcal seqdef smoke test (no comparison with R)",
)


def test_load_rdata_and_seqdef_runs():
    pyreadr = pytest.importorskip("pyreadr")
    root = pathlib.Path(__file__).resolve().parents[1]
    rdata_path = root / "DATA" / "TraMineR" / "actcal.RData"
    if not rdata_path.exists():
        pytest.skip("TraMineR demo data not found")
    res = pyreadr.read_r(str(rdata_path))
    assert len(res) > 0
    # Take the first object which should be a DataFrame-like
    df = next(iter(res.values()))
    # Use first few columns as states (stringify) and a small subset of rows
    import pandas as pd

    if not isinstance(df, (pd.DataFrame,)):
        pytest.skip("Unexpected RData content type")
    sub = df.iloc[:20, : min(10, df.shape[1])].astype(str)
    # Lazy import to avoid circular
    from traminex import seqdef

    x = seqdef(sub)
    # Basic sanity
    assert x.nseq == sub.shape[0]
    assert x.seqlen == sub.shape[1]
