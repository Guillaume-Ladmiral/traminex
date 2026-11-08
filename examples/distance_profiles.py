# Example: TRATE + OM distances and DHD
import numpy as np

from traminex import seqdef
from traminex.distances import om_trate_pairwise, pairwise

rng = np.random.default_rng(123)
mat = rng.integers(0, 5, size=(30, 15))
x = seqdef(mat)

# TRATE convenience (normalized by length)
D_trate = om_trate_pairwise(
    x.data,
    indel=1.0,
    trate_smoothing=1e-6,
    trate_mode="symmetric",
    trate_normalize="none",
    norm="length",
)
print("TRATE OM matrix shape:", D_trate.shape)

# DHD mean with smoothing
D_dhd = pairwise(x.data, method="DHD", variant="mean", smooth_window=3)
print("DHD matrix shape:", D_dhd.shape)
