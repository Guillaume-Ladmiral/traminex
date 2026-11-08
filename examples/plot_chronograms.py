# Example: per-cluster chronograms with representative overlays
import numpy as np

from traminex import cluster_chronograms_with_representatives, seqdef
from traminex.distances import pairwise

rng = np.random.default_rng(42)
mat = rng.integers(0, 6, size=(60, 25))
x = seqdef(mat)
D = pairwise(x.data, method="HAM")
labels = np.repeat([0, 1, 2], repeats=20)

fig, axes = cluster_chronograms_with_representatives(x, labels, D, include_missing=False, cols=3)
fig.suptitle("Chronograms with representative overlays (shared scale)")
fig.subplots_adjust(top=0.78, bottom=0.20, right=0.83)
fig.savefig("chronograms_reps.png", dpi=150)
print("Saved chronograms_reps.png")
