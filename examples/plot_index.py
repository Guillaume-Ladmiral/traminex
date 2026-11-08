# Example: sequence index plot with legend and grouping
import numpy as np

from traminex import cluster_agglomerative, seqdef, sequence_index_plot, sequence_index_plot_grouped
from traminex.distances import pairwise

# Generate random sequences
rng = np.random.default_rng(0)
mat = rng.integers(0, 5, size=(50, 20))
x = seqdef(mat)
# assign a palette (optional)
x.cpal = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

# Basic index plot with legend
# Give all 50 row identifiers enough vertical space to remain readable.
fig, ax = sequence_index_plot(x, legend=True, figsize=(10, 10))
fig.suptitle("Index plot with legend")
fig.savefig("index_plot.png", dpi=150)

# Grouped plot by clustering labels
D = pairwise(x.data, method="HAM")
res = cluster_agglomerative(D, n_clusters=3)
fig2, ax2, order = sequence_index_plot_grouped(x, res.labels_, legend=True)
fig2.suptitle("Grouped index plot (3 clusters)")
fig2.savefig("index_grouped.png", dpi=150)
print("Saved index_plot.png and index_grouped.png")
