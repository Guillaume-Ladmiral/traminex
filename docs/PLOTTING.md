# Plotting and visualization

traminex owns the verified analytical series and the matplotlib renderers that
re-express the TraMineR plots; the semantic parity ledger
([API_PARITY_LEDGER.md](API_PARITY_LEDGER.md)) records each as RE_EXPRESSED or
DIVERGENCE, none as EXACT. The dedicated grammar-of-graphics library is the sibling
[`PROJECT_PY_SEQPLOT`](INDEX.qmd#data-examples-references), package
`py_seqplot`: a plotnine port of ggseqplot that depends on traminex and
reimplements no sequence numerics.

## Matplotlib/reference layer

Install this layer with `traminex[plot]`. It is not pulled into a core install.
The maintained plotting install profile imports from a pure wheel, forces the
noninteractive `Agg` backend, creates a figure, and closes it cleanly.

The `seqplot(x, type=...)` dispatcher covers the main TraMineR routes:

- `i`/`I`: sequence index plots. **Deliberate divergence:** R's `type="i"`
  defaults to the first ten sequences and `type="I"` to all; traminex renders
  all sequences for both spellings, with subsetting/sampling available through
  the renderer's explicit options. A weighted plot draws each row with a height
  equal to its sequence's weight, as R does;
- `d`: state distribution;
- `Ht`/`dH`: transversal entropy;
- `f`: sequence-frequency bands;
- `mt`/`ms`: mean-time and modal-state plots;
- `r`/`rf`: representatives and relative-frequency groups;
- `pc`: parallel coordinate plot (`seqpcplot`) of each sequence's distinct
  successive states, as R routes it.

As R's `seqplot`, every type takes `cpal` (one colour per state; ignored for
`Ht`, as R ignores it), `ylab` (one label, or one per group level) and
`with_legend`: `"auto"` (the default) or `True` draws one state legend for the
figure, `False` none, and `Ht` and `pc` never get one; legends `seqplot` did
not draw (in a figure the caller gave through `ax`) are left alone. The
renderers keep their own legend defaults, called directly or with `main` or
`group` (`sequence_index_plot` draws none unless `legend=True`). `plot_SAMM`
draws one state legend for the figure, as R's `plot.SAMM`, which is
`seqplot(group=)`.

Explicit functions include `sequence_index_plot`, `state_distribution_plot`,
`seqfplot`, `seqmtplot`, `seqmsplot`, `seqplot_tentrop`, `seqplot_tentrop_m`,
`seqrplot`, `seqrfplot`, `seqplot_rf`, `seqsurvplot`, `seqsplot`, `plot_SAMM`,
`plot_seqimplic`, `seqalignplot`, `seqdiffplot`, `seqplotMD`, `seqpcplot`, and
`seqedplot`.

## R4 verified routes

Sequence index plots accept `sampling`, `sample_meth`, `weighted`, and
`random_state` after sorting. `seqplot(type="i")` and `seqplot(type="I")`
forward those controls. Systematic row selection matches TraMineR 2.2-14;
factor-like sort keys retain R's stable level order and place missing values last.
`sampling=None` deliberately keeps all rows for backward compatibility instead
of R's code-level first-ten default. Random sampling uses an isolated NumPy
generator and never changes global RNG state.

`seqplot_tentrop` draws grouped transversal-entropy lines and
`seqplot_tentrop_m` draws the same verified `seqstatd` series across two or more
sequence objects. Both expose source-aligned title, palette, line, limit, time
tick, and legend controls. `seqplot_rf` is the explicit TraMineRextras-name
route over the existing `seqrf` analytical result and supports medoid,
dissimilarity-to-medoid, or paired panels. Custom color sequences and common
matplotlib patch properties such as `alpha` and `linewidth` are applied only at
the rendering boundary and are never forwarded into `seqrf` analytics.

`plot_SAMM(result, type=...)` selects transition rows by origin state and
delegates every supported `seqplot` type. `plot_seqimplic(result, ...)` applies
the R plotting sign reversal and negative clipping to a copy; it never changes
the stored implication array. It retains source sequence time labels and
`legend_prop` reserves the R layout fraction for the legend. `seqsplot(...)`
creates group facets or
`per_state=True` state facets and calls `seqsurv` for every Kaplan–Meier
estimate. Its state/group ordering, missing-state extension, palettes, titles,
layout, axes, messages, and empty panels are structurally verified. Per-call
`xtstep` and `tick_last` override sequence metadata, including final times that
occur only in censored spells. `legend_prop` reserves the R layout fraction and
`use_layout` explicitly toggles Matplotlib constrained layout.

Convenience plots such as `chronogram`, `transition_heatmap`, grouped index
plots, and cluster chronograms remain available. They are Python additions and
should not be confused with an exact R plot contract. A chronogram's rows are
`seqstatd`'s states: the real states, plus the missing state when
`include_missing=True` (a zero row if no value is missing); the void marker
never has a row.

Plot verification is structural: plotted data, ordering, layers, axes, and
labels are checked against R-golden numerical series. Pixel-perfect rendering
is not used as the parity criterion. Renderers canonicalize declared external
marker codes at entry (`STSList._canonicalized`), so both marker
representations of the same data draw identically; `missing_color` styles the
nr state wherever no explicit palette entry covers it. A user colour list gives
one colour per real state, as R's `cpal` (the entry after them colours the nr
state, else `missing_color`; void is white). A shorter list, which R's
documentation does not allow, is recycled over the real states, so each state
keeps one colour; R recycles it cell by cell in its index, frequency and
representative plots, drawing one state in several colours. `seqedplot`'s event
palette follows the same rule: a shorter list is recycled over the events, where
R draws the events beyond its palette in no colour. A plotting call that
raises closes any figure it opened, so failures do not accumulate open
figures; figures the caller supplied are left open. Time ticks follow R's grid:
a tick every `xtstep` positions from the first, plus the last with
`tick_last`. As R's `axis()`, the labels are drawn left to right and one that
would come closer than the width of an "m" to the last label drawn is left
out; the rule runs when the axis is drawn, so the labels shown depend on the
figure's size, as R's depend on the device's.

## Dedicated plotnine layer

Use `py_seqplot` when a composable grammar-of-graphics object is required. It
provides all nine `ggseqplot` renderers (`ggseqdplot`, `ggseqmsplot`,
`ggseqmtplot`, `ggseqtrplot`, `ggseqiplot`, `ggseqfplot`, `ggseqeplot`,
`ggseqrplot`, and `ggseqrfplot`) over traminex's public objects and series.

Numerical or state-sequence changes belong in traminex. Plotnine rendering and
ggseqplot parity belong in `PROJECT_PY_SEQPLOT`.

The supported analytical boundary is frozen in
[`VISUALIZATION_SERIES_CONTRACT.md`](VISUALIZATION_SERIES_CONTRACT.md). It
specifies object shapes, ordering, state/time metadata, weighting,
missing/void denominators, representative/RF inputs, nonmutation, and refusal
behavior. Private underscore-prefixed helpers are not a cross-project API.
