"""The deterministic components of dissCompare / seqCompare / seqLRT / seqBIC
compared to live R on a corpus where the options change the answer (W-02).

These four are stochastic at their default ``s=100``, which is why earlier
cycles pinned them only lightly. Two surfaces are exactly comparable anyway:

* ``s=0`` disables sampling completely -- no ``set.seed``, no ``permute``, no
  ``sample`` (``dissCompare.R:198``, ``seqCompare.R:198``). Every returned
  statistic is deterministic and is compared here at ``%.17g``.
* ``seqxcomp`` (``seqCompare.R:312``) is a pure function of its arguments, so
  the test statistic, its p-value, the BIC term and the Bayes factor -- plus
  the three sums of squares they are built from -- are compared directly.

The existing ``group_compare_r3`` goldens mine a *constant* dissimilarity
matrix, deliberately, so that their sampling cases do not depend on
R-versus-NumPy random bytes. The side effect is that their deterministic cases
cannot discriminate the options they name: ``s0_squared_false`` and
``s0_squared_true`` are byte-identical there. This module uses a varied matrix
built from real sequences, and asserts that separation explicitly.

The sampling route's *values* are not cross-language reproducible, so what is
pinned for ``s>0`` is its shape and column labelling, which are. Python-side
seeded reproducibility is proven separately, below.

Generator: ``scripts/r_goldens/make_goldens_group_compare_deterministic.R``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR, load_golden_lines

from traminex import dissCompare, seqBIC, seqCompare, seqdef, seqLRT
from traminex.group_compare import _seqxcomp

# --------------------------------------------------------------------------- #
# Corpus, rebuilt from the fixtures rather than from R's RNG.
# --------------------------------------------------------------------------- #
ALPHABET = load_golden_lines("group_compare_det_alphabet")
WEIGHTS = np.array([float(v) for v in load_golden_lines("group_compare_det_weights")])
GROUP = load_golden_lines("group_compare_det_group")
SETS = load_golden_lines("group_compare_det_sets")
R1 = np.array([int(v) for v in load_golden_lines("group_compare_det_r1")], dtype=int)
R2 = np.array([int(v) for v in load_golden_lines("group_compare_det_r2")], dtype=int)


def _codes() -> np.ndarray:
    return pd.read_csv(R_GOLDENS_DIR / "group_compare_det_codes.csv").to_numpy(dtype=int)


def _diss() -> np.ndarray:
    return pd.read_csv(R_GOLDENS_DIR / "group_compare_det_diss.csv").to_numpy(dtype=float)


def _sequences():
    states = np.asarray(ALPHABET, dtype=object)[_codes()]
    return seqdef(pd.DataFrame(states), alphabet=list(ALPHABET), weights=WEIGHTS)


DISS = _diss()


def test_the_reference_dissimilarity_matrix_is_not_degenerate():
    """The whole point of this generator: earlier goldens used a constant matrix."""
    assert DISS.shape == (18, 18)
    assert len(np.unique(DISS)) > 10


def test_group_sizes_per_set_match_r():
    expected = {}
    for line in (R_GOLDENS_DIR / "group_compare_det_sizes.csv").read_text().splitlines():
        if line.startswith("#"):
            continue
        name, n1, n2 = line.split(",")
        expected[name] = (int(n1), int(n2))

    group = np.asarray(GROUP)
    sets = np.asarray(SETS)
    assert (int(np.sum(group == "g1")), int(np.sum(group == "g2"))) == expected["ALL"]
    for level in ("early", "late"):
        in_set = sets == level
        got = (
            int(np.sum(in_set & (group == "g1"))),
            int(np.sum(in_set & (group == "g2"))),
        )
        assert got == expected[level], level


# --------------------------------------------------------------------------- #
# 1. seqxcomp: the statistic, and the sums of squares behind it.
# --------------------------------------------------------------------------- #
def _worker_goldens():
    rows = []
    for line in (R_GOLDENS_DIR / "group_compare_det_seqxcomp.csv").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split(",")
        rows.append(
            {
                "squared": parts[0] == "TRUE",
                "lrt_pow": int(parts[1]),
                "weight_mode": parts[2],
                "stats": np.array([float(v) for v in parts[3:7]]),
                "ss": np.array([float(v) for v in parts[7:10]]),
                "n0": int(parts[10]),
            }
        )
    return rows


WORKER = _worker_goldens()


def test_the_worker_grid_actually_varies():
    """Guards every case below against passing on a degenerate corpus."""
    assert len(WORKER) == 12
    assert len({tuple(row["stats"]) for row in WORKER}) == 12


@pytest.mark.parametrize(
    "case", WORKER, ids=lambda c: f"sq={c['squared']}-pow={c['lrt_pow']}-{c['weight_mode']}"
)
def test_seqxcomp_matches_r_including_its_sums_of_squares(case):
    weighted = case["weight_mode"] != "none"
    weights = WEIGHTS if weighted else None

    got = _seqxcomp(
        R1,
        R2,
        DISS,
        weights,
        is_lrt=True,
        is_bic=True,
        squared=case["squared"],
        weighted=weighted,
        weight_by="by.group" if case["weight_mode"] == "by.group" else "global",
        lrt_pow=case["lrt_pow"],
    )

    np.testing.assert_allclose(got, case["stats"], rtol=0, atol=R_GOLDEN_ATOL)

    # And the intermediates, so a divergence can be localised rather than
    # merely observed: SS, SS1, SS2 are what the likelihood ratio is built from.
    from traminex.diss import disscenter

    n1, n2 = R1.size, R2.size
    if weighted:
        w1, w2 = WEIGHTS[R1].copy(), WEIGHTS[R2].copy()
        if case["weight_mode"] == "by.group":
            w1 = n1 / w1.sum() * w1
            w2 = n2 / w2.sum() * w2
        wcat = np.concatenate([w1, w2])
    else:
        w1, w2 = np.ones(n1), np.ones(n2)
        wcat = np.ones(n1 + n2)
    both = np.concatenate([R1, R2])
    power = case["lrt_pow"]
    ss = float(
        np.sum(
            wcat
            * np.asarray(
                disscenter(DISS[np.ix_(both, both)], weights=wcat, squared=case["squared"])
            )
            ** power
        )
    )
    ss1 = float(
        np.sum(
            w1
            * np.asarray(disscenter(DISS[np.ix_(R1, R1)], weights=w1, squared=case["squared"]))
            ** power
        )
    )
    ss2 = float(
        np.sum(
            w2
            * np.asarray(disscenter(DISS[np.ix_(R2, R2)], weights=w2, squared=case["squared"]))
            ** power
        )
    )
    np.testing.assert_allclose([ss, ss1, ss2], case["ss"], rtol=0, atol=1e-9)
    assert n1 + n2 == case["n0"]


# --------------------------------------------------------------------------- #
# 2. The full deterministic (s = 0) grid.
# --------------------------------------------------------------------------- #
def _grid() -> dict[str, pd.DataFrame]:
    frame = pd.read_csv(R_GOLDENS_DIR / "group_compare_det_grid.csv", dtype=str)
    out: dict[str, pd.DataFrame] = {}
    for case, block in frame.groupby("case", sort=False):
        pivot = block.pivot(index="row", columns="column", values="value")
        # preserve R's column order rather than pandas' alphabetical sort
        order = list(dict.fromkeys(block["column"]))
        rows = list(dict.fromkeys(block["row"]))
        out[str(case)] = pivot.loc[rows, order].astype(float)
    return out


GRID = _grid()

_SQUARED = {"FALSE": False, "TRUE": True, "LRTonly": "LRTonly"}
_WEIGHTED = {"FALSE": False, "TRUE": True, "by.group": "by.group"}

DISS_CASES = [key for key in GRID if key.startswith("diss|")]
SEQ_CASES = [key for key in GRID if "|method=" in key]


def test_the_grid_covers_the_whole_option_product():
    assert len(DISS_CASES) == 3 * 3 * 3 * 2  # stat x squared x weighted x inc.df


def test_squared_actually_changes_the_answer_on_this_corpus():
    """The r3 goldens cannot show this; their two squared rows are identical."""
    false_case = GRID["diss|stat=all|squared=FALSE|weighted=TRUE|inc.df=FALSE"]
    true_case = GRID["diss|stat=all|squared=TRUE|weighted=TRUE|inc.df=FALSE"]
    assert not np.allclose(false_case.to_numpy(), true_case.to_numpy(), rtol=0, atol=1e-9)


@pytest.mark.parametrize("case", DISS_CASES)
def test_disscompare_s0_matches_r_across_the_option_grid(case):
    _, stat, squared, weighted, inc_df = (part.split("=")[-1] for part in case.split("|"))
    expected = GRID[case]

    got = dissCompare(
        DISS,
        GROUP,
        s=0,
        stat=stat,
        squared=_SQUARED[squared],
        weighted=_WEIGHTED[weighted],
        weights=WEIGHTS,
        inc_df=inc_df == "TRUE",
    )

    assert list(got.columns) == list(expected.columns)
    np.testing.assert_allclose(
        got.to_numpy(dtype=float), expected.to_numpy(), rtol=0, atol=R_GOLDEN_ATOL
    )


@pytest.mark.parametrize("stat", ["all", "LRT", "BIC"])
def test_disscompare_s0_with_sets_matches_r_shape_labels_and_first_set(stat):
    """R's ``set=`` route is right only for the first level -- see below."""
    expected = GRID[f"diss_sets|stat={stat}"]

    got = dissCompare(DISS, GROUP, set=SETS, s=0, stat=stat, weights=WEIGHTS, inc_df=True)

    assert list(got.index) == list(expected.index) == ["early", "late"]
    assert list(got.columns) == list(expected.columns)
    np.testing.assert_allclose(
        got.to_numpy(dtype=float)[0], expected.to_numpy()[0], rtol=0, atol=R_GOLDEN_ATOL
    )


@pytest.mark.parametrize("stat", ["all", "LRT", "BIC"])
@pytest.mark.parametrize("level", ["early", "late"])
def test_each_set_matches_r_run_on_that_set_alone(level, stat):
    """The correct per-set answer, derived from R without the defective path.

    ``dissCompare.R:119`` computes ``which(group[set == lev] == lev.g[1])`` --
    positions *within* the set -- and then indexes the *full* dissimilarity
    matrix with them at line 179. For the first set level those positions
    happen to equal the absolute rows, so R is right by accident; for every
    later level it compares the wrong cases. Calling R on the explicitly
    subsetted matrix takes a single-set path where the defect cannot fire, and
    that is the expectation used here.
    """
    rows = np.flatnonzero(np.asarray(SETS) == level)
    expected = GRID[f"diss_set_alone|set={level}|stat={stat}"]

    got = dissCompare(DISS, GROUP, set=SETS, s=0, stat=stat, weights=WEIGHTS, inc_df=True)
    row = got.loc[[level]]

    reference = dissCompare(
        DISS[np.ix_(rows, rows)],
        [GROUP[i] for i in rows],
        s=0,
        stat=stat,
        weights=WEIGHTS[rows],
        inc_df=True,
    )
    np.testing.assert_allclose(
        reference.to_numpy(dtype=float), expected.to_numpy(), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_allclose(
        row.to_numpy(dtype=float), expected.to_numpy(), rtol=0, atol=R_GOLDEN_ATOL
    )


def test_the_upstream_set_indexing_defect_is_visible_and_not_copied():
    """Pins the divergence itself, so it cannot be "fixed" back into R's answer.

    R agrees with traminex on the first set level and disagrees on the second,
    because only the second one's within-set positions differ from its absolute
    rows. If a later change made traminex reproduce R here, this fails.
    """
    r_values = GRID["diss_sets|stat=all"].to_numpy()
    got = dissCompare(DISS, GROUP, set=SETS, s=0, stat="all", weights=WEIGHTS, inc_df=True)

    np.testing.assert_allclose(got.to_numpy()[0], r_values[0], rtol=0, atol=R_GOLDEN_ATOL)
    assert not np.allclose(got.to_numpy()[1], r_values[1], rtol=0, atol=1e-6)


@pytest.mark.parametrize("case", SEQ_CASES)
def test_sequence_route_s0_matches_r_across_methods_and_weighting(case):
    entry, method_part, weighted_part = case.split("|")
    method = method_part.split("=")[-1]
    weighted = _WEIGHTED[weighted_part.split("=")[-1]]
    expected = GRID[case]

    kwargs: dict[str, object] = {
        "group": GROUP,
        "s": 0,
        "weighted": weighted,
        "method": method,
    }
    if method == "OM":
        # R's `sm` is spelled `subm` here (docs/API_CONSISTENCY.md section 1b).
        kwargs.update(subm="TRATE", indel=1.0)
    caller = {"seqCompare": seqCompare, "seqLRT": seqLRT, "seqBIC": seqBIC}[entry]
    if entry == "seqCompare":
        kwargs["stat"] = "all"

    got = caller(_sequences(), **kwargs)

    assert list(got.columns) == list(expected.columns)
    np.testing.assert_allclose(
        got.to_numpy(dtype=float), expected.to_numpy(), rtol=0, atol=R_GOLDEN_ATOL
    )


def test_two_object_sequence_route_s0_matches_r():
    sequences = _sequences()
    group = np.asarray(GROUP)
    expected = GRID["seqCompare|two_object"]

    got = seqCompare(
        _subset(sequences, group == "g1"),
        _subset(sequences, group == "g2"),
        s=0,
        method="LCS",
        stat="all",
    )

    assert list(got.columns) == list(expected.columns)
    np.testing.assert_allclose(
        got.to_numpy(dtype=float), expected.to_numpy(), rtol=0, atol=R_GOLDEN_ATOL
    )


def test_set_labelled_sequence_route_s0_matches_r():
    expected = GRID["seqCompare|sets"]

    got = seqCompare(
        _sequences(), group=GROUP, set=SETS, s=0, method="LCS", stat="all", inc_df=True
    )

    assert list(got.index) == list(expected.index)
    assert list(got.columns) == list(expected.columns)
    np.testing.assert_allclose(
        got.to_numpy(dtype=float), expected.to_numpy(), rtol=0, atol=R_GOLDEN_ATOL
    )


def _subset(sequences, mask):
    from traminex.group_compare import _subset_sequence

    return _subset_sequence(sequences, np.flatnonzero(mask))


# --------------------------------------------------------------------------- #
# 3. Sampling route: shape and labels are deterministic even though the values
#    are not, so they are compared to R directly.
# --------------------------------------------------------------------------- #
def _shapes():
    out = {}
    for line in (R_GOLDENS_DIR / "group_compare_det_shapes.csv").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        tag, nrow, ncol, columns, rownames = line.split(",", 4)
        out[tag] = (int(nrow), int(ncol), columns.split("|"), rownames.split("|"))
    return out


SHAPES = _shapes()


@pytest.mark.parametrize("bfopt", [None, 1, 2])
@pytest.mark.parametrize("inc_df", [False, True])
def test_sampling_route_shape_and_labels_match_r(bfopt, inc_df):
    tag = f"sample|BFopt={'NULL' if bfopt is None else bfopt}" + ("|inc.df" if inc_df else "")
    nrow, ncol, columns, _ = SHAPES[tag]

    got = dissCompare(DISS, GROUP, s=4, seed=99, BFopt=bfopt, weights=WEIGHTS, inc_df=inc_df)

    assert got.shape == (nrow, ncol)
    assert list(got.columns) == columns


@pytest.mark.parametrize("stat", ["all", "LRT", "BIC"])
@pytest.mark.parametrize("sampled", [True, False])
def test_stat_selection_shape_matches_r_in_both_routes(stat, sampled):
    tag = f"{'sample' if sampled else 's0'}|stat={stat}"
    nrow, ncol, columns, _ = SHAPES[tag]

    got = dissCompare(DISS, GROUP, s=4 if sampled else 0, seed=99, stat=stat, weights=WEIGHTS)

    assert got.shape == (nrow, ncol)
    assert list(got.columns) == columns


def test_sampling_route_with_sets_matches_r_shape_and_row_labels():
    nrow, ncol, columns, rownames = SHAPES["sample|sets"]

    got = dissCompare(DISS, GROUP, set=SETS, s=4, seed=99, weights=WEIGHTS, inc_df=True)

    assert got.shape == (nrow, ncol)
    assert list(got.columns) == columns
    assert list(got.index) == rownames


# --------------------------------------------------------------------------- #
# 4. Seeded reproducibility of the stochastic route -- a Python-side property,
#    since the draws are not cross-language reproducible.
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("caller", [dissCompare], ids=["dissCompare"])
def test_an_integer_seed_reproduces_the_sampled_result_exactly(caller):
    first = caller(DISS, GROUP, s=5, seed=4242, weights=WEIGHTS)
    second = caller(DISS, GROUP, s=5, seed=4242, weights=WEIGHTS)

    pd.testing.assert_frame_equal(first, second, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL)


def test_different_seeds_give_different_sampled_results():
    """Guards the reproducibility test against passing because s is ignored."""
    first = dissCompare(DISS, GROUP, s=5, seed=4242, weights=WEIGHTS)
    second = dissCompare(DISS, GROUP, s=5, seed=99, weights=WEIGHTS)

    assert not np.allclose(first.to_numpy(), second.to_numpy(), rtol=0, atol=1e-12)


def test_sequence_route_is_reproducible_under_an_integer_seed():
    first = seqLRT(_sequences(), group=GROUP, s=5, seed=808, method="LCS")
    second = seqLRT(_sequences(), group=GROUP, s=5, seed=808, method="LCS")

    pd.testing.assert_frame_equal(first, second, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL)


def test_a_supplied_generator_is_consumed_rather_than_reseeded():
    """An explicit Generator must advance, so two calls differ; a seed must not."""
    rng = np.random.default_rng(2024)
    first = dissCompare(DISS, GROUP, s=5, seed=rng, weights=WEIGHTS)
    second = dissCompare(DISS, GROUP, s=5, seed=rng, weights=WEIGHTS)

    assert not np.allclose(first.to_numpy(), second.to_numpy(), rtol=0, atol=1e-12)

    fresh = dissCompare(DISS, GROUP, s=5, seed=np.random.default_rng(2024), weights=WEIGHTS)
    np.testing.assert_allclose(fresh.to_numpy(), first.to_numpy(), rtol=0, atol=0)


def test_sampling_does_not_disturb_the_global_numpy_random_state():
    np.random.seed(7)
    before = np.random.random()
    np.random.seed(7)
    dissCompare(DISS, GROUP, s=5, seed=4242, weights=WEIGHTS)
    after = np.random.random()

    assert before == after


def test_s0_ignores_the_seed_entirely():
    """s=0 takes no draws, so the seed cannot matter."""
    first = dissCompare(DISS, GROUP, s=0, seed=1, weights=WEIGHTS)
    second = dissCompare(DISS, GROUP, s=0, seed=999_999, weights=WEIGHTS)

    pd.testing.assert_frame_equal(first, second, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL)
