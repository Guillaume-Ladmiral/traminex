"""The Numba backend must agree with the reference kernels at every dispatch site.

Numba is a **third** distance backend, not a variant of the two ``TRAMINEX_BINDING``
selects: it is active by default whenever ``numba`` imports, it was active for every
suite run and in CI, and it changes which kernel computes a result.  Against that,
the entire cross-backend agreement suite used to be one LCS case in
``tests/test_benchmarks_cli.py`` -- 1 of the 7 production dispatch sites -- while
the compiled backend had 211 byte-identity cases (cycle-4b F-008).

The seven sites in ``traminex/distances/python_kernels.py`` are OM (scalar and
matrix ``subm``), LCS, LCP, RLCP, ``DUR`` with ``metric="chi2"``, ``JSD_TRANS``,
and DHD *inside the ``_is_legacy_dhd`` branch* -- the default TraMineR-faithful
DHD is not Numba-routed.  All eight configurations below are covered.  Since
2026-11-07 two more sites run the NMSMST and SVRspell counts compiled (2026-10-25
study, finding F-06); their routes are byte-identical, so the sites number nine.

Both sides run in subprocesses selected by the **supported** ``TRAMINEX_NUMBA``
switch rather than by monkeypatching the private ``_HAVE_NUMBA`` flag, so the test
exercises the same mechanism a user has.

Tolerances are measured, not assumed. On the 4-state fixture seven routes are
byte-identical (``np.array_equal``, atol = 0) and only ``JSD_TRANS`` differs, at
1.11e-16, from log/division reassociation inside the parallel kernel. ``DUR``
chi2 is byte-identical there only because the fixture has 4 states: on 9 states
it differs by 2.22e-16 (the ``DUR_chi2_nine_states`` route; up to 4.44e-16 over
a 306-case battery, 2026-09-22 audit S1b-14). Both are given 1e-15 -- well below
the smallest distance the fixtures produce and just above the measured
difference, so a real divergence still fails.

The difference grows with the input, so ``docs/DISTANCES.md`` states no fixed
figure (TX-AUD-20261010-068). The ``_large`` routes use 20 states and 60
sequences of length 60: ``JSD_TRANS`` differs there by 8.88e-16 (4 rounding
units; up to 9.99e-16 over five seeds) and ``DUR`` chi2 by 2.22e-16, measured
2026-10-22. They are given 2e-15 and 1e-15.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]

pytest.importorskip("numba", reason="the Numba backend is not installed")

#: route -> (kwargs, absolute tolerance). atol 0 means byte-identical.
ROUTES: dict[str, tuple[dict, float]] = {
    "OM_scalar": ({"method": "OM", "indel": 1.0, "subm": 2.0}, 0.0),
    "OM_matrix": ({"method": "OM", "indel": 1.0, "subm": "MATRIX"}, 0.0),
    "LCS": ({"method": "LCS"}, 0.0),
    "LCP": ({"method": "LCP"}, 0.0),
    "RLCP": ({"method": "RLCP"}, 0.0),
    "DUR_chi2": ({"method": "DUR", "metric": "chi2"}, 0.0),
    "DUR_chi2_nine_states": ({"method": "DUR", "metric": "chi2", "_states": 9}, 1e-15),
    "JSD_TRANS": ({"method": "JSD_TRANS"}, 1e-15),
    "JSD_TRANS_large": ({"method": "JSD_TRANS", "_states": 20}, 2e-15),
    "DUR_chi2_large": ({"method": "DUR", "metric": "chi2", "_states": 20}, 1e-15),
    "DHD_legacy": ({"method": "DHD", "variant": "mean"}, 0.0),
    "NMSMST": ({"method": "NMSMST"}, 0.0),
    "NMSMST_tpow": ({"method": "NMSMST", "tpow": 0.5}, 0.0),
    "SVRspell": ({"method": "SVRspell"}, 0.0),
    "SVRspell_tpow": ({"method": "SVRspell", "tpow": 0.5}, 0.0),
}

_COMPUTE = """
import json, sys
import numpy as np
sys.path.insert(0, {ROOT!r})
from traminex.distances import pairwise
import traminex.distances.python_kernels as kernels

rng = np.random.default_rng(20260823)
corpus = rng.integers(0, 4, size=(24, 10)).astype(np.int32)
corpus_nine = np.random.default_rng(20260823).integers(0, 9, size=(24, 10)).astype(np.int32)
corpus_large = np.random.default_rng(20260823).integers(0, 20, size=(60, 60)).astype(np.int32)
matrix_subm = np.full((4, 4), 2.0)
np.fill_diagonal(matrix_subm, 0.0)
matrix_subm[0, 1] = matrix_subm[1, 0] = 1.5

ROUTES = {ROUTES!r}
result = {{"dispatch_active": bool(kernels._HAVE_NUMBA), "routes": {{}}}}
for name, options in ROUTES.items():
    options = dict(options)
    if options.get("subm") == "MATRIX":
        options["subm"] = matrix_subm
    source = {{4: corpus, 9: corpus_nine, 20: corpus_large}}[options.pop("_states", 4)]
    values = np.asarray(pairwise(source.copy(), **options), dtype=float)
    result["routes"][name] = [float(v).hex() for v in values.ravel()]
    result.setdefault("shapes", {{}})[name] = list(values.shape)
print(json.dumps(result))
"""


def _run(numba_mode: str) -> dict:
    script = _COMPUTE.format(ROOT=str(ROOT), ROUTES={k: v[0] for k, v in ROUTES.items()})
    completed = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
        env={
            **{k: v for k, v in __import__("os").environ.items() if k != "TRAMINEX_NUMBA"},
            "TRAMINEX_NUMBA": numba_mode,
            "TRAMINEX_BINDING": "python",
        },
    )
    return json.loads(completed.stdout.strip().splitlines()[-1])


@pytest.fixture(scope="module")
def backends() -> tuple[dict, dict]:
    accelerated = _run("auto")
    reference = _run("off")
    assert accelerated["dispatch_active"] is True, "TRAMINEX_NUMBA=auto did not enable dispatch"
    assert reference["dispatch_active"] is False, "TRAMINEX_NUMBA=off did not disable dispatch"
    return accelerated, reference


def _values(payload: dict, route: str) -> np.ndarray:
    flat = np.array([float.fromhex(v) for v in payload["routes"][route]], dtype=float)
    return flat.reshape(payload["shapes"][route])


@pytest.mark.parametrize("route", sorted(ROUTES), ids=sorted(ROUTES))
def test_numba_agrees_with_the_reference_kernel(route, backends):
    accelerated, reference = backends
    got = _values(accelerated, route)
    expected = _values(reference, route)
    tolerance = ROUTES[route][1]

    assert got.shape == expected.shape
    if tolerance == 0.0:
        assert np.array_equal(got, expected), (
            f"{route}: expected byte-identical results, "
            f"max|delta| = {float(np.max(np.abs(got - expected)))!r}"
        )
    else:
        np.testing.assert_allclose(got, expected, rtol=0.0, atol=tolerance)


def test_every_numba_dispatch_site_is_covered():
    """Guard the denominator: a new ``if _HAVE_NUMBA:`` site must gain a route here."""
    source = (ROOT / "traminex" / "distances" / "python_kernels.py").read_text(encoding="utf-8")
    # One definition guard (``if _NUMBA_IMPORTED:``) plus the production dispatch
    # sites. Count only the dispatch conditions.
    dispatch_sites = sum(
        1
        for line in source.splitlines()
        if "_HAVE_NUMBA" in line and line.strip().startswith(("if ", "elif "))
    )
    assert dispatch_sites == 9, (
        f"python_kernels.py has {dispatch_sites} Numba dispatch sites but this module "
        f"covers 9; add the new route to ROUTES."
    )


def test_the_supported_switch_rejects_an_unknown_value():
    completed = subprocess.run(
        [sys.executable, "-c", "import traminex.distances.python_kernels"],
        capture_output=True,
        text=True,
        cwd=ROOT,
        env={
            **{k: v for k, v in __import__("os").environ.items()},
            "TRAMINEX_NUMBA": "sometimes",
        },
    )
    assert completed.returncode != 0
    assert "TRAMINEX_NUMBA must be one of: auto, off" in completed.stderr


def test_the_switch_is_read_once_at_import():
    """Documented behaviour: both backend switches are read when traminex is imported.

    Assigning to ``os.environ`` afterwards has no effect, which is why
    ``docs/USAGE_new_apis.md`` tells callers to set them before the import. Pinned
    so the semantics cannot silently become a per-call lookup.
    """
    import os as _os

    script = (
        f"import sys, os; sys.path.insert(0, {str(ROOT)!r});\n"
        "import traminex.distances.python_kernels as k;\n"
        "before = k._HAVE_NUMBA;\n"
        "os.environ['TRAMINEX_NUMBA'] = 'off';\n"
        "print(before, k._HAVE_NUMBA)\n"
    )
    completed = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
        env={k: v for k, v in _os.environ.items() if k != "TRAMINEX_NUMBA"},
    )
    before, after = completed.stdout.split()
    assert before == "True", completed.stdout
    assert after == "True", "TRAMINEX_NUMBA must not take effect after import"


def test_distances_md_states_no_fixed_figure_for_the_differences():
    text = " ".join((ROOT / "docs" / "DISTANCES.md").read_text(encoding="utf-8").split())
    assert "measured at most 2.2e-16 and 4.4e-16" not in text
    assert "grows with the number of states and the sequence length" in text


def test_the_compiled_weighted_pair_sum_is_the_numpy_one_bit_for_bit(monkeypatch):
    """TX-AUD-20261010-083: dissassoc's permutation test spent most of its time here."""
    import traminex.diss as diss
    import traminex.distances.python_kernels as kernels

    rng = np.random.default_rng(20261107)
    for trial in range(300):
        n = int(rng.integers(0, 40))
        D = rng.uniform(0.0, 10.0, size=(n, n))
        D = D + D.T
        if trial % 3 == 0:
            D = np.asfortranarray(D.T)
        w = rng.uniform(0.0, 3.0, size=n)
        if trial % 5 == 0 and n:
            w[int(rng.integers(0, n))] = rng.choice([0.0, -1.0, np.inf, np.nan])
        size = int(rng.integers(0, n + 1))
        members = rng.choice(n, size=size, replace=trial % 4 == 0) if n else np.zeros(0, int)
        with np.errstate(all="ignore"):
            monkeypatch.setattr(kernels, "_HAVE_NUMBA", True)
            compiled = diss._c_weighted_pair_sum(D, members, w)
            monkeypatch.setattr(kernels, "_HAVE_NUMBA", False)
            reference = diss._c_weighted_pair_sum(D, members, w)
        assert np.array(compiled).tobytes() == np.array(reference).tobytes(), trial


def test_dissassoc_is_the_same_with_and_without_numba(monkeypatch):
    import traminex as tx
    import traminex.distances.python_kernels as kernels

    rng = np.random.default_rng(20261108)
    points = rng.uniform(0.0, 5.0, size=(60, 2))
    D = np.sqrt(((points[:, None, :] - points[None, :, :]) ** 2).sum(axis=2))
    group = rng.choice(["a", "b", "c"], size=60)
    weights = rng.uniform(0.5, 2.0, size=60)
    results = []
    for mode in (True, False):
        monkeypatch.setattr(kernels, "_HAVE_NUMBA", mode)
        results.append(
            tx.dissassoc(
                D, group=group, weights=weights, R=50, weight_permutation="diss", random_state=3
            )
        )
    compiled, reference = results
    for key in ("t0", "t", "pval"):
        assert (
            np.asarray(compiled.perms[key], float).tobytes()
            == np.asarray(reference.perms[key], float).tobytes()
        )
