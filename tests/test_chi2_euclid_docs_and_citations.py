"""CHI2/EUCLID documentation and stale citations (TX-AUD-20261019-027, TX-AUD-20261011-053).

``docs/NEW_DISTANCE_METHODS.md`` listed ``weights`` for EUCLID, which R ignores (``seqdist-CHI2.R:6-8`` sets
``weighted <- FALSE``), called the CHI2/EUCLID ``norm`` a boolean, said ``norm`` only divides by the square root of
the number of intervals, and used the deprecated ``norm=True`` in its example. Several places cited "seqdist.R line
205" (a comment) for the check at lines 206-208. ``ROADMAP.md`` gave the 2026-08-13 audit as the current one, and
``docs/INDEX.md``, which calls itself the authoritative map, left out three tracked directories and two pages.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pytest

import traminex as tx

REPO = Path(__file__).resolve().parents[1]
SEQMAT = np.array(
    [[0, 0, 1, 1, 2, 2], [0, 1, 1, 2, 2, 2], [2, 2, 1, 0, 0, 1], [0, 0, 0, 1, 1, 0]], dtype=np.int32
)


def _documented_chi2(seqmat: np.ndarray, step: int, norm: bool, euclid: bool) -> np.ndarray:
    """The distance as the document describes it, for a ``step`` dividing the length."""
    n, length = seqmat.shape
    states = int(seqmat.max()) + 1
    intervals = [list(range(start, start + step)) for start in range(0, length, step)]
    total = np.zeros((n, n))
    for columns in intervals:
        counts = np.array(
            [[np.sum(row[columns] == k) for k in range(states)] for row in seqmat], float
        )
        shares = counts / counts.sum(axis=1, keepdims=True)
        overall = counts.sum(axis=0) / counts.sum()
        present = overall > 0
        if euclid:
            scale = np.full(states, 2.0 if norm else 1.0)
        else:
            scale = overall.copy()
            if norm and present.sum() > 1:
                smallest = np.sort(overall[present])
                scale *= 1 / smallest[0] + 1 / smallest[1]
        for i in range(n):
            for j in range(n):
                diff = shares[i, present] - shares[j, present]
                total[i, j] += float(np.sum(diff**2 / scale[present]))
    distances = np.sqrt(total)
    return distances / np.sqrt(len(intervals)) if norm else distances


@pytest.mark.parametrize("method", ["CHI2", "EUCLID"])
@pytest.mark.parametrize("step", [1, 2])
@pytest.mark.parametrize("norm", ["none", "auto"])
def test_the_documented_chi2_and_euclid_formulas_are_the_code_s(
    method: str, step: int, norm: str
) -> None:
    actual = np.asarray(tx.seqdist(SEQMAT, method=method, step=step, norm=norm), dtype=float)
    expected = _documented_chi2(SEQMAT, step, norm == "auto", method == "EUCLID")
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-15)


@pytest.mark.parametrize("method", ["CHI2", "EUCLID"])
def test_the_default_norm_is_none_and_a_logical_norm_is_deprecated(method: str) -> None:
    default = np.asarray(tx.seqdist(SEQMAT, method=method, step=1), dtype=float)
    none = np.asarray(tx.seqdist(SEQMAT, method=method, step=1, norm="none"), dtype=float)
    np.testing.assert_array_equal(default, none)
    auto = np.asarray(tx.seqdist(SEQMAT, method=method, step=1, norm="auto"), dtype=float)
    for logical, same in ((True, auto), (False, none)):
        with pytest.warns(DeprecationWarning, match=f"norm={logical} is deprecated"):
            value = np.asarray(tx.seqdist(SEQMAT, method=method, step=1, norm=logical), dtype=float)
        np.testing.assert_array_equal(value, same)


def test_euclid_ignores_weights_and_chi2_uses_them() -> None:
    weights = np.array([1.0, 5.0, 0.5, 2.0])
    for method, changes in (("EUCLID", False), ("CHI2", True)):
        plain = np.asarray(tx.seqdist(SEQMAT, method=method, step=1), dtype=float)
        weighted = np.asarray(
            tx.seqdist(SEQMAT, method=method, step=1, weights=weights), dtype=float
        )
        assert (not np.array_equal(plain, weighted)) is changes, method


def test_new_distance_methods_md_states_the_options() -> None:
    text = (REPO / "docs" / "NEW_DISTANCE_METHODS.md").read_text(encoding="utf-8")
    euclid_row = next(line for line in text.splitlines() if line.startswith("| `EUCLID` |"))
    assert "`weights` is accepted and ignored" in euclid_row
    assert "seqdist-CHI2.R:6-8" in text
    assert "`norm` (bool)" not in text
    assert "norm=True" not in text
    assert 'D_chi2 = tx.seqdist(seqmat, method="CHI2", step=1, norm="auto")' in text
    assert "1/p₁ + 1/p₂" in text
    assert "line 205" not in text and "lines 206-208" in text


def test_the_norm_refusal_cites_the_check_s_lines() -> None:
    assert "line 205" not in (REPO / "docs" / "DISTANCES_NORMALIZATION.md").read_text(
        encoding="utf-8"
    )
    with pytest.raises(ValueError, match=r"seqdist\.R:206-208"):
        tx.seqdist(SEQMAT, method="NMS", norm="auto")


def test_roadmap_and_distances_md_status_wording() -> None:
    roadmap = (REPO / "ROADMAP.md").read_text(encoding="utf-8")
    assert "latest on 2026-10-19" in roadmap
    distances = (REPO / "docs" / "DISTANCES.md").read_text(encoding="utf-8")
    assert "R-faithful default (smoothing=0)" not in distances


def _tracked() -> list[str]:
    try:
        out = subprocess.run(
            ["git", "ls-files"], cwd=REPO, capture_output=True, text=True, check=True
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        pytest.skip("not a git checkout")
    return out.splitlines()


def test_index_names_every_tracked_directory_and_docs_page() -> None:
    index = (REPO / "docs" / "INDEX.md").read_text(encoding="utf-8")
    tracked = _tracked()
    directories = sorted({path.split("/")[0] for path in tracked if "/" in path})
    assert [d for d in directories if f"{d}/" not in index] == []
    pages = [Path(p) for p in tracked if p.startswith("docs/") and p.count("/") == 1]
    pages = [
        p
        for p in pages
        if p.suffix == ".md" or (p.suffix == ".qmd" and not (REPO / p.with_suffix(".md")).exists())
    ]
    assert [p.name for p in pages if p.name not in index] == []
    # The census of 2026-10-26 (17 on 2026-10-22, before .codebase-memory/ and graphify-out/ left git);
    # fewer means git listed too little.
    assert len(directories) >= 15
