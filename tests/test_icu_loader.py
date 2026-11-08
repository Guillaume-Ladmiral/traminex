"""ICU documentation and the portable version parser (TX-AUD-20261010-030).

``seqtab`` orders its labels with ICU whether or not frequencies are tied (with
ICU absent, counts 3, 2, 1 raised ``RuntimeError``), while docs/DEV_BUILD.md
listed it under "tied frequencies". The loader read the ICU version only from a
Linux ``.so.N`` name, so a macOS ``libicui18n.78.dylib`` failed with an unexplained
message; an absolute path is now also read through its symlinks. (A Windows
``icuinN.dll`` form added then could not be reached and was removed,
TX-REP-20261019-005.) And the order follows the
loaded ICU: conda's ICU 78 orders ``₧ ₨``, ICU 74 (and R with ICU 74.2) ``₨ ₧``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx
import traminex._categorical_order as order

REPO = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    "name, version",
    [
        ("libicui18n.so.74", "74"),
        ("/opt/conda/lib/libicui18n.so.78.3", "78"),
        ("libicui18n.78.dylib", "78"),
        ("/usr/local/opt/icu4c/lib/libicui18n.78.1.dylib", "78"),
    ],
)
def test_the_icu_version_is_read_from_each_library_name_form(name: str, version: str) -> None:
    assert order._icu_version(name) == version


def test_an_unknown_library_name_is_refused_by_name() -> None:
    with pytest.raises(RuntimeError, match="libicui18n.dylib"):
        order._icu_version("libicui18n.dylib")


# TX-REP-20261019-005 (batch-5 review M-4): find_library("icui18n") gives, on macOS, the
# path dyld finds, typically the unversioned libicui18n.dylib symlink, and on Windows only
# an icui18n[.dll] name; the icuinN.dll form could not be reached.


@pytest.mark.parametrize(
    "link, target, version",
    [
        ("libicui18n.dylib", "libicui18n.74.2.dylib", "74"),
        ("libicui18n.so", "libicui18n.so.74.2", "74"),
    ],
)
def test_an_absolute_library_path_is_read_through_its_symlink(
    tmp_path: Path, link: str, target: str, version: str
) -> None:
    (tmp_path / target).write_bytes(b"")
    (tmp_path / link).symlink_to(target)

    assert order._icu_version(str(tmp_path / link)) == version


@pytest.mark.parametrize("name", ["icuin78.dll", "C:\\icu\\bin\\icui18n.dll"])
def test_a_windows_library_name_is_refused_with_the_supported_form(name: str) -> None:
    with pytest.raises(RuntimeError, match=r"libicui18n\.so\.N"):
        order._icu_version(name)


def test_seqtab_needs_icu_without_ties(monkeypatch) -> None:
    def absent():
        raise RuntimeError("simulated: ICU could not be loaded")

    monkeypatch.setattr(order, "_icu_root_collation_api", absent)
    seqs = tx.seqdef(
        np.array([list("AAB")] * 3 + [list("ABB")] * 2 + [list("BBB")], dtype=object),
        alphabet=["A", "B"],
    )
    with pytest.raises(RuntimeError, match="simulated"):
        tx.seqtab(seqs)


def test_the_documents_describe_seqtab_and_the_icu_version() -> None:
    dev_build = (REPO / "docs" / "DEV_BUILD.md").read_text(encoding="utf-8")
    readme = (REPO / "README.md").read_text(encoding="utf-8")
    assert "seqtab` (tied frequencies)" not in dev_build
    # Owner decision of 2026-10-22 (TX-AUD-20261019-053, option E): the documents say the order follows the loaded
    # ICU, whose versions change over time, and give the version-proof rule (an explicit alphabet=).
    for text in (dev_build, readme):
        flat = " ".join(text.split())
        assert "ICU versions change over time" in flat and "pass `alphabet=`" in flat
