"""The test run shows the loaded ICU, and the documents give the version-proof rule (TX-AUD-20261019-053).

Owner decision of 2026-10-22 (option E). The automatic state order follows the ICU library that is loaded, as R's
follows its own, and ICU versions change over time on both sides; a few characters (``₨``, ``₧``) order differently
between versions. An explicit ``alphabet=`` is kept by R and traminex under any ICU version, so the documents give
that rule; the test run prints the loaded ICU version, and nothing fails or warns when it differs from R's.
"""

from __future__ import annotations

import ctypes.util
from pathlib import Path

import conftest

import traminex as tx
from traminex._categorical_order import _icu_version

REPO = Path(__file__).resolve().parents[1]


def _header() -> str:
    lines = conftest.pytest_report_header(None)
    return " ".join(lines) if isinstance(lines, list) else str(lines)


def test_the_test_run_names_the_loaded_icu() -> None:
    library = ctypes.util.find_library("icui18n")
    text = _header()
    assert library is not None and library in text
    assert f"major version {_icu_version(library)}" in text


def test_the_header_does_not_stop_the_run_when_icu_is_missing_or_unreadable(monkeypatch) -> None:
    monkeypatch.setattr(ctypes.util, "find_library", lambda name: None)
    assert "not found" in _header()
    monkeypatch.setattr(ctypes.util, "find_library", lambda name: "libsomething.so")
    assert "version unreadable" in _header()


def test_the_documents_give_the_alphabet_rule() -> None:
    def flat(text: str) -> str:
        return " ".join(text.split())

    readme = flat((REPO / "README.md").read_text(encoding="utf-8"))
    dev_build = flat((REPO / "docs" / "DEV_BUILD.md").read_text(encoding="utf-8"))
    seqdef_doc = flat(tx.seqdef.__doc__ or "")
    for text in (readme, dev_build):
        assert "pass `alphabet=`" in text
        assert "need the same ICU major version as R" not in text
        assert "needs R's ICU major version" not in text
    assert "pass ``alphabet=``" in seqdef_doc


def test_the_header_survives_a_failing_library_lookup(monkeypatch) -> None:
    # Follow-up review m-10: the lookup itself sits inside the guard, as the docstring promises.
    def fail(name: str) -> str:
        raise OSError("simulated lookup failure")

    monkeypatch.setattr(ctypes.util, "find_library", fail)
    assert "lookup failed" in _header()


def test_an_unreadable_version_is_not_described_as_the_order_in_use(monkeypatch) -> None:
    # Follow-up review m-9: traminex refuses to load a library whose version it cannot read.
    monkeypatch.setattr(ctypes.util, "find_library", lambda name: "libsomething.so")
    text = _header()
    assert "version unreadable" in text and "the automatic state order follows it" not in text


def test_the_alphabet_rule_is_scoped_to_the_state_order() -> None:
    # Follow-up review m-4: alphabet= fixes the state order only; other text orders follow the loaded ICU too.
    for path in (REPO / "README.md", REPO / "docs" / "DEV_BUILD.md"):
        text = " ".join(path.read_text(encoding="utf-8").split())
        assert "other text orders" in text and "follow the loaded ICU as well" in text, path


def test_a_quiet_run_still_names_the_loaded_icu_once() -> None:
    # 2026-10-25 study, finding F-03: make test-python and make test-pybind11 pass -q, which hides the header.
    import subprocess
    import sys

    target = "tests/test_icu_version_report.py::test_the_documents_give_the_alphabet_rule"
    for flags in (["-q"], []):
        command = [sys.executable, "-m", "pytest", *flags, "-p", "no:cacheprovider", target]
        run = subprocess.run(
            command, cwd=REPO, capture_output=True, text=True, timeout=300, check=False
        )
        assert run.returncode == 0, run.stdout + run.stderr
        icu_lines = [line for line in run.stdout.splitlines() if line.startswith("ICU: ")]
        assert len(icu_lines) == 1, flags
