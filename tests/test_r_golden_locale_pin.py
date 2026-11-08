"""Every R-golden generator pins its locale (TX-AUD-20261010-019).

Run under ``LC_ALL=C``, the generators changed 24 tracked fixtures: R sorted
text by bytes, read the scripts' non-ASCII literals as bytes (so ``é`` sorted
first even where ``icuSetCollate(locale = "root")`` was set), and quoted
messages with ``'`` instead of ``‘ ’``. Under ``C.UTF-8`` every fixture was
byte-identical to the tracked one. Each script now starts by switching to
``C.UTF-8`` with the ICU root collation (stopping if the locale is missing),
and the provenance check runs the scripts with ``LC_ALL=C.UTF-8`` and records
it. The locale does not set R's message language: from a shell with
``LANGUAGE=de`` the scripts still wrote 45 fixtures in German (batch-5 review
I-1), so the pin also sets the language to English, and it requires an R built
with ICU, without which ``icuSetCollate`` only warns.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import check_r_golden_provenance  # noqa: E402

PIN = (
    'stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); '
    'invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")'
)
SCRIPTS = sorted((REPO / "scripts" / "r_goldens").glob("*.R"))


def _first_code_line(path: Path) -> str:
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            return stripped
    return ""


@pytest.mark.parametrize("script", SCRIPTS, ids=[path.name for path in SCRIPTS])
def test_every_script_pins_its_locale(script: Path) -> None:
    assert _first_code_line(script) == PIN


#: The only later switch allowed: byte collation, the same from any shell
#: (make_goldens_audit20260922_s6.R records C-collation orders on purpose).
BYTE_COLLATION = 'invisible(Sys.setlocale("LC_COLLATE", "C"))'


@pytest.mark.parametrize("script", SCRIPTS, ids=[path.name for path in SCRIPTS])
def test_no_script_switches_back_to_a_session_collation(script: Path) -> None:
    # A later Sys.setlocale to a UTF-8 collation resets R's ICU collator: it undid
    # the pin in the three scripts that switched LC_COLLATE themselves when R was
    # started from a C shell.
    code = [
        line.strip()
        for line in script.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]
    later = [line for line in code[1:] if "Sys.setlocale" in line]
    assert all(line == BYTE_COLLATION for line in later), later


def test_the_execute_environment_pins_the_locale(monkeypatch) -> None:
    monkeypatch.setenv("LC_ALL", "C")
    monkeypatch.setenv("LANGUAGE", "fr")
    env = check_r_golden_provenance._child_environment(None)
    assert env["LC_ALL"] == "C.UTF-8"
    assert env["LANGUAGE"] == "en"
