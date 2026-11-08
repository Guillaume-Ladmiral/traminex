"""release-check names its prerequisites (TX-AUD-20261010-046).

In an environment without the development install, ``make release-check``
stopped with a bare ``ModuleNotFoundError`` (``tokenize_rt`` in format-check,
``traminex`` in docs-check). The install that provides both is documented
(docs/DEV_BUILD.md: the editable install with the ``dev`` extra, whose
``black[jupyter]`` brings ``tokenize-rt``), but the gate did not say so.
``scripts/check_dev_environment.py`` now runs first and names every missing
module and that install command.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import check_dev_environment  # noqa: E402

INSTALL = "-m pip install -e '.[plot,perf,extras,tree,dev]'"


def test_missing_modules_are_named_with_the_install_command(monkeypatch, capsys) -> None:
    real = check_dev_environment.find_spec

    def find_spec(name, *args, **kwargs):
        return None if name in {"tokenize_rt", "traminex"} else real(name, *args, **kwargs)

    monkeypatch.setattr(check_dev_environment, "find_spec", find_spec)
    assert check_dev_environment.main() == 1
    out = capsys.readouterr().out
    assert "tokenize_rt" in out and "traminex" in out
    assert INSTALL in out


def test_a_provisioned_environment_passes(capsys) -> None:
    assert check_dev_environment.main() == 0
    assert "development environment OK" in capsys.readouterr().out


def test_release_check_runs_it_first() -> None:
    makefile = (REPO / "Makefile").read_text(encoding="utf-8")
    line = next(line for line in makefile.splitlines() if line.startswith("release-check:"))
    assert line.split(":", 1)[1].split()[0] == "dev-environment-check"
    assert "import IPython; import tokenize_rt" not in makefile
