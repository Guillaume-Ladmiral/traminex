"""The build documents describe the environment the gates use (TX-AUD-20261010-084).

docs/DEV_BUILD.md said the ``datasci`` environment installs traminex editable from this
tree; it has neither traminex nor ``tokenize_rt`` (which ``black[jupyter]`` needs for
``make docs-check``), and no document gave the private environment the gates run in. The
recipe is a virtual environment layered on ``datasci`` with ``--system-site-packages``,
which leaves ``datasci`` unchanged. DOCS_POLICY.md, INDEX.md and CONTRIBUTING.md sent the
checks to the ``datasci`` interpreter.
"""

from __future__ import annotations

from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def _flat(name: str) -> str:
    return " ".join((REPO / name).read_text(encoding="utf-8").split())


def test_dev_build_gives_the_private_environment_recipe() -> None:
    text = _flat("docs/DEV_BUILD.md")

    assert "`datasci` environment installs traminex editable" not in text
    assert "-m venv --system-site-packages" in text
    assert "-m pip install --no-deps --no-build-isolation -e ." in text
    assert "tokenize_rt" in text


def test_no_document_sends_the_checks_to_the_datasci_interpreter() -> None:
    assert "with the `datasci` interpreter" not in _flat("docs/DOCS_POLICY.md")
    assert "run with the `datasci` Python" not in _flat("docs/INDEX.md")
    assert "datasci/bin/python scripts/" not in _flat("docs/INDEX.md")
    assert "datasci` environment for maintained local gates" not in _flat("CONTRIBUTING.md")
