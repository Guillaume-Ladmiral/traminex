"""``traminex.__version__`` is the version of the files actually imported (TX-AUD-20261011-032).

It returned the installed distribution's version whenever one was installed, so a copied
or second tree imported ahead of it reported the installed version, not its own. The
installed version now counts only when that distribution owns the imported
``__init__.py``: a wheel install puts it beside its metadata (with or without a RECORD), an
editable install names its tree in ``direct_url.json`` (PEP 610). Otherwise the tree's own
``pyproject.toml`` decides.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from importlib.metadata import PackageNotFoundError, PathDistribution
from pathlib import Path

import traminex

REPO = Path(__file__).resolve().parents[1]


def _init_file(root: Path) -> Path:
    init = root / "traminex" / "__init__.py"
    init.parent.mkdir(parents=True)
    init.write_text("", encoding="utf-8")
    return init.resolve()


def _distribution(site: Path, record: list[str] | None, direct_url: dict | None = None):
    info = site / "traminex-1.2.3.dist-info"
    info.mkdir(parents=True)
    (info / "METADATA").write_text(
        "Metadata-Version: 2.1\nName: traminex\nVersion: 1.2.3\n", encoding="utf-8"
    )
    if record is not None:
        (info / "RECORD").write_text("".join(f"{path},,\n" for path in record), encoding="utf-8")
    if direct_url is not None:
        (info / "direct_url.json").write_text(json.dumps(direct_url), encoding="utf-8")
    return PathDistribution(info)


def test_a_wheel_install_owns_the_init_file_its_record_lists(tmp_path: Path) -> None:
    site = tmp_path / "site-packages"
    installed = _init_file(site)
    distribution = _distribution(site, ["traminex/__init__.py"])

    assert traminex._distribution_owns(distribution, installed)
    assert not traminex._distribution_owns(distribution, _init_file(tmp_path / "copy"))


def test_an_install_without_a_record_owns_the_files_beside_it(tmp_path: Path) -> None:
    # Batch-17 review M-3: an installer that writes no RECORD lost the installed version.
    site = tmp_path / "site-packages"
    installed = _init_file(site)
    distribution = _distribution(site, None)

    assert distribution.files is None
    assert traminex._distribution_owns(distribution, installed)
    assert not traminex._distribution_owns(distribution, _init_file(tmp_path / "copy"))


def test_an_editable_install_owns_the_tree_its_direct_url_names(tmp_path: Path) -> None:
    checkout = _init_file(tmp_path / "checkout")
    direct_url = {"url": (tmp_path / "checkout").as_uri(), "dir_info": {"editable": True}}
    distribution = _distribution(
        tmp_path / "site-packages", ["_editable_skbc_traminex.py"], direct_url
    )

    assert traminex._distribution_owns(distribution, checkout)
    assert not traminex._distribution_owns(distribution, _init_file(tmp_path / "worktree"))


def test_a_non_editable_direct_url_claims_no_tree(tmp_path: Path) -> None:
    checkout = _init_file(tmp_path / "checkout")
    direct_url = {"url": (tmp_path / "checkout").as_uri(), "dir_info": {}}
    distribution = _distribution(tmp_path / "site-packages", [], direct_url)

    assert not traminex._distribution_owns(distribution, checkout)


def test_without_a_distribution_the_tree_pyproject_decides(monkeypatch) -> None:
    def missing(name):
        raise PackageNotFoundError(name)

    monkeypatch.setattr(traminex._metadata, "distribution", missing)

    assert traminex._resolve_version() == "0.2.3"


def test_the_checkout_reports_its_own_version() -> None:
    assert Path(traminex.__file__).resolve() == REPO / "traminex" / "__init__.py"
    assert traminex.__version__ == "0.2.3"


def test_a_copied_tree_reports_its_own_pyproject_version(tmp_path: Path) -> None:
    copy = tmp_path / "copy"
    shutil.copytree(
        REPO / "traminex",
        copy / "traminex",
        ignore=shutil.ignore_patterns("__pycache__", "*.so"),
    )
    pyproject = (REPO / "pyproject.toml").read_text(encoding="utf-8")
    assert pyproject.count('version = "0.2.3"') == 1
    (copy / "pyproject.toml").write_text(
        pyproject.replace('version = "0.2.3"', 'version = "9.9.9"'), encoding="utf-8"
    )
    code = (
        "import sys\n"
        "sys.meta_path[:] = [f for f in sys.meta_path"
        " if 'editable' not in (type(f).__module__ + repr(f)).lower()]\n"
        f"sys.path.insert(0, {str(copy)!r})\n"
        "import importlib.metadata, traminex\n"
        "print(traminex.__file__)\n"
        "print(traminex.__version__)\n"
        "print(importlib.metadata.version('traminex'))\n"
    )

    result = subprocess.run(
        [sys.executable, "-c", code], cwd=tmp_path, capture_output=True, text=True, check=True
    )

    imported, version, installed = result.stdout.split()
    assert imported == str(copy / "traminex" / "__init__.py")
    assert installed == "0.2.3"
    assert version == "9.9.9"
