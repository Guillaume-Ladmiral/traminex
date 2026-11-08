"""The stub gate checks that a stub re-export names the runtime object (TX-AUD-20261010-009).

``from .stats import seqtransn as seqtab`` in the root stub left every name and
signature list intact, so the gate passed while the stub described a different
function under a public name. The checker now resolves each re-export to its
source module and original name and requires the runtime attribute to be that
object.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import check_stub_signatures  # noqa: E402


def _copy_with(tmp_path: Path, old: str, new: str) -> Path:
    stub_root = tmp_path / "traminex"
    shutil.copytree(REPO / "stubs" / "traminex", stub_root)
    init = stub_root / "__init__.pyi"
    text = init.read_text(encoding="utf-8")
    assert text.count(old) == 1
    init.write_text(text.replace(old, new), encoding="utf-8")
    return stub_root


@pytest.mark.parametrize(
    "old, new, wrong",
    [
        ("    seqtab as seqtab,\n", "    seqtransn as seqtab,\n", "traminex.stats.seqtransn"),
        (
            "    seqrep as seqrep,\n",
            "    _label_value as seqrep,\n",
            "traminex.cluster._label_value",
        ),
    ],
)
def test_a_root_reexport_of_another_object_fails(tmp_path, capsys, old, new, wrong) -> None:
    stub_root = _copy_with(tmp_path, old, new)
    status = check_stub_signatures.main(["--stub-root", str(stub_root)])
    output = capsys.readouterr().out
    assert status == 1
    exported = new.split(" as ")[1].rstrip(",\n")
    assert f"traminex: '{exported}' re-exports {wrong}" in output


def test_the_tracked_stubs_pass(capsys) -> None:
    assert check_stub_signatures.main([]) == 0
    assert "Re-exports checked:" in capsys.readouterr().out


def test_a_renamed_reexport_in_a_stub_without_a_literal_all_fails(tmp_path, capsys) -> None:
    """``x as y`` is checked when ``y`` is public at run time, with or without ``__all__``.

    Only ``name as name`` was read in a stub without a literal ``__all__``, so
    ``indel_from_subm as seqcost`` in ``distances/__init__.pyi`` (``__all__: tuple[str, ...]``)
    passed though ``traminex.distances.seqcost`` is another function (TX-REP-20261019-004,
    batch-5 review M-1).
    """
    stub_root = tmp_path / "traminex"
    shutil.copytree(REPO / "stubs" / "traminex", stub_root)
    stub = stub_root / "distances" / "__init__.pyi"
    text = stub.read_text(encoding="utf-8")
    assert text.count("    seqcost as seqcost,\n") == 1
    stub.write_text(
        text.replace("    seqcost as seqcost,\n", "    indel_from_subm as seqcost,\n"),
        encoding="utf-8",
    )

    status = check_stub_signatures.main(["--stub-root", str(stub_root)])

    assert status == 1
    assert (
        "traminex.distances: 'seqcost' re-exports traminex.distances.costs.indel_from_subm"
        in capsys.readouterr().out
    )


def test_a_renamed_reexport_of_a_name_the_runtime_lacks_fails(tmp_path, capsys) -> None:
    """``x as y`` with no runtime ``y`` was skipped in a stub without a literal ``__all__``
    (batch-16 review M-4); the ``y as y`` form of a missing name was already reported."""
    stub_root = tmp_path / "traminex"
    shutil.copytree(REPO / "stubs" / "traminex", stub_root)
    stub = stub_root / "distances" / "__init__.pyi"
    text = stub.read_text(encoding="utf-8")
    assert text.count("from .omloc import omloc as omloc\n") == 1
    stub.write_text(
        text.replace(
            "from .omloc import omloc as omloc\n",
            "from .omloc import omloc as omloc\nfrom .costs import seqcost as seqcost_alias\n",
        ),
        encoding="utf-8",
    )

    status = check_stub_signatures.main(["--stub-root", str(stub_root)])

    assert status == 1
    assert (
        "traminex.distances: 'seqcost_alias' re-exports traminex.distances.costs.seqcost, "
        "but the runtime module has no 'seqcost_alias'" in capsys.readouterr().out
    )
