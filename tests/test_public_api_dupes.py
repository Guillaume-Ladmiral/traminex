import ast
import importlib
import inspect
import subprocess
import sys
from pathlib import Path

import pytest

from traminex.api_manifest import PUBLIC_API

R3_PUBLIC_ADDITIONS = {
    "createdatadiscrete",
    "group_p",
    "seqerulesdisc",
    "seqgen_missing",
    "seqsamm",
    "seqsammeha",
    "seqsammseq",
    "sorti",
    "sortv",
}
R4_PUBLIC_ADDITIONS = {
    "plot_SAMM",
    "plot_seqimplic",
    "seqplot_rf",
    "seqplot_tentrop",
    "seqplot_tentrop_m",
    "seqsplot",
}


def test_public_api_no_duplicates():
    mod = importlib.import_module("traminex")
    names = list(getattr(mod, "__all__", []))
    assert len(names) == len(set(names)), "__all__ contains duplicate entries"


def test_public_api_matches_manifest():
    mod = importlib.import_module("traminex")
    names = list(getattr(mod, "__all__", []))
    assert names == list(PUBLIC_API), "Runtime __all__ must match traminex.api_manifest.PUBLIC_API"


def test_public_tree_visualization_callables_have_runtime_docstrings():
    mod = importlib.import_module("traminex")
    names = {
        "disstree2dotp",
        "disstree_plot",
        "disstreedisplay",
        "seqtree2dot",
        "seqtreedisplay",
    }

    assert names <= set(PUBLIC_API)
    missing = sorted(name for name in names if not inspect.getdoc(getattr(mod, name)))
    assert missing == []


def test_seqeapplysub_top_level_export_identity_and_signature():
    mod = importlib.import_module("traminex")
    events = importlib.import_module("traminex.events")

    assert list(PUBLIC_API).count("seqeapplysub") == 1
    assert list(mod.__all__).count("seqeapplysub") == 1
    assert mod.seqeapplysub is events.seqeapplysub
    assert str(inspect.signature(mod.seqeapplysub)) == (
        "(subseq: 'SubSeqElist', method: 'Optional[Union[int, str]]' = None, "
        "constraint: 'Optional[SeqEConstraint]' = None, rules: 'bool' = False) "
        "-> 'pd.DataFrame'"
    )


def test_r3_closed_on_exact_220_name_surface_and_seqsep_stays_stable():
    mod = importlib.import_module("traminex")
    stats = importlib.import_module("traminex.stats")

    public_names = set(PUBLIC_API)
    landed_r4 = public_names & R4_PUBLIC_ADDITIONS
    assert len(PUBLIC_API) == len(public_names) == 220 + len(landed_r4)
    assert len(public_names - R3_PUBLIC_ADDITIONS - landed_r4) == 211
    assert R3_PUBLIC_ADDITIONS <= public_names
    assert list(PUBLIC_API).count("seqsep") == 1
    assert list(mod.__all__).count("seqsep") == 1
    assert mod.seqsep is stats.seqsep
    assert str(inspect.signature(mod.seqsep)) == "(seqdata, sl: 'int' = 1, sep: 'str' = '-')"


def test_public_api_lazy_resolution_smoke():
    import sys

    sys.modules.pop("traminex", None)
    mod = importlib.import_module("traminex")
    assert "seqdef" not in mod.__dict__
    _ = mod.seqdef
    assert "seqdef" in mod.__dict__


def test_lazy_resolution_cannot_shadow_function_named_like_submodule():
    code = """
import traminex as tx
_ = tx.seqdef
assert callable(tx.seqformat), repr(tx.seqformat)
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("name", ["disstree", "seqdiff", "seqemlt", "seqformat", "seqtree"])
def test_submodule_first_import_cannot_shadow_public_callable(name):
    code = f"""
import importlib
importlib.import_module('traminex.{name}')
import traminex
assert callable(getattr(traminex, '{name}')), repr(getattr(traminex, '{name}'))
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr


def test_stub_public_api_matches_manifest():
    stub_path = Path(__file__).resolve().parents[1] / "stubs" / "traminex" / "__init__.pyi"
    src = stub_path.read_text(encoding="utf-8")
    tree = ast.parse(src, filename=str(stub_path))
    stub_all = None
    for node in tree.body:
        if isinstance(node, ast.Assign):
            if any(isinstance(t, ast.Name) and t.id == "__all__" for t in node.targets):
                if isinstance(node.value, (ast.List, ast.Tuple)):
                    vals = []
                    ok = True
                    for elt in node.value.elts:
                        if isinstance(elt, ast.Constant) and isinstance(elt.value, str):
                            vals.append(elt.value)
                        else:
                            ok = False
                            break
                    if ok:
                        stub_all = vals
                        break
    assert stub_all is not None, "Could not parse __all__ from stubs/traminex/__init__.pyi"
    assert stub_all == list(PUBLIC_API), "Stub __all__ must match traminex.api_manifest.PUBLIC_API"


def test_public_api_excludes_removed_legacy_names():
    mod = importlib.import_module("traminex")
    names = set(getattr(mod, "__all__", []))
    # RP7a (Visualization) ships the missing base-R plots under their R names
    # (seqfplot now; seqmtplot/seqmsplot/seqHtplot/seqdHplot/seqpcplot/seqplotMD/
    # seqplot dispatcher as they land) -- each is dropped from this guard in the
    # commit that adds it. Kept here: genuinely-removed legacy stubs + R aliases
    # traminex does NOT expose (it keeps descriptive names sequence_index_plot/
    # state_distribution_plot/state_legend; seqgbar is an internal helper).
    removed = {
        "states_to_events",
        "seqiplot",
        "seqIplot",
        "seqdplot",
        "seqgbar",
        "seqlegend",
        "OMloc",
        "seqOMloc",
        "representative_sequence_index",
    }
    present = sorted(removed.intersection(names))
    assert not present, f"Legacy names unexpectedly present in __all__: {present}"
