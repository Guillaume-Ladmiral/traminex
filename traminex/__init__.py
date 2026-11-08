"""Top-level package exports for traminex.

Public names are defined in ``traminex.api_manifest.PUBLIC_API`` and resolved lazily
on first access to reduce import-time coupling and optional dependency blast radius.
"""

from __future__ import annotations

import importlib
import importlib.util
import sys
import types
from importlib import metadata as _metadata
from pathlib import Path
from typing import Dict, Iterable, Mapping

from .api_manifest import PUBLIC_API, PUBLIC_API_SET

__all__ = list(PUBLIC_API)


def _distribution_owns(distribution: _metadata.Distribution, init_file: Path) -> bool:
    """Whether ``distribution`` installed the ``__init__.py`` at ``init_file`` (resolved).

    A wheel install puts the package beside its metadata, where ``locate_file``
    finds it, whether or not the installer wrote a RECORD (batch-17 review M-3).
    An editable install puts only its redirect files there and names its source
    tree in ``direct_url.json`` (PEP 610).
    """
    installed = Path(str(distribution.locate_file("traminex/__init__.py")))
    if installed.resolve() == init_file:
        return True
    import json
    from urllib.parse import urlparse
    from urllib.request import url2pathname

    try:
        direct_url = json.loads(distribution.read_text("direct_url.json") or "{}")
    except ValueError:
        return False
    url = urlparse(str(direct_url.get("url", "")))
    if url.scheme != "file" or not direct_url.get("dir_info", {}).get("editable"):
        return False
    tree = Path(url2pathname(url.path))
    return (tree / "traminex" / "__init__.py").resolve() == init_file


def _resolve_version() -> str:
    """Return the version of the imported files, single-sourced from ``pyproject.toml``.

    ``pyproject.toml`` ``[project].version`` is the one declared source. When the
    installed (editable or wheel) distribution owns the imported ``__init__.py``,
    its metadata carries that version and ``importlib.metadata.version`` reads it
    back -- so there is no hand-maintained literal to drift. Any other tree (an
    uninstalled checkout, a copy or a second worktree imported ahead of the
    installed one) reads the ``pyproject.toml`` beside it (TX-AUD-20261011-032),
    and finally a clearly unreleased sentinel.
    """
    init_file = Path(__file__).resolve()
    try:
        owned = _distribution_owns(_metadata.distribution("traminex"), init_file)
    except _metadata.PackageNotFoundError:
        owned = False
    if owned:
        return _metadata.version("traminex")
    import re

    pyproject = init_file.parent.parent / "pyproject.toml"
    try:
        text = pyproject.read_text(encoding="utf-8")
    except OSError:
        return "0+unknown"
    m = re.search(r'(?m)^version\s*=\s*"([^"]+)"\s*$', text)
    return m.group(1) if m else "0+unknown"


__version__ = _resolve_version()


_LAZY_EXPORTS: Dict[str, tuple[str, ...]] = {
    "align_p2": ("SeqAlign", "seqalign"),
    "cluster": (
        "cluster_agglomerative",
        "cluster_representatives",
        "cluster_summary",
        "medoid_index",
        "seqrep",
        "seqrep_grp",
    ),
    "config_dataclasses": (
        "DHDConfig",
        "DURConfig",
        "LCSLCPConfig",
        "OMConfig",
        "TRATEConfig",
        "default_configs",
    ),
    "diss": (
        "dissassoc",
        "dissassoc_effect",
        "dissassoc_test",
        "disscenter",
        "disscenter_double_center",
        "dissdomassoc",
        "dissvar",
        "dissvar_grp",
        "dissvar_summary",
        "double_center_matrix",
        "gower_matrix",
        "weighted_inertia",
        "weighted_pairwise_sum",
    ),
    "discrete_rules": ("createdatadiscrete", "seqerulesdisc"),
    "diss_indic": ("dissindic",),
    "diss_merge": ("dissmergegroups",),
    "diss_mfa": ("dissmfacw",),
    "diss_p1": ("dissrep",),
    "domassoc": ("seqdomassoc",),
    "group_compare": ("dissCompare", "seqCompare", "seqLRT", "seqBIC"),
    "seqemlt": ("seqemlt",),
    "disstree": (
        "disstree",
        "disstree2dot",
        "disstree2dot_explore",
        "disstree2dotp",
        "disstree_assign",
        "disstree_feature_importance",
        "disstree_fit",
        "disstree_get_rules",
        "disstree_leaf_summary",
        "disstree_plot",
        "disstree_predict",
        "disstreedisplay",
        "disstreeleaf",
    ),
    "distances": (
        "from_condensed",
        "get_lcp_length",
        "get_rlcp_length",
        "om_trate_pairwise",
        "omloc",
        "sanitize_distance_matrix",
        "seqcost",
        "seqdist",
        "seqLLCP",
        "seqLLCS",
        "seqRLCP",
        "seqsubm",
        "to_condensed",
    ),
    "events": (
        "EventSequence",
        "is_eseq",
        "is_seqelist",
        "is_stslist",
        "is_subseqelist",
        "seqecmpgroup",
        "seqedist",
        "events_to_state_sequences",
        "events_transition_degrees",
        "events_transition_matrix",
        "frequent_event_pairs",
        "frequent_event_pairs_timegap",
        "frequent_events",
        "frequent_ngrams",
        "seqeconstraint",
        "seqecreate",
        "seqeapplysub",
        "seqefsub",
        "seqelength",
        "seqentrans",
        "seqerules",
        "seqetm",
        "seqe2stm",
        "seqeid",
        "seqeweight",
        "seqecontain",
    ),
    "format": ("FCE_to_TSE", "HSPELL_to_STS", "STS_to_SPELL", "STS_to_SRS", "TSE_to_STS", "seqdef"),
    "multidomain": ("seqMD_gridsearch_weights", "seqMD_optimize_weights"),
    "multidomain_cat": ("seqMD", "seqdistmc"),
    "patterns": ("frequent_subsequences",),
    "plotting": (
        "chronogram",
        "seqalignplot",
        "seqdiffplot",
        "seqplotMD",
        "plot_subseqelist",
        "plot_subseqelistchisq",
        "seqpcplot",
        "seqeordplot",
        "plot_eseq",
        "plot_seqelist",
        "seqplot",
        "cluster_chronograms",
        "cluster_chronograms_with_representatives",
        "duration_boxplot",
        "overlay_representatives_on_index_plot",
        "seqfplot",
        "seqfrequency_plot",
        "seqmtplot",
        "seqmsplot",
        "seqHtplot",
        "seqplot_tentrop",
        "seqplot_tentrop_m",
        "seqdHplot",
        "seqrplot",
        "seqrfplot",
        "seqplot_rf",
        "plot_SAMM",
        "plot_seqimplic",
        "seqsplot",
        "seqsurvplot",
        "seqedplot",
        "sequence_index_plot",
        "sequence_index_plot_grouped",
        "state_distribution_plot",
        "state_legend",
        "state_legend_figure",
        "transition_heatmap",
    ),
    "repseq_rf": ("dissrf", "seqrf"),
    "samm": ("seqsamm", "seqsammeha", "seqsammseq"),
    "seqdiff": ("seqdiff",),
    "seqformat": ("STS_to_TSE", "seqdef_from_long", "seqformat", "seqstart", "toPersonPeriod"),
    "seqtree": (
        "seqtree",
        "seqtree2dot",
        "seqtree_cv_fit",
        "seqtree_fit",
        "seqtree_leaf_summary",
        "seqtree_plot",
        "seqtree_predict",
        "seqtreedisplay",
    ),
    "sequences": (
        "SeqMD",
        "STSList",
        "alphabet",
        "cpal",
        "seqgen_missing",
        "seqMD_distances",
        "stlab",
    ),
    "stats": (
        "entropy",
        "group_p",
        "rowmode",
        "seqhasmiss",
        "seqibad",
        "seqindic",
        "seqindic_dyn",
        "seqintegr",
        "seqipos",
        "seqistatd",
        "seqivardur",
        "seqlength",
        "seqlogp",
        "seqmaintokens",
        "seqnum",
        "seqdss",
        "seqdur",
        "seqivolatility",
        "seqprecarity",
        "seqprecstart",
        "seqidegrad",
        "seqinsecurity",
        "seqST",
        "seqasnum",
        "seqstatd",
        "seqtab",
        "seqtabstocc",
        "seqtransn",
        "seqtransn_count",
        "seqtrate",
        "sorti",
        "sortv",
        "turbulence_exact",
        "turbulence_proxy",
    ),
    "stats_p1": (
        "seqfind",
        "seqfpos",
        "seqfposend",
        "seqmeant",
        "seqmodst",
        "seqstatf",
        "seqstatl",
    ),
    "stats_p2": (
        "seqgen",
        "seqici",
        "seqient",
        "seqrecode",
        "seqconc",
        "seqsep",
        "seqdecomp",
        "seqdim",
        "seqmpos",
        "seqcomp",
        "seqsubsn",
        "seqpm",
        "seqimplic",
        "seqsha",
        "seqgranularity",
        "seqpcfilter",
        "seqauto",
        "seqcta",
        "seqsurv",
        "seqpolyads",
    ),
}


def _iter_lazy_items(exports: Mapping[str, Iterable[str]]):
    for module_name, names in exports.items():
        for name in names:
            yield name, module_name


def _build_lazy_index() -> Dict[str, str]:
    out: Dict[str, str] = {}
    for name, module_name in _iter_lazy_items(_LAZY_EXPORTS):
        if name in out:
            raise RuntimeError(f"Duplicate lazy export mapping for '{name}'")
        out[name] = module_name
    missing = sorted(PUBLIC_API_SET - set(out))
    extra = sorted(set(out) - PUBLIC_API_SET)
    if missing or extra:
        raise RuntimeError(
            "Lazy export mapping mismatch with api manifest. " f"Missing={missing}, Extra={extra}"
        )
    return out


_LAZY_ATTR_TO_MODULE = _build_lazy_index()
_MODULE_NAME_EXPORT_COLLISIONS = {
    module_name for module_name, names in _LAZY_EXPORTS.items() if module_name in names
}


def _restore_loaded_submodule_collisions() -> None:
    """Keep public callables from being shadowed by imported submodules."""
    for export_name in _MODULE_NAME_EXPORT_COLLISIONS:
        module = sys.modules.get(f"{__name__}.{export_name}")
        if module is None:
            continue
        value = getattr(module, export_name, None)
        if value is not None:
            globals()[export_name] = value


#: Optional third-party distributions and the extra that installs each, mirroring
#: ``[project.optional-dependencies]`` in ``pyproject.toml``. A missing one used to
#: surface as a bare ``ImportError: No module named 'matplotlib'`` at attribute
#: resolution, which never connected back to the extras ``README.md`` advertises
#: (cycle-4b F-014). The scikit-learn path was already exemplary --
#: ``RuntimeError: scikit-learn is required for disstree_fit`` -- and this matches
#: that wording.
#: Every key must have at least one real import site in the package, or the
#: ``ImportError`` handler below can never observe its name and the advertised
#: message is unreachable. ``lifelines`` and ``graphviz`` were mapped here while
#: the package imported neither (cycle-6 F-004); the gate is
#: ``tests/test_supported_versions.py::
#: test_every_optional_dependency_extra_has_a_real_import_site``.
#: Only an import that runs unguarded while a lazy module is being imported
#: reaches that handler, which today is matplotlib. seaborn is imported inside
#: the plotting functions, which read its entry through
#: ``visualization.utils._import_seaborn``. scikit-learn and numba are imported
#: inside ``try`` blocks -- their features fail with their own message
#: (``scikit-learn is required for ...``) or fall back to the reference kernels --
#: so entries for them could never be shown and were removed (audit S7-16).
_OPTIONAL_DEPENDENCY_EXTRAS = {
    "matplotlib": "plot",
    "seaborn": "plot",
}


def __getattr__(name: str):
    module_name = _LAZY_ATTR_TO_MODULE.get(name)
    if module_name is None:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
    try:
        module = importlib.import_module(f".{module_name}", __name__)
    except ImportError as exc:
        missing = (exc.name or "").split(".")[0]
        extra = _OPTIONAL_DEPENDENCY_EXTRAS.get(missing)
        if extra is None:
            raise
        raise ImportError(
            f"{missing} is required for {name}; install the '{extra}' extra "
            f"(pip install 'traminex[{extra}]')"
        ) from exc
    try:
        value = getattr(module, name)
    except AttributeError as exc:
        raise AttributeError(
            f"module '{module.__name__}' does not define expected export '{name}'"
        ) from exc
    globals()[name] = value
    # Importing one lazy module may import another whose name is also a public
    # callable (for example ``format`` imports the ``seqformat`` submodule).
    _restore_loaded_submodule_collisions()
    return value


#: Lazy modules that import an optional distribution when they are imported, so
#: none of their exports resolves while it is absent. ``help(traminex)`` and
#: ``inspect.getmembers(traminex)`` call ``getattr`` on every name ``dir()``
#: returns and skip only ``AttributeError``; the informative ``ImportError``
#: above made both fail on a core install (audit S7-05). ``__dir__`` therefore
#: omits these exports while the distribution cannot be found. ``__all__`` keeps
#: every manifest name, so ``from traminex import *`` still raises that
#: ``ImportError`` without the extra. ``tests/test_optional_extra_diagnostics.py``
#: checks that exactly the unresolvable names are omitted for each extra.
_LAZY_MODULE_OPTIONAL_IMPORTS: Dict[str, tuple[str, ...]] = {"plotting": ("matplotlib",)}


def _distribution_importable(name: str) -> bool:
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):
        return False


def __dir__():
    unavailable = {
        export
        for module_name, required in _LAZY_MODULE_OPTIONAL_IMPORTS.items()
        if not all(_distribution_importable(name) for name in required)
        for export in _LAZY_EXPORTS[module_name]
        if export not in globals()
    }
    return sorted((set(globals()) | set(__all__)) - unavailable)


class _PublicFacadeModule(types.ModuleType):
    """Prevent importlib from replacing colliding manifest exports with modules."""

    def __setattr__(self, name: str, value) -> None:
        if name in _MODULE_NAME_EXPORT_COLLISIONS and isinstance(value, types.ModuleType):
            public_callable = getattr(value, name, None)
            if public_callable is not None:
                super().__setattr__(f"_loaded_submodule_{name}", value)
                value = public_callable
        super().__setattr__(name, value)


sys.modules[__name__].__class__ = _PublicFacadeModule
