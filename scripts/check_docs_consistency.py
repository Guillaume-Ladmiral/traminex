#!/usr/bin/env python3
"""Lightweight documentation consistency checks.

This script is intentionally strict for active docs and lenient for archived
historical files under archive/.
"""

from __future__ import annotations

import ast
import importlib
import inspect
import json
import re
import subprocess
import sys
import textwrap
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


_INCLUDE_DIRECTIVE = re.compile(r"\{\{<\s*include\s+([^\s>]+)\s*>\}\}")


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _read_rendered(path: Path) -> str:
    """Read a doc as Quarto renders it, following ``{{< include X >}}``.

    Seven of the ``docs/`` pairs are a ``.qmd`` wrapper around the ``.md`` that
    carries the prose. A content check that reads only the wrapper sees front
    matter and an include directive, so it has to resolve the include or it
    reports the prose as missing. Reading what is actually rendered is also
    stricter than reading the file: it is the page a reader meets.
    """
    text = _read(path)
    for name in _INCLUDE_DIRECTIVE.findall(text):
        included = path.parent / name.strip()
        if included.is_file():
            text += "\n" + _read(included)
    return text


def _extract_pyproject_version(pyproject: Path) -> str:
    txt = _read(pyproject)
    m = re.search(r'(?m)^version\s*=\s*"([^"]+)"\s*$', txt)
    if not m:
        raise RuntimeError("Could not find [project] version in pyproject.toml")
    return m.group(1)


def _init_version_is_single_sourced(init_py: Path) -> list[str]:
    """Verify ``traminex/__init__.py`` derives ``__version__`` from metadata.

    The version must not be a hand-maintained string literal (the H1-class drift
    this guard exists to prevent); it must be read from
    ``importlib.metadata.version`` so ``pyproject.toml`` stays the single source.
    """
    src = _read(init_py)
    problems: list[str] = []
    if re.search(r'(?m)^\s*__version__\s*=\s*["\']', src):
        problems.append(
            "traminex/__init__.py hardcodes a __version__ string literal; "
            "the version must be sourced from importlib.metadata.version()"
        )
    if "metadata.version(" not in src:
        problems.append(
            "traminex/__init__.py does not source __version__ from " "importlib.metadata.version()"
        )
    return problems


def _check_doc_version_strings(pyproject_version: str) -> list[str]:
    """Enforce that doc "current version" declarations match pyproject.toml.

    README / CONTRIBUTING / INDEX each declare the package version in a fixed,
    parseable form; this keeps them from drifting away from the single declared
    source (the audit's H1 defect). If a marker is reworded so its pattern no
    longer matches, this fails loudly rather than letting the version go
    unchecked.
    """
    checks = [
        (
            REPO / "README.md",
            r"Current package version: `([^`]+)`",
            "README.md current-version line",
        ),
        (
            REPO / "CONTRIBUTING.md",
            r"(?m)^\*\*Version\*\*:\s*(\S+)\s*$",
            "CONTRIBUTING.md Version footer",
        ),
        (
            REPO / "docs" / "INDEX.md",
            r"pyproject\.toml \(package metadata, version ([^;]+);",
            "docs/INDEX.md pyproject metadata line",
        ),
    ]
    problems: list[str] = []
    for path, pattern, label in checks:
        if not path.exists():
            problems.append(f"Missing doc for version check: {path.relative_to(REPO)}")
            continue
        m = re.search(pattern, _read(path))
        if not m:
            problems.append(
                f"{label}: version marker not found (expected to declare {pyproject_version})"
            )
            continue
        found = m.group(1).strip()
        if found != pyproject_version:
            problems.append(f"{label}: declares {found} but pyproject.toml is {pyproject_version}")
    return problems


def _extract_literal_str_list(node: ast.AST) -> list[str] | None:
    if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
        out: list[str] = []
        for elt in node.elts:
            if isinstance(elt, ast.Constant) and isinstance(elt.value, str):
                out.append(elt.value)
            else:
                return None
        return out
    return None


def _extract_assigned_str_list(path: Path, var_name: str) -> list[str]:
    src = _read(path)
    mod = ast.parse(src, filename=str(path))
    for node in mod.body:
        value: ast.AST | None = None
        if isinstance(node, ast.Assign):
            if any(isinstance(t, ast.Name) and t.id == var_name for t in node.targets):
                value = node.value
        elif isinstance(node, ast.AnnAssign):
            if isinstance(node.target, ast.Name) and node.target.id == var_name:
                value = node.value
        if value is None:
            continue
        parsed = _extract_literal_str_list(value)
        if parsed is not None:
            return parsed
    raise RuntimeError(f"Could not parse {var_name} list from {path}")


def _runtime_all_uses_manifest(init_py: Path) -> bool:
    src = _read(init_py)
    mod = ast.parse(src, filename=str(init_py))
    has_manifest_import = False
    has_manifest_all = False
    for node in mod.body:
        if isinstance(node, ast.ImportFrom) and node.module == "api_manifest":
            if any(alias.name == "PUBLIC_API" for alias in node.names):
                has_manifest_import = True
        if isinstance(node, ast.Assign):
            if not any(isinstance(t, ast.Name) and t.id == "__all__" for t in node.targets):
                continue
            if (
                isinstance(node.value, ast.Call)
                and isinstance(node.value.func, ast.Name)
                and node.value.func.id == "list"
                and len(node.value.args) == 1
                and isinstance(node.value.args[0], ast.Name)
                and node.value.args[0].id == "PUBLIC_API"
            ):
                has_manifest_all = True
    return has_manifest_import and has_manifest_all


def _extract_doc_api_refs(path: Path) -> set[str]:
    txt = _read(path)
    refs: set[str] = set()

    for m in re.finditer(r"\btx\.([A-Za-z_][A-Za-z0-9_]*)\b", txt):
        refs.add(m.group(1))

    for m in re.finditer(r"from\s+traminex\s+import\b\s*", txt):
        refs.update(_traminex_import_names(txt, m.start(), m.end()))

    return refs


def _traminex_import_names(txt: str, start: int, names_start: int) -> set[str]:
    """The names one ``from traminex import`` statement imports, read with ``ast``.

    The statement runs to the end of its line (continued past a line ending in a
    backslash), or to the closing parenthesis of a parenthesised list: each ``)`` is tried
    in turn until the text parses, so one inside a comment does not end it (batch-16
    review M-2). Aliases (``name as alias``) and comments are not API names
    (TX-AUD-20261019-055: a multi-line block was read as the name '('). Text that does
    not parse, such as an inline mention followed by prose, keeps the plain reading:
    the leading identifier of each comma-separated item on the line.
    """
    if txt.startswith("(", names_start):
        ends = [names_start + match.end() for match in re.finditer(r"\)", txt[names_start:])]
        ends = ends[:50] or [len(txt)]
    else:
        ends = []
        position = names_start
        while True:
            newline = txt.find("\n", position)
            end = newline if newline != -1 else len(txt)
            ends.append(end)
            if newline == -1 or not txt[start:end].rstrip().endswith("\\"):
                break
            position = newline + 1
    for end in ends:
        try:
            module = ast.parse(txt[start:end].strip())
        except SyntaxError:
            continue
        if len(module.body) == 1 and isinstance(module.body[0], ast.ImportFrom):
            return {alias.name for alias in module.body[0].names if alias.name != "*"}
    names = set()
    for item in txt[names_start : ends[0]].split("#", 1)[0].split(","):
        match = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)", item)
        if match:
            names.add(match.group(1))
    return names


#: Public functions whose ``**kwargs`` are passed on to another function, a renderer or
#: matplotlib, so docs-check cannot name a closed keyword set for them. ``seqdist``,
#: ``seqdef`` and the sequence builders below validate theirs and are checked
#: (TX-AUD-20261019-054); a test requires the three lists to equal the public
#: ``**kwargs`` callables.
_KWARGS_FORWARDERS = (
    "omloc",
    "plot_SAMM",
    "plot_eseq",
    "plot_seqelist",
    "plot_seqimplic",
    "seqBIC",
    "seqCompare",
    "seqLRT",
    "seqdiffplot",
    "seqdistmc",
    "seqeordplot",
    "seqindic_dyn",
    "seqplot",
    "seqplotMD",
    "seqplot_rf",
    "seqplot_tentrop",
    "seqplot_tentrop_m",
    "seqpolyads",
    "seqrep",
    "seqrep_grp",
    "seqrfplot",
    "seqrplot",
    "seqsplot",
    "seqsubm",
)

#: Builders whose ``**kwargs`` take only ``SEQ_KWARGS`` (``validate_seq_kwargs`` refuses any
#: other name; batch-16 review M-1).
_KWARGS_SEQUENCE_BUILDERS = ("HSPELL_to_STS", "TSE_to_STS", "seqdef_from_long")


def _popped_keywords(function: object, mapping: str) -> set[str]:
    """The public names ``function`` pops from its ``mapping`` argument, read from its source.

    A ``for`` loop over a literal tuple whose body pops from ``validation_kwargs`` counts
    too (``seqdist``'s transport keys). Reading the names, rather than copying them here,
    keeps one definition; a restructured function makes the set shrink and docs-check fail
    loudly, never pass a wrong keyword.
    """
    tree = ast.parse(textwrap.dedent(inspect.getsource(function)))
    names: set[str] = set()
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "pop"
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == mapping
            and node.args
            and isinstance(node.args[0], ast.Constant)
            and isinstance(node.args[0].value, str)
        ):
            names.add(node.args[0].value)
        if isinstance(node, ast.For) and isinstance(node.iter, ast.Tuple):
            pops_validation = any(
                isinstance(inner, ast.Attribute)
                and inner.attr == "pop"
                and isinstance(inner.value, ast.Name)
                and inner.value.id == "validation_kwargs"
                for inner in ast.walk(node)
            )
            if pops_validation:
                names.update(
                    element.value
                    for element in node.iter.elts
                    if isinstance(element, ast.Constant) and isinstance(element.value, str)
                )
    return {name for name in names if not name.startswith("_")}


def _seqdist_accepted_keywords(method: str, parameters) -> set[str] | None:
    """The keywords ``seqdist(method=method)`` accepts, or None when its set is not closed.

    The run time accepts: the keywords ``seqdist`` pops itself; for one of R's 16 methods,
    every R argument of ``R_SEQDIST_ARGUMENT_USERS`` (an argument the method does not use
    is dropped, as R ignores it, TX-AUD-20261010-071) and the names
    ``resolve_r_seqdist_arguments`` pops itself (``check_costs``, ``check_max_size``,
    ``opt_args``); then the per-method list of ``validate_pairwise_kwargs``. A method that
    list does not cover is not checked (batch-16 review I-1).
    """
    from traminex import distances
    from traminex.distances import _validation

    m = method.upper()
    allowed = _validation._allowed_pairwise_kwargs(dispatcher=True).get(m)
    if allowed is None:
        return None
    accepted = set(parameters) | _popped_keywords(distances.seqdist, "kwargs") | allowed
    if m in _validation.R_SEQDIST_METHODS:
        accepted |= set(_validation.R_SEQDIST_ARGUMENT_USERS)
        accepted |= _popped_keywords(_validation.resolve_r_seqdist_arguments, "kwargs")
    return accepted


def _literal_argument(node: ast.Call, keyword: str, position: int | None) -> str | None:
    """The string literal given for ``keyword`` (by name, or at ``position`` unless the
    parameter is keyword-only), else None."""
    for item in node.keywords:
        if item.arg == keyword:
            value = item.value
            break
    else:
        if position is None or len(node.args) <= position:
            return None
        if any(isinstance(arg, ast.Starred) for arg in node.args[: position + 1]):
            return None
        value = node.args[position]
    return value.value if isinstance(value, ast.Constant) and isinstance(value.value, str) else None


def _unknown_validated_keywords(target: object, node: ast.Call, signature) -> list[str]:
    """Keywords outside the closed set of a ``**kwargs`` function that validates its own.

    ``seqdist`` accepts the keywords ``_seqdist_accepted_keywords`` lists; without a
    literal method (an expansion may carry one), or for a method its per-method list does
    not cover, nothing is checked. ``seqdef`` with the STS format (the default, unless an
    expansion may carry another) accepts no keyword outside its signature (format.py); the
    sequence builders accept their signature and ``SEQ_KWARGS``.
    """
    import traminex as tx

    explicit = [item.arg for item in node.keywords if item.arg is not None]
    # A ``**options`` or ``*args`` expansion may carry the method or the format: nothing
    # is assumed for it (batch-16 review I-1).
    expanded = any(item.arg is None for item in node.keywords) or any(
        isinstance(arg, ast.Starred) for arg in node.args
    )
    if target is tx.seqdist:
        method = _literal_argument(node, "method", 1)
        if method is None and "method" not in explicit and len(node.args) < 2 and not expanded:
            method = "OM"
        if method is None:
            return []
        accepted = _seqdist_accepted_keywords(method, signature.parameters)
        if accepted is None:
            return []
        return [
            f"keyword '{name}' is not accepted by seqdist(method={method!r})"
            for name in explicit
            if name not in accepted
        ]
    if target is tx.seqdef:
        given = any(item.arg == "informat" for item in node.keywords)
        if not given and expanded:
            return []
        informat = _literal_argument(node, "informat", None) if given else "STS"
        if informat != "STS":
            return []
        return [
            f"keyword '{name}' is not accepted by seqdef with the STS format"
            for name in explicit
            if name not in signature.parameters
        ]
    if target in tuple(getattr(tx, name) for name in _KWARGS_SEQUENCE_BUILDERS):
        from traminex._sequence_builders import SEQ_KWARGS

        return [
            f"keyword '{name}' is not accepted (validate_seq_kwargs takes {sorted(SEQ_KWARGS)})"
            for name in explicit
            if name not in signature.parameters and name not in SEQ_KWARGS
        ]
    return []


def _check_python_api_calls(path: Path) -> list[str]:
    """Validate explicit traminex API calls in Python code fences.

    The check binds each documented call shape against the installed traminex
    signature without executing it.  Dynamic ``*args``/``**kwargs`` calls cannot
    be bound statically, but their explicit keywords are still checked. For a
    function that takes ``**kwargs`` any keyword binds, so ``seqdist``, ``seqdef`` and the
    ``_KWARGS_SEQUENCE_BUILDERS``, which refuse unknown keywords at run time, are checked
    against their closed keyword sets; the forwarders in ``_KWARGS_FORWARDERS`` are not
    (TX-AUD-20261019-054).
    """
    import traminex as tx

    txt = _read(path)
    display_path = path.relative_to(REPO) if path.is_relative_to(REPO) else path
    problems: list[str] = []
    # Plain Markdown fences and executable Quarto chunks (audit S7-21).
    fence_pattern = re.compile(r"```(?:python|\{python[^}\n]*\})[^\n]*\n(.*?)```", flags=re.S)

    for match in fence_pattern.finditer(txt):
        code = match.group(1)
        first_doc_line = txt.count("\n", 0, match.start(1)) + 1
        try:
            module = ast.parse(code, filename=str(path))
        except SyntaxError as exc:
            line = first_doc_line + (exc.lineno or 1) - 1
            problems.append(f"{display_path}:{line}: invalid Python syntax: {exc.msg}")
            continue

        module_aliases: dict[str, object] = {"tx": tx}
        imported_aliases: dict[str, tuple[str, object | None]] = {}
        for statement in module.body:
            if isinstance(statement, ast.Import):
                for alias in statement.names:
                    if alias.name == "traminex" or alias.name.startswith("traminex."):
                        importlib.import_module(alias.name)
                        if not alias.asname:
                            # ``import traminex.submodule`` binds the root
                            # package name.
                            module_aliases[alias.name.split(".", maxsplit=1)[0]] = tx
                            continue
                        # ``import traminex.a.b as c`` compiles to IMPORT_FROM,
                        # i.e. ``getattr(getattr(traminex, "a"), "b")``, with
                        # ``sys.modules`` only as the fallback. For the public
                        # functions named like their submodule (``seqtree``,
                        # ``disstree``, ...) that binds the function, not the
                        # submodule (audit S7-06).
                        bound: object = tx
                        qualified = "traminex"
                        for part in alias.name.split(".")[1:]:
                            qualified = f"{qualified}.{part}"
                            try:
                                bound = getattr(bound, part)
                            except AttributeError:
                                bound = sys.modules[qualified]
                        module_aliases[alias.asname] = bound
            elif (
                isinstance(statement, ast.ImportFrom)
                and statement.module is not None
                and (statement.module == "traminex" or statement.module.startswith("traminex."))
            ):
                imported = importlib.import_module(statement.module)
                for alias in statement.names:
                    if alias.name == "*":
                        continue
                    target = getattr(imported, alias.name, None)
                    local_name = alias.asname or alias.name
                    imported_aliases[local_name] = (alias.name, target)

        for node in ast.walk(module):
            if not isinstance(node, ast.Call):
                continue

            api_name: str | None = None
            target: object | None = None
            if (
                isinstance(node.func, ast.Attribute)
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id in module_aliases
            ):
                api_name = node.func.attr
                target = getattr(module_aliases[node.func.value.id], api_name, None)
            elif isinstance(node.func, ast.Name) and node.func.id in imported_aliases:
                api_name, target = imported_aliases[node.func.id]

            if api_name is None:
                continue
            if target is None or not callable(target):
                line = first_doc_line + node.lineno - 1
                problems.append(
                    f"{display_path}:{line}: traminex.{api_name} is not a callable traminex API"
                )
                continue
            try:
                signature = inspect.signature(target)
            except (TypeError, ValueError):
                continue

            explicit_keywords = {kw.arg: object() for kw in node.keywords if kw.arg is not None}
            has_dynamic_arguments = any(isinstance(arg, ast.Starred) for arg in node.args) or any(
                kw.arg is None for kw in node.keywords
            )
            try:
                if has_dynamic_arguments:
                    accepts_var_keywords = any(
                        param.kind is inspect.Parameter.VAR_KEYWORD
                        for param in signature.parameters.values()
                    )
                    if not accepts_var_keywords:
                        unexpected = next(
                            (
                                name
                                for name in explicit_keywords
                                if name not in signature.parameters
                            ),
                            None,
                        )
                        if unexpected is not None:
                            raise TypeError(f"got an unexpected keyword argument '{unexpected}'")
                else:
                    signature.bind(
                        *(object() for _ in node.args),
                        **explicit_keywords,
                    )
            except TypeError as exc:
                line = first_doc_line + node.lineno - 1
                problems.append(f"{display_path}:{line}: traminex.{api_name}{signature}: {exc}")
                continue
            for message in _unknown_validated_keywords(target, node, signature):
                line = first_doc_line + node.lineno - 1
                problems.append(f"{display_path}:{line}: traminex.{api_name}: {message}")

    return problems


def _display_path(path: Path) -> Path:
    return path.relative_to(REPO) if path.is_relative_to(REPO) else path


_AUDIT_RUN_CITATION = re.compile(r"reports/audit/([A-Za-z0-9_.-]+)/([A-Za-z0-9_./-]*)")


def _tracked_files(repo: Path, prefix: str) -> set[str] | None:
    """Return tracked paths under ``prefix``, or ``None`` outside a git checkout."""
    result = subprocess.run(
        ["git", "-C", str(repo), "ls-files", "-z", "--", prefix],
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        return None
    return {name.decode() for name in result.stdout.split(b"\0") if name}


def _check_cited_audit_runs_are_tracked(paths: list[Path], tracked: set[str]) -> list[str]:
    """Require every cited ``reports/audit/<run>/...`` path to be in version control.

    ``.gitignore`` excludes ``reports/audit/`` by default and ``docs/DOCS_POLICY.md``
    tracks only a run's report and evidence index, so a document can cite a run
    directory that exists only on the author's machine (audit S7-13). A cited
    file must be tracked; a cited directory must contain a tracked file.
    """
    problems: list[str] = []
    for path in paths:
        if not path.exists():
            continue
        for match in _AUDIT_RUN_CITATION.finditer(_read(path)):
            rest = match.group(2).rstrip(".")
            cited = f"reports/audit/{match.group(1)}/{rest}"
            if rest and not rest.endswith("/"):
                resolved = cited in tracked
            else:
                directory = cited if cited.endswith("/") else cited + "/"
                resolved = any(name.startswith(directory) for name in tracked)
            if not resolved:
                problems.append(
                    f"{_display_path(path)} cites {cited}, which is not tracked; cite the "
                    "run's tracked report or track it as docs/DOCS_POLICY.md describes"
                )
    return problems


def _vendored_oracle_versions(repo: Path = REPO) -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for package in ("TraMineR", "TraMineRextras"):
        try:
            description = _read(repo / "SOURCES" / package / "DESCRIPTION")
        except OSError:
            versions[package] = None
            continue
        match = re.search(r"(?m)^Version:\s*(\S+)\s*$", description)
        versions[package] = match.group(1) if match else None
    return versions


def _check_oracle_versions_named(path: Path, repo: Path = REPO) -> list[str]:
    """Require an active parity doc to name the vendored oracle versions.

    ``docs/DOCS_POLICY.md`` says active parity docs must use the current oracle,
    but nothing checked it (audit S7-21). The versions are read from the
    vendored ``DESCRIPTION`` files, so moving the oracle fails this check until
    the document follows.
    """
    text = _read(path)
    problems: list[str] = []
    for package, version in _vendored_oracle_versions(repo).items():
        if version is None:
            problems.append(f"cannot read the vendored {package} version from SOURCES/")
        elif f"{package} {version}" not in text:
            problems.append(
                f"{_display_path(path)} does not name the vendored oracle "
                f"{package} {version} (docs/DOCS_POLICY.md)"
            )
    return problems


def _check_active_command_policy(path: Path) -> list[str]:
    """Reject executable Windows launch syntax in maintained documentation."""
    text = _read(path)
    patterns = (
        r"(?im)^```(?:powershell|pwsh)\b",
        r"(?im)^\s*(?:powershell|pwsh)(?:\s|$)",
        r"(?i)(?:^|[\s`])\S+\.ps1(?:\s|$)",
        r"(?i)(?:^|[\s`])python\.exe(?:\s|$)",
    )
    if any(re.search(pattern, text) for pattern in patterns):
        return [
            f"{_display_path(path)} contains an active Windows/PowerShell launch command; "
            "maintained commands must use Ubuntu shell syntax"
        ]
    return []


def _check_hardcoded_module_claims(path: Path) -> list[str]:
    """Reject the known transient 'all N package modules' documentation claim."""
    # "package", "source" or both ("all 78 package source modules", S7-14).
    if re.search(r"(?i)\ball\s+\d+\s+(?:(?:package|source)\s+){0,2}modules\b", _read(path)):
        return [
            f"{_display_path(path)} contains a hard-coded package-module count; "
            "describe the configured module set or derive the count mechanically"
        ]
    return []


def _check_stale_active_status(path: Path) -> list[str]:
    """Reject the obsolete pre-R6 audit status and require its evidence boundary."""
    text = _read(path)
    has_pre_r6_status = "Pre-R6 full audit" in text
    stale = re.search(r"(?is)Pre-R6 full audit.{0,160}\bnot started\b", text)
    if stale or (
        has_pre_r6_status and ("2026-08-13" not in text or "final verification" not in text)
    ):
        return [
            f"{_display_path(path)} contains a stale or unbounded pre-R6 full-audit status; "
            "name the 2026-08-13 audit and make remediation completion depend on final verification"
        ]
    return []


def _check_unqualified_speedup_claims(path: Path) -> list[str]:
    """Reject universal or context-free numeric performance promises."""
    if re.search(r"(?i)(?:massive|expected)\s+speedups?", _read(path)):
        return [
            f"{_display_path(path)} contains an unqualified speedup claim; "
            "refer readers to the reproducible benchmark and qualify workload, hardware, and versions"
        ]
    return []


def _check_api_consistency_missing_keywords(path: Path) -> list[str]:
    """Keep the documented missing-state keywords bound to live signatures."""
    import traminex as tx

    expected_names = {
        "seqindic",
        "seqistatd",
        "seqnum",
        "seqstatd",
        "seqtabstocc",
        "seqtrate",
    }
    # Dispatchers that accept a missing-state keyword through ``**kwargs`` are
    # invisible to a signature check, and were therefore undocumented while still
    # being *mandatory* on missing-bearing data (cycle-4b F-007). Derive them from
    # the live package rather than hardcoding, so a new one cannot drift in
    # undocumented.
    kwargs_dispatchers = _kwargs_missing_keyword_dispatchers()
    expected_names |= set(kwargs_dispatchers)
    declarations = dict(
        re.findall(
            r"(?m)^\s*-\s+([A-Za-z_][A-Za-z0-9_]*)\(x,\s*" r"(include_missing|with_missing)=",
            _read_rendered(path),
        )
    )
    problems: list[str] = []
    if set(declarations) != expected_names:
        missing = sorted(expected_names - set(declarations))
        extra = sorted(set(declarations) - expected_names)
        problems.append(
            f"{_display_path(path)} missing-state declaration set differs from the live contract; "
            f"missing={missing}, extra={extra}"
        )
    for name, keyword in declarations.items():
        parameters = inspect.signature(getattr(tx, name)).parameters
        if keyword in parameters:
            continue
        if name in kwargs_dispatchers:
            if keyword not in kwargs_dispatchers[name]:
                problems.append(
                    f"{_display_path(path)} documents {name}({keyword}=...) but the "
                    f"dispatcher consumes {sorted(kwargs_dispatchers[name])}"
                )
            elif keyword not in (inspect.getdoc(getattr(tx, name)) or ""):
                problems.append(
                    f"{_display_path(path)} documents {name}({keyword}=...) but the "
                    f"function's own docstring never mentions it"
                )
            continue
        problems.append(f"{_display_path(path)} documents unsupported {name} keyword: {keyword}")
    return problems


def _check_every_promoted_r_object_has_a_disposition() -> list[str]:
    """No public export may name a defined R object without a decision row.

    traminex promotes several objects that R defines but does not export. The
    decision CSV's ``r_kind`` domain was ``export``/``s3``, which the NAMESPACE can
    describe, so a promoted internal was structurally unrepresentable and seven
    public analytical surfaces carried no authoritative disposition (cycle-4b
    F-004). ``r_kind = internal`` closes that, and this check keeps it closed:
    every registered promoted internal must have a row, with that kind, mapping to
    a real public export.
    """
    import csv

    from traminex.api_manifest import PUBLIC_API_SET

    sys.path.insert(0, str(REPO))
    from scripts._r_namespace import promoted_internals

    problems: list[str] = []
    decisions = REPO / "docs" / "API_PARITY_DECISIONS.csv"
    with decisions.open(newline="", encoding="utf-8") as handle:
        rows = {row["r_name"]: row for row in csv.DictReader(handle)}

    try:
        verified = promoted_internals(REPO / "SOURCES")
    except Exception as exc:  # noqa: BLE001 - surfaced as a gate failure, not a crash
        return [f"promoted-internal registry no longer matches the vendored source: {exc}"]

    for name, package in sorted(verified.items()):
        row = rows.get(name)
        if row is None:
            problems.append(
                f"docs/API_PARITY_DECISIONS.csv has no row for the promoted R internal {name!r}; "
                "add one with r_kind=internal"
            )
            continue
        if row["r_kind"] != "internal":
            problems.append(
                f"docs/API_PARITY_DECISIONS.csv row {name!r} has r_kind={row['r_kind']!r}; "
                "a promoted R internal must be recorded as r_kind=internal"
            )
        if row["r_package"] != package:
            problems.append(
                f"docs/API_PARITY_DECISIONS.csv row {name!r} names package "
                f"{row['r_package']!r}; the object is defined in {package!r}"
            )
        if row["python_name"] not in PUBLIC_API_SET:
            problems.append(
                f"docs/API_PARITY_DECISIONS.csv row {name!r} maps to "
                f"{row['python_name']!r}, which is not a public export"
            )
    return problems


def _check_backend_switch_domains() -> list[str]:
    """Pin the two backend switches' domains and keep them distinct in the docs.

    ``docs/USAGE_new_apis.{md,qmd}`` presented ``TRAMINEX_BINDING`` as giving
    "automatic pybind11/Numba/Python selection". It accepts exactly three values
    and cannot disable Numba, which is a separate, default-active backend
    (cycle-4b F-015/F-008). Both domains are checked against the live package so
    the prose cannot drift from the code again.
    """
    from traminex.distances import python_kernels

    problems: list[str] = []

    binding_source = _read(REPO / "traminex" / "distances" / "__init__.py")
    if 'TRAMINEX_BINDING must be one of: auto, python, pybind11"' not in binding_source:
        problems.append(
            "traminex/distances/__init__.py no longer states the three-value "
            "TRAMINEX_BINDING domain; update the docs check and the documentation"
        )
    if tuple(python_kernels.NUMBA_MODES) != ("auto", "off"):
        problems.append(
            f"TRAMINEX_NUMBA domain changed to {python_kernels.NUMBA_MODES}; "
            "update docs/DISTANCES.md and docs/USAGE_new_apis.{md,qmd}"
        )

    for name in ("USAGE_new_apis.md", "USAGE_new_apis.qmd"):
        text = _read_rendered(REPO / "docs" / name)
        if "pybind11/Numba/Python selection" in text:
            problems.append(
                f"docs/{name} still presents TRAMINEX_BINDING as selecting Numba; "
                "it accepts only auto, python, pybind11"
            )
        if "TRAMINEX_NUMBA" not in text:
            problems.append(f"docs/{name} does not name the TRAMINEX_NUMBA switch")
    # DISTANCES.qmd is the file Quarto renders into the published site, so an
    # omission there is the one a reader actually meets (cycle-6 F-005).
    for name in ("DISTANCES.md", "DISTANCES.qmd"):
        if "TRAMINEX_NUMBA" not in _read_rendered(REPO / "docs" / name):
            problems.append(f"docs/{name} does not name the TRAMINEX_NUMBA switch")
    problems.extend(_check_doc_pair_switch_symmetry())
    return problems


#: Environment switches that change observable behaviour. A ``.md``/``.qmd`` pair
#: covering the same topic must not name one of these in only one member.
_DOCUMENTED_SWITCHES = (
    "TRAMINEX_BINDING",
    "TRAMINEX_NUMBA",
    "TRAMINEX_PARITY",
    "TRAMINEX_CALIBRATED_PARAMS",
    "GVPATH",
)


def _check_doc_pair_switch_symmetry() -> list[str]:
    """Require both members of a standalone doc pair to name the same switches.

    ``docs/INDEX.md`` declares that for a paired file "the ``.qmd`` is the
    Quarto-site source and the ``.md`` is the GitHub-readable companion; they
    cover the same topic". When this check was written three of seven pairs were
    include-wrappers and four were standalone duplicates, one of which had
    drifted on a user-visible backend switch while every gate passed (cycle-6
    F-005). Since then every paired ``.qmd`` has become an include-wrapper, so
    today the loop below skips every pair; it stays to catch a standalone pair
    reintroduced later (audit S7-21).

    Presence, not count, is compared: a file may reasonably mention a switch
    more or fewer times than its companion, but naming it in only one member
    means the rendered site and the readable companion disagree about what
    changes a user's result.
    """
    docs = REPO / "docs"
    problems: list[str] = []
    for qmd in sorted(docs.glob("*.qmd")):
        markdown = qmd.with_suffix(".md")
        if not markdown.is_file():
            continue
        qmd_text = _read(qmd)
        included = {name.strip() for name in _INCLUDE_DIRECTIVE.findall(qmd_text)}
        if markdown.name in included:
            continue  # An include-wrapper cannot drift from the file it includes.
        markdown_text = _read(markdown)
        for switch in _DOCUMENTED_SWITCHES:
            in_md = switch in markdown_text
            in_qmd = switch in qmd_text
            if in_md == in_qmd:
                continue
            present, absent = (markdown.name, qmd.name) if in_md else (qmd.name, markdown.name)
            problems.append(
                f"docs/{present} names {switch} but its standalone companion "
                f"docs/{absent} does not; the pair documents the same topic, so "
                "either document the switch in both or make one an include-wrapper"
            )
    return problems


def _check_parameter_rename_table(path: Path) -> list[str]:
    """Keep the documented semantic parameter renames bound to live signatures.

    A ported function's argument names are part of its public contract, and an
    ``EXACT`` parity row implies a user can call it as documented in R. The
    mechanical rules (dots, camelCase, keyword collisions) are self-evident; the
    renames no rule predicts are only discoverable from this table, so the table
    must not drift away from the code (cycle-4b F-006).

    The R side of each row is verified against the live R namespaces by the parity
    audit, not here: this gate runs without R.
    """
    import traminex as tx

    text = _read_rendered(path)
    problems: list[str] = []
    if "Parameter naming: how R argument names map to Python" not in text:
        problems.append(f"{_display_path(path)} has no parameter-naming section")
        return problems
    for rule in ("dots become underscores", "camelCase becomes snake_case", "trailing underscore"):
        if rule not in text:
            problems.append(f"{_display_path(path)} parameter-naming section omits: {rule}")

    rows = re.findall(
        r"(?m)^\s*\|\s*`([A-Za-z_][A-Za-z0-9_.]*)`\s*\|\s*`([A-Za-z_][A-Za-z0-9_.]*)`"
        r"\s*\|\s*`([A-Za-z_][A-Za-z0-9_]*)`\s*\|\s*$",
        text,
    )
    if not rows:
        problems.append(f"{_display_path(path)} parameter-naming section has no rename rows")
    for r_name, r_formal, python_parameter in rows:
        function = getattr(tx, r_name, None)
        if not callable(function):
            problems.append(
                f"{_display_path(path)} parameter-naming table names {r_name}, "
                "which is not a public traminex callable"
            )
            continue
        try:
            parameters = inspect.signature(function).parameters
        except (TypeError, ValueError):
            continue
        if python_parameter not in parameters:
            problems.append(
                f"{_display_path(path)} parameter-naming table maps {r_name}({r_formal}) to "
                f"'{python_parameter}', which is not a parameter of traminex.{r_name}"
            )
    return problems


def _kwargs_missing_keyword_dispatchers() -> dict[str, set[str]]:
    """Public ``**kwargs`` dispatchers that read a missing-state keyword.

    A dispatcher with no named missing-state parameter can still require one --
    ``seqdist`` raises ``ValueError: with_missing must be true when sequences
    contain missing values`` -- and neither its signature nor its stub reveals it.
    Enumerated from the manifest and the function source so the documentation gate
    covers this class instead of a hardcoded list.
    """
    import traminex as tx
    from traminex.api_manifest import PUBLIC_API

    keywords = ("with_missing", "include_missing")
    found: dict[str, set[str]] = {}
    for name in PUBLIC_API:
        function = getattr(tx, name, None)
        if not callable(function):
            continue
        try:
            parameters = inspect.signature(function).parameters
        except (TypeError, ValueError):
            continue
        if any(keyword in parameters for keyword in keywords):
            continue
        if not any(p.kind is inspect.Parameter.VAR_KEYWORD for p in parameters.values()):
            continue
        try:
            source = inspect.getsource(function)
        except (OSError, TypeError):
            continue
        consumed = {kw for kw in keywords if re.search(rf'["\']{kw}["\']', source)}
        if consumed:
            found[name] = consumed
    return found


def _check_nms_coverage_claims() -> list[str]:
    """Pin the soft-NMS source and bound claims to the routes tests exercise."""
    new_methods = _read(REPO / "docs" / "NEW_DISTANCE_METHODS.md")
    roadmap = _read(REPO / "ROADMAP.md")
    index = _read(REPO / "docs" / "INDEX.md")
    problems: list[str] = []
    if "NMSMSTSoftdistanceII.cpp" not in new_methods:
        problems.append(
            "docs/NEW_DISTANCE_METHODS.md does not identify NMSMSTSoftdistanceII.cpp "
            "as the non-identity proximity kernel"
        )
    if "Default and risk-weighted option routes are pinned" not in new_methods:
        problems.append(
            "docs/NEW_DISTANCE_METHODS.md does not bound its R-golden coverage claim to "
            "default and risk-weighted option routes"
        )
    if "default and risk-weighted option routes are golden-pinned" not in roadmap:
        problems.append("ROADMAP.md overstates or omits the bounded distance R-golden coverage")
    if "Default and risk-weighted option routes are pinned" not in index:
        problems.append("docs/INDEX.md overstates or omits the bounded distance R-golden coverage")
    return problems


def _check_generated_docs_artifacts(repo: Path = REPO) -> list[str]:
    """Reject generated docs only when tracked or not intentionally ignored."""
    docs_dir = repo / "docs"
    candidates = [docs_dir / "_site", docs_dir / ".quarto", *sorted(docs_dir.glob("*.html"))]
    candidates = [path for path in candidates if path.exists()]
    if not candidates:
        return []

    in_git_repo = (
        subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "--is-inside-work-tree"],
            capture_output=True,
            text=True,
            check=False,
        ).returncode
        == 0
    )
    ignore_rules = set()
    ignore_file = docs_dir / ".gitignore"
    if ignore_file.is_file():
        ignore_rules = {
            line.strip()
            for line in _read(ignore_file).splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        }

    problems: list[str] = []
    for path in candidates:
        relative = path.relative_to(repo)
        if in_git_repo:
            tracked = subprocess.run(
                ["git", "-C", str(repo), "ls-files", "--", relative.as_posix()],
                capture_output=True,
                text=True,
                check=False,
            ).stdout.splitlines()
            for tracked_path in tracked:
                problems.append(f"Tracked generated docs artifact: {tracked_path}")
            if tracked:
                continue
            ignored = (
                subprocess.run(
                    ["git", "-C", str(repo), "check-ignore", "-q", "--", relative.as_posix()],
                    capture_output=True,
                    check=False,
                ).returncode
                == 0
            )
        elif path == docs_dir / ".quarto":
            ignored = "/.quarto/" in ignore_rules
        elif path == docs_dir / "_site":
            ignored = "/_site/" in ignore_rules
        else:
            ignored = "/*.html" in ignore_rules
        if not ignored:
            problems.append(f"Generated docs artifact is not ignored: {relative.as_posix()}")
    return problems


def _check_audit_disclosures(repo: Path = REPO) -> list[str]:
    """Keep user-facing docs aligned with audited ordering and helper contracts."""
    problems: list[str] = []
    parity = _read(repo / "docs" / "API_PARITY.md")
    # The default alphabet is R's seqstatl order on every constructor (owner decision of
    # 2026-10-19); the section must say so, name the constructors and the explicit-alphabet rule.
    required_alphabet_phrases = (
        "`seqstatl` order",
        "alphabet=",
        "keeps the caller's order",
        "per-state rows",
        "seqformat",
        "seqdef_from_long",
        "TSE_to_STS",
        "HSPELL_to_STS",
        "events_to_state_sequences",
    )
    missing = [phrase for phrase in required_alphabet_phrases if phrase not in parity]
    if missing:
        problems.append(
            "docs/API_PARITY.md does not disclose the seqdef alphabet-order impact: "
            + ", ".join(missing)
        )

    seqtree = _read(repo / "docs" / "SEQTREE.md")
    if "seqtree_leaf_summary(result)" not in seqtree:
        problems.append("docs/SEQTREE.md does not show the canonical leaf-summary signature")
    if "per-sequence normalized entropy" not in seqtree:
        problems.append("docs/SEQTREE.md does not describe the implemented entropy feature")
    if "global entropy (same for all rows" in seqtree:
        problems.append("docs/SEQTREE.md retains the false global-entropy claim")

    examples = _read(repo / "docs" / "TraMineR_examples.qmd")
    if "seqtree_leaf_summary(X, res)" in examples:
        problems.append("docs/TraMineR_examples.qmd uses the deprecated leaf-summary call")
    return problems


def main() -> int:
    errors: list[str] = []

    pyproject_version = _extract_pyproject_version(REPO / "pyproject.toml")
    errors.extend(_init_version_is_single_sourced(REPO / "traminex" / "__init__.py"))
    errors.extend(_check_doc_version_strings(pyproject_version))

    active_docs = [
        REPO / "README.md",
        REPO / "CHANGELOG.md",
        REPO / "CONTRIBUTING.md",
        REPO / "ROADMAP.md",
        REPO / "docs" / "README.md",
        REPO / "docs" / "DEV_BUILD.md",
        REPO / "docs" / "INDEX.md",
        REPO / "docs" / "API_PARITY.md",
        REPO / "docs" / "API_CONSISTENCY.md",
        REPO / "docs" / "API_CONSISTENCY.qmd",
        REPO / "docs" / "CONSUMED_SURFACE.md",
        REPO / "docs" / "DOCS_POLICY.md",
        REPO / "docs" / "TraMineR_examples.qmd",
        REPO / "docs" / "USAGE_new_apis.md",
        REPO / "docs" / "USAGE_new_apis.qmd",
        REPO / "scripts" / "r_goldens" / "README.md",
        REPO / "reports" / "parity" / "README.md",
        REPO / "stubs" / "README.md",
    ]

    forbidden_version_tokens = [
        "0.1.1",
        "0.1.0",
        "0.2.x",
    ]
    forbidden_alias_tokens = [
        "states_to_events",
        "tx.seqturbulence(",
        "tx.turbulence(",
    ]
    forbidden_local_sibling_links = ["../../PROJECT_PY_GGPLOTSEQ/README.md"]

    for path in active_docs:
        if not path.exists():
            errors.append(f"Missing active doc: {path.relative_to(REPO)}")
            continue
        txt = _read(path)
        for tok in forbidden_version_tokens:
            if tok in txt:
                errors.append(f"{path.relative_to(REPO)} contains stale version token: {tok}")
        for tok in forbidden_alias_tokens:
            if tok in txt:
                errors.append(f"{path.relative_to(REPO)} contains removed alias token: {tok}")
        errors.extend(_check_active_command_policy(path))
        errors.extend(_check_hardcoded_module_claims(path))

    errors.extend(_check_stale_active_status(REPO / "ROADMAP.md"))
    for path in (
        REPO / "docs" / "API_PARITY.md",
        REPO / "docs" / "INDEX.md",
        REPO / "scripts" / "r_goldens" / "README.md",
    ):
        errors.extend(_check_oracle_versions_named(path))
    for path in (REPO / "ROADMAP.md", REPO / "examples" / "p0_features_demo.py"):
        errors.extend(_check_unqualified_speedup_claims(path))

    for path in (
        REPO / "docs" / "API_CONSISTENCY.md",
        REPO / "docs" / "API_CONSISTENCY.qmd",
    ):
        errors.extend(_check_api_consistency_missing_keywords(path))
        errors.extend(_check_parameter_rename_table(path))
    errors.extend(_check_backend_switch_domains())
    errors.extend(_check_every_promoted_r_object_has_a_disposition())
    errors.extend(_check_nms_coverage_claims())
    errors.extend(_check_audit_disclosures())

    # Every curated Markdown/Quarto source under docs/ is active. Validate
    # explicit top-level API call shapes even when a file is not in the smaller
    # version/alias-token inventory above.
    docs_dir = REPO / "docs"
    docs_for_python_calls = sorted(docs_dir.glob("*.md")) + sorted(docs_dir.glob("*.qmd"))
    for path in docs_for_python_calls:
        errors.extend(_check_python_api_calls(path))
        txt = _read(path)
        for link in forbidden_local_sibling_links:
            if link in txt:
                errors.append(
                    f"{path.relative_to(REPO)} contains non-portable sibling link: {link}"
                )

    tracked_audit = _tracked_files(REPO, "reports/audit")
    if tracked_audit is not None:
        cited_docs = sorted(set(active_docs) | set(docs_for_python_calls))
        errors.extend(_check_cited_audit_runs_are_tracked(cited_docs, tracked_audit))

    # Ignored render state is normal. Generated docs must remain untracked and
    # explicitly excluded from distribution inputs.
    errors.extend(_check_generated_docs_artifacts(REPO))

    manifest_path = REPO / "traminex" / "api_manifest.py"
    init_path = REPO / "traminex" / "__init__.py"
    stub_init_path = REPO / "stubs" / "traminex" / "__init__.pyi"
    orphan_stub_path = REPO / "stubs" / "traminex" / "visualization" / "export.pyi"

    if not manifest_path.exists():
        errors.append("Missing canonical API manifest: traminex/api_manifest.py")
        manifest_api = []
    else:
        try:
            manifest_api = _extract_assigned_str_list(manifest_path, "PUBLIC_API")
        except Exception as exc:
            errors.append(f"Failed to parse PUBLIC_API from traminex/api_manifest.py: {exc}")
            manifest_api = []

    if manifest_api and len(manifest_api) != len(set(manifest_api)):
        errors.append("traminex/api_manifest.py contains duplicate PUBLIC_API entries")

    if not _runtime_all_uses_manifest(init_path):
        errors.append("traminex/__init__.py __all__ is not wired to PUBLIC_API manifest")

    if stub_init_path.exists() and manifest_api:
        try:
            stub_all = _extract_assigned_str_list(stub_init_path, "__all__")
            if stub_all != manifest_api:
                errors.append(
                    "stubs/traminex/__init__.pyi __all__ does not match traminex/api_manifest.py PUBLIC_API"
                )
        except Exception as exc:
            errors.append(f"Failed to parse stubs/traminex/__init__.pyi __all__: {exc}")

    if orphan_stub_path.exists():
        errors.append("Orphan stub present: stubs/traminex/visualization/export.pyi")

    parity_doc = REPO / "docs" / "API_PARITY.md"
    if parity_doc.exists() and manifest_api:
        txt = _read(parity_doc)
        required_semantic_refs = {
            "API_PARITY_DECISIONS.csv",
            "API_PARITY_LEDGER.md",
            "traminex_semantic_parity.json",
        }
        missing_refs = sorted(ref for ref in required_semantic_refs if ref not in txt)
        if missing_refs:
            errors.append(
                "docs/API_PARITY.md missing semantic parity source/view references: "
                + ", ".join(missing_refs)
            )

        semantic_json = REPO / "reports" / "audit" / "traminex_semantic_parity.json"
        if not semantic_json.exists():
            errors.append(
                "Missing semantic parity report: reports/audit/traminex_semantic_parity.json"
            )
        else:
            try:
                semantic = json.loads(_read(semantic_json))
                semantic_current = semantic["target_public_api"]["current_count"]
                if semantic_current != len(manifest_api):
                    errors.append(
                        "semantic parity current API count mismatch: "
                        f"report={semantic_current} vs manifest={len(manifest_api)}"
                    )
                if semantic.get("source") != "docs/API_PARITY_DECISIONS.csv":
                    errors.append(
                        "semantic parity report does not name docs/API_PARITY_DECISIONS.csv "
                        "as its source"
                    )
                if semantic.get("decision_count") != len(semantic.get("decisions", [])):
                    errors.append("semantic parity decision count does not match its decision rows")
            except (KeyError, TypeError, ValueError) as exc:
                errors.append(
                    f"Failed to validate reports/audit/traminex_semantic_parity.json: {exc}"
                )

    # Validate that top-level API references in docs resolve to manifest entries.
    if manifest_api:
        manifest_set = set(manifest_api)
        docs_for_api_refs = [
            REPO / "README.md",
            REPO / "docs" / "API_PARITY.md",
            REPO / "docs" / "TraMineR_examples.qmd",
            REPO / "docs" / "USAGE_new_apis.md",
            REPO / "docs" / "USAGE_new_apis.qmd",
        ]
        examples_dir = REPO / "examples"
        if examples_dir.exists():
            docs_for_api_refs.extend(sorted(examples_dir.glob("*.ipynb")))
            docs_for_api_refs.extend(sorted(examples_dir.glob("*.py")))
            docs_for_api_refs.extend(sorted(examples_dir.glob("*.md")))
        allow_extra = {"__version__"}
        for path in docs_for_api_refs:
            if not path.exists():
                continue
            refs = _extract_doc_api_refs(path)
            unknown = sorted((refs - manifest_set) - allow_extra)
            if unknown:
                errors.append(
                    f"{path.relative_to(REPO)} references non-public API names: {', '.join(unknown)}"
                )

    if errors:
        print("Documentation consistency check FAILED:")
        for err in errors:
            print(f"- {err}")
        return 1

    print("Documentation consistency check PASSED.")
    print(f"- Version: {pyproject_version}")
    print(f"- Active docs checked: {len(active_docs)}")
    if manifest_path.exists():
        print(
            f"- Public API entries: {len(_extract_assigned_str_list(manifest_path, 'PUBLIC_API'))}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
