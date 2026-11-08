"""Shared parser for the R package NAMESPACE files used by audit tooling."""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator
from pathlib import Path


class NamespaceParseError(ValueError):
    """Raised when a NAMESPACE directive cannot be parsed consistently."""


def _strip_r_comments(text: str) -> str:
    """Remove R comments while preserving hashes inside quoted tokens."""
    result: list[str] = []
    quote: str | None = None
    escaped = False
    in_comment = False

    for character in text:
        if in_comment:
            if character in "\r\n":
                result.append(character)
                in_comment = False
            continue

        if quote is not None:
            result.append(character)
            if escaped:
                escaped = False
            elif character == "\\":
                escaped = True
            elif character == quote:
                quote = None
            continue

        if character in {'"', "'", "`"}:
            quote = character
            result.append(character)
        elif character == "#":
            in_comment = True
        else:
            result.append(character)

    if quote is not None:
        raise NamespaceParseError("unterminated quoted token in NAMESPACE")
    return "".join(result)


def _directive_bodies(text: str, directive: str) -> Iterator[str]:
    pattern = re.compile(rf"\b{re.escape(directive)}\s*\(")
    for match in pattern.finditer(text):
        start = match.end()
        depth = 1
        quote: str | None = None
        escaped = False
        for index in range(start, len(text)):
            character = text[index]
            if quote is not None:
                if escaped:
                    escaped = False
                elif character == "\\":
                    escaped = True
                elif character == quote:
                    quote = None
                continue
            if character in {'"', "'", "`"}:
                quote = character
            elif character == "(":
                depth += 1
            elif character == ")":
                depth -= 1
                if depth == 0:
                    yield text[start:index]
                    break
        else:
            raise NamespaceParseError(f"unterminated {directive}(...) directive")


def _split_arguments(content: str) -> list[str]:
    arguments: list[str] = []
    start = 0
    depth = 0
    quote: str | None = None
    escaped = False

    for index, character in enumerate(content):
        if quote is not None:
            if escaped:
                escaped = False
            elif character == "\\":
                escaped = True
            elif character == quote:
                quote = None
            continue
        if character in {'"', "'", "`"}:
            quote = character
        elif character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
        elif character == "," and depth == 0:
            arguments.append(content[start:index].strip())
            start = index + 1
    arguments.append(content[start:].strip())
    return [argument for argument in arguments if argument]


def _canonical_token(token: str) -> str:
    token = token.strip()
    if len(token) >= 2 and token[0] in {'"', "'", "`"} and token[-1] == token[0]:
        token = token[1:-1]
    if not token:
        raise NamespaceParseError("empty name in NAMESPACE directive")
    return token


def _add_name(names: dict[str, str], name: str, kind: str) -> None:
    existing = names.get(name)
    if existing is not None and existing != kind:
        raise NamespaceParseError(f"namespace name is both export and S3 method: {name}")
    names.setdefault(name, kind)


def parse_namespace(path: Path) -> dict[str, str]:
    """Return active export and canonical S3 names in deterministic source order."""
    text = _strip_r_comments(path.read_text(encoding="utf-8", errors="ignore"))
    names: dict[str, str] = {}

    for content in _directive_bodies(text, "export"):
        for token in _split_arguments(content):
            name = _canonical_token(token)
            _add_name(names, name, "export")
            if name.endswith("<-"):
                _add_name(names, name[:-2], "export")

    for content in _directive_bodies(text, "S3method"):
        arguments = _split_arguments(content)
        if len(arguments) < 2:
            raise NamespaceParseError("S3method(...) requires a generic and class")
        generic = _canonical_token(arguments[0])
        class_name = _canonical_token(arguments[1])
        _add_name(names, f"{generic}.{class_name}", "s3")

    return names


#: R objects that traminex promotes to its public API although R defines them
#: without exporting them.  The NAMESPACE file cannot describe these -- it lists
#: only exports and S3 methods -- so a promoted internal was structurally
#: unrepresentable in the decision CSV, whose ``r_kind`` domain was
#: ``export``/``s3``. Seven public analytical surfaces therefore carried no
#: authoritative parity disposition at all (cycle-4b F-004).
#:
#: Each entry maps the R object name to (package, source file relative to
#: ``SOURCES/``). :func:`promoted_internals` verifies both halves of the claim
#: against the vendored tree: the object is *defined* in that file and is *absent*
#: from that package's NAMESPACE.
#:
#: Six further promotions -- ``plot.SAMM``, ``plot.seqimplic``,
#: ``plot.subseqelist``, ``plot.subseqelistchisq``, ``plot.eseq`` and
#: ``plot.seqelist`` -- are S3 methods that the NAMESPACE *does* register, so they
#: are already carried with ``r_kind = s3`` and are deliberately not listed here.
PROMOTED_INTERNALS: dict[str, tuple[str, str]] = {
    "STS_to_SPELL": ("TraMineR", "TraMineR/R/seqformat-STS_to_SPELL.R"),
    "STS_to_SRS": ("TraMineR", "TraMineR/R/seqformat-STS_to_SRS.R"),
    "STS_to_TSE": ("TraMineR", "TraMineR/R/seqformat-STS_to_TSE.R"),
    "entropy": ("TraMineR", "TraMineR/R/entropy.R"),
    "is.subseqelist": ("TraMineR", "TraMineR/R/subseqelist.R"),
    "seqasnum": ("TraMineR", "TraMineR/R/seqasnum.R"),
    "seqerules": ("TraMineR", "TraMineR/R/seqerules.R"),
}


def _defines_object(source: Path, name: str) -> bool:
    """Whether an R source file assigns ``name`` as a function at top level."""
    pattern = re.compile(
        rf"""^["']?{re.escape(name)}["']?\s*(?:<-|=)\s*function""",
        re.MULTILINE,
    )
    return bool(pattern.search(source.read_text(encoding="utf-8", errors="ignore")))


def promoted_internals(sources_root: Path) -> dict[str, str]:
    """Return ``{r_name: package}`` for the verified promoted internals.

    Raises :class:`NamespaceParseError` when an entry no longer holds -- the
    object is not defined where recorded, or the package has started exporting
    it, in which case it belongs in the ordinary ``export`` vocabulary instead.
    """
    verified: dict[str, str] = {}
    for name, (package, relative) in PROMOTED_INTERNALS.items():
        source = sources_root / relative
        if not source.is_file():
            raise NamespaceParseError(f"promoted internal {name}: missing source {relative}")
        if not _defines_object(source, name):
            raise NamespaceParseError(f"promoted internal {name}: not defined in {relative}")
        exported = parse_namespace(sources_root / package / "NAMESPACE")
        if name in exported:
            raise NamespaceParseError(
                f"promoted internal {name}: {package} now exports it; "
                "move the row to r_kind=export"
            )
        verified[name] = package
    return verified


#: Top-level ``name <- function`` / ``name = function`` assignment in an R source.
_R_FUNCTION_DEFINITION = re.compile(
    r"""^(["']?)([A-Za-z.][A-Za-z0-9._]*)\1\s*(?:<-|=)\s*function""",
    re.MULTILINE,
)


def vendored_object_definitions(sources_root: Path) -> dict[str, set[str]]:
    """Return ``{r_name: {package, ...}}`` for every top-level function in the sources.

    Derived from the vendored tree rather than from a live R session, so the
    ledger validator keeps working in an environment with no R -- which is the
    environment ``make generated-check`` runs in. It is deliberately a *superset*
    of what a live ``ls(asNamespace(...))`` reports for our purpose: any name it
    misses could not be matched against a Python public name anyway, and callers
    narrow it by NAMESPACE membership.
    """
    definitions: dict[str, set[str]] = {}
    for package in ("TraMineR", "TraMineRextras"):
        directory = sources_root / package / "R"
        if not directory.is_dir():
            raise NamespaceParseError(f"missing vendored R sources: {directory}")
        for source in sorted(directory.glob("*.R")):
            text = source.read_text(encoding="utf-8", errors="ignore")
            for match in _R_FUNCTION_DEFINITION.finditer(text):
                definitions.setdefault(match.group(2), set()).add(package)
    if not definitions:
        raise NamespaceParseError("no R function definitions found in the vendored sources")
    return definitions


def unrecorded_internal_promotions(
    public_api: Iterable[str],
    sources_root: Path,
) -> list[tuple[str, str, str]]:
    """Python public names that shadow an unexported R internal without a record.

    The ledger validator enforces the R -> Python direction only: every R
    namespace name must have exactly one row. It cannot see the reverse, so a
    traminex public name with no row is invisible -- and 64 of the 226 legitimately
    have none, being traminex-only extensions. A *future* promotion of an
    unexported R internal would land silently in that same bucket instead of being
    forced into :data:`PROMOTED_INTERNALS` (cycle-6 section 8.1).

    A public name is checked under both its own spelling and its dot-for-underscore
    R spelling. A match is accounted for when the name appears in either NAMESPACE
    -- as an export or as a registered S3 method -- or when it is already a recorded
    promotion. Anything else is returned as ``(python_name, r_name, packages)``.
    """
    definitions = vendored_object_definitions(sources_root)
    namespace: dict[str, tuple[str, str]] = {}
    for package in ("TraMineR", "TraMineRextras"):
        for name, kind in parse_namespace(sources_root / package / "NAMESPACE").items():
            namespace.setdefault(name, (package, kind))

    unrecorded: list[tuple[str, str, str]] = []
    for python_name in public_api:
        for candidate in (python_name, python_name.replace("_", ".")):
            if candidate not in definitions:
                continue
            if candidate in namespace or candidate in PROMOTED_INTERNALS:
                break
            unrecorded.append((python_name, candidate, ", ".join(sorted(definitions[candidate]))))
            break
    return unrecorded
