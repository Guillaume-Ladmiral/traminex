"""Leaf labels, classification rules, and assignment for dissimilarity trees."""

from __future__ import annotations

import re
from typing import Any, List, Optional, Sequence, Union, cast

import numpy as np
import pandas as pd

from ._rformat import r_as_character, r_number_text


def _label_text(value: Any) -> str:
    """A factor level or profile value as R's ``as.character`` renders it.

    R's split labels are ``levels(factor(pred))`` and ``disstree.assign``
    compares a profile value with them through ``%in%``, which coerces with
    ``as.character`` too -- so a logical reads ``TRUE``/``FALSE`` on both sides.
    """
    text = r_as_character(value)
    return "NA" if text is None else text


def _r_format(x: float, digits: int) -> str:
    """R's ``format(x, digits = digits)`` for one number (``disstree.R:63-70``).

    R keeps fixed notation unless it is wider than scientific notation, so
    123456 prints as ``123456`` where C's ``%.5g`` gives ``1.2346e+05``; the
    printed threshold is what ``disstree_assign`` evaluates.
    """
    return r_number_text(float(x), digits)


def _root(tree: Any) -> Any:
    return tree.root if hasattr(tree, "root") else tree


def _leaf_ids(root: Any) -> np.ndarray:
    """Map the root's zero-based observation indexes to terminal node IDs."""
    node_root = _root(root)
    if not hasattr(node_root, "info"):
        raise TypeError("tree should contain a DissTreeNode root")
    indexes = np.asarray(node_root.info["ind"], dtype=int)
    size = int(indexes.max()) + 1 if indexes.size else 0
    result = np.full(size, -1, dtype=np.int64)
    stack = [node_root]
    while stack:
        node = stack.pop()
        if node.kids is None:
            result[np.asarray(node.info["ind"], dtype=int)] = int(node.id)
        else:
            stack.append(node.kids[1])
            stack.append(node.kids[0])
    return result


def _split_side_labels(split: Any, varname: str, collapse: str) -> tuple[str, str]:
    """Render the left and right predicates of one fitted split."""
    if split.breaks is not None:
        threshold = _r_format(float(split.breaks), 5)
        sides = [f"{varname} <= {threshold}", f"{varname} > {threshold}"]
    else:
        labels = list(split.labels or [])
        allocation = np.asarray(split.index, dtype=int)
        sides = []
        for side in (1, 2):
            values = [
                _label_text(label) for label, group in zip(labels, allocation) if group == side
            ]
            sides.append(f"{varname} [{collapse.join(values)}]")
    if split.na_group is not None:
        side_index = int(split.na_group) - 1
        sides[side_index] = f"{sides[side_index]} with NA"
    return sides[0], sides[1]


def _leaf_label_map(tree: Any, collapse: str) -> dict[int, Optional[str]]:
    """Build path labels for every terminal node with an iterative DFS.

    A root-only tree has no path: R labels its leaf ``NA`` (``None``).
    """
    if not hasattr(tree, "root") or not hasattr(tree, "data"):
        raise TypeError("labelled leaves require a fitted dissimilarity tree")
    mapping: dict[int, Optional[str]] = {}
    stack: list[tuple[Any, tuple[str, ...]]] = [(tree.root, ())]
    columns = list(tree.data.columns)
    while stack:
        node, path = stack.pop()
        if node.kids is None:
            mapping[int(node.id)] = " & ".join(path) if path else None
            continue
        split = node.split
        if split is None:
            raise ValueError("non-terminal tree node has no split")
        varname = str(columns[int(split.varindex)])
        left_label, right_label = _split_side_labels(split, varname, collapse)
        stack.append((node.kids[1], (*path, right_label)))
        stack.append((node.kids[0], (*path, left_label)))
    return mapping


def disstreeleaf(tree: Any, label: bool = False, collapse: str = ", ") -> np.ndarray:
    """Return terminal node IDs or R-compatible root-to-leaf labels."""
    ids = _leaf_ids(tree)
    if not label or not hasattr(tree, "root"):
        return ids
    labels = _leaf_label_map(tree, collapse)
    return np.asarray([labels[int(node_id)] for node_id in ids], dtype=object)


def _leaf_labels_in_r_level_order(tree: Any, collapse: str) -> List[Optional[str]]:
    """Terminal-path labels in R's factor-level order.

    ``disstree.R`` labels each node when its parent is expanded (both kids
    first, then the first kid's subtree), and ``factor(factor(...))`` keeps
    that level order, so the rule order is this insertion order restricted to
    the terminal nodes -- not a lexicographic sort.
    """
    labels = _leaf_label_map(tree, collapse)
    order: list[int] = []
    stack: list[Any] = [tree.root]
    while stack:
        node = stack.pop()
        if node.kids is None:
            continue
        order.append(int(node.kids[0].id))
        order.append(int(node.kids[1].id))
        stack.append(node.kids[1])
        stack.append(node.kids[0])
    ordered = [labels[node_id] for node_id in order if node_id in labels]
    if not ordered:  # single-leaf tree: the root itself is terminal
        ordered = list(labels.values())
    return list(dict.fromkeys(ordered))


def disstree_get_rules(tree: Any, collapse: str = "; ") -> List[Optional[str]]:
    """Return terminal-path labels in R's level order as R condition syntax.

    A root-only tree gives R's one ``NA`` rule, ``[None]``, which assigns no profile.
    """
    rules: list[Optional[str]] = []
    for label in _leaf_labels_in_r_level_order(tree, collapse):
        if label is None:
            rules.append(None)
            continue
        rule = label.replace("[", "%in% c('")
        rule = rule.replace(collapse, "', '")
        rule = rule.replace("]", "')")
        rules.append(rule)
    return rules


# The covariate name is everything before the operator, so R names with dots
# (N.Eastern), backquoted names (`my var`) and the bare non-syntactic names that
# disstree.get.rules writes (my var <= 5) are all read; the threshold of a
# comparison is its last token.
_RULE_IN = re.compile(r"^(.+?)\s+%in%\s+c\((.*)\)$")
_RULE_CMP = re.compile(r"^(.+)\s+(<=|>)\s+(\S+)$")
# R parses each rule (disstree.R:166-167): a level is an R string literal, in single or double quotes, in which a
# backslash escapes the next character. Text that is not such a list of literals (an unescaped quote inside a level,
# or a doubled quote '') is not valid R, and R's disstree.assign stops on it.
_R_STRING = r"'(?:[^'\\]|\\.)*'|\"(?:[^\"\\]|\\.)*\""
_RULE_VALUES = re.compile(rf"^\s*(?:{_R_STRING})(?:\s*,\s*(?:{_R_STRING}))*\s*$", re.S)


# R's string escapes (the parser's StringValue): the single-character ones, octal up to three digits
# (at most \377), \x with one or two hex digits, \u and \U with up to four and eight (or in braces).
# Any other escape, a missing digit or a nul stops R (batch-12 review M-3).
_R_SIMPLE_ESCAPES = {
    "a": "\a",
    "b": "\b",
    "f": "\f",
    "n": "\n",
    "r": "\r",
    "t": "\t",
    "v": "\v",
    "\\": "\\",
    '"': '"',
    "'": "'",
    "`": "`",
    " ": " ",
    "\n": "\n",
}
_R_ESCAPE = re.compile(
    r"\\(?:([0-7]{1,3})|x([0-9A-Fa-f]{0,2})|u(\{[0-9A-Fa-f]*\}?|[0-9A-Fa-f]{0,4})"
    r"|U(\{[0-9A-Fa-f]*\}?|[0-9A-Fa-f]{0,8})|(.))",
    re.S,
)


def _r_unescape(body: str, term: str) -> str:
    """A string literal's text as R's parser reads it; ValueError where R stops."""

    def replace(match: "re.Match[str]") -> str:
        octal, hexadecimal, short, long, other = match.groups()
        if other is not None:
            if other not in _R_SIMPLE_ESCAPES:
                raise ValueError(f"cannot parse rule term: {term!r}")
            return _R_SIMPLE_ESCAPES[other]
        if octal is not None or hexadecimal is not None:
            digits, base = (octal, 8) if octal is not None else (hexadecimal, 16)
            if not digits or int(digits, base) > 0o377:
                raise ValueError(f"cannot parse rule term: {term!r}")
            value = int(digits, base)
            if value >= 0x80:
                # A raw byte: R's string is not valid UTF-8 and equals no text; a lone
                # surrogate (Python's escape for an undecodable byte) equals none either.
                return chr(0xDC00 + value)
        else:
            digits = cast(str, short if short is not None else long)
            braced = digits.startswith("{")
            if braced:
                if (
                    not digits.endswith("}")
                    or len(digits) == 2
                    or len(digits) > (6 if short else 10)
                ):
                    raise ValueError(f"cannot parse rule term: {term!r}")
                digits = digits[1:-1]
            if not digits or int(digits, 16) > 0x10FFFF:
                raise ValueError(f"cannot parse rule term: {term!r}")
            value = int(digits, 16)
        if value == 0:
            raise ValueError(f"cannot parse rule term: {term!r}")  # R: nul character not allowed
        return chr(value)

    return _R_ESCAPE.sub(replace, body)


def _split_rule(rule: str) -> list[str]:
    """Split a rule into its terms at each ``' & '`` outside an R string literal."""
    terms: list[str] = []
    start = index = 0
    quote: Optional[str] = None
    while index < len(rule):
        char = rule[index]
        if quote is not None:
            if char == "\\":
                index += 1
            elif char == quote:
                quote = None
        elif char in "'\"":
            quote = char
        elif rule.startswith(" & ", index):
            terms.append(rule[start:index])
            start = index = index + 3
            continue
        index += 1
    terms.append(rule[start:])
    return terms


def _covariate_name(text: str) -> str:
    """A rule's covariate name, without R's enclosing backquotes."""
    name = text.strip()
    if len(name) > 1 and name.startswith("`") and name.endswith("`"):
        return name[1:-1]
    return name


def _parse_term(term: str):
    """Parse the restricted grammar emitted by :func:`disstree_get_rules`."""
    text = term.strip()
    includes_missing = text.endswith(" with NA")
    if includes_missing:
        text = text[: -len(" with NA")].rstrip()

    membership = _RULE_IN.match(text)
    if membership is not None:
        raw_values = membership.group(2)
        if raw_values.strip() and _RULE_VALUES.match(raw_values) is None:
            raise ValueError(f"cannot parse rule term: {term!r}")
        values = tuple(
            _r_unescape(literal[1:-1], term)
            for literal in re.findall(_R_STRING, raw_values, flags=re.S)
        )
        return ("in", _covariate_name(membership.group(1)), values, includes_missing)

    comparison = _RULE_CMP.match(text)
    if comparison is not None:
        try:
            threshold = float(comparison.group(3))
        except ValueError as exc:
            raise ValueError(f"cannot parse rule term: {term!r}") from exc
        return (
            "cmp",
            _covariate_name(comparison.group(1)),
            comparison.group(2),
            threshold,
            includes_missing,
        )
    raise ValueError(f"cannot parse rule term: {term!r}")


def _is_na_rule(rule: Any) -> bool:
    """Whether a rule is R's ``NA`` (``None`` or a float NaN)."""
    return rule is None or (isinstance(rule, float) and np.isnan(rule))


def _eval_term(term, row: pd.Series) -> bool:
    """Evaluate one parsed predicate against one profile row."""
    kind = term[0]
    value = row[term[1]]
    missing = bool(pd.isna(value))
    if kind == "in":
        _, _, choices, includes_missing = term
        return bool(includes_missing) if missing else _label_text(value) in choices
    _, _, operator, threshold, includes_missing = term
    if missing:
        return bool(includes_missing)
    numeric = float(value)
    return numeric <= threshold if operator == "<=" else numeric > threshold


def disstree_assign(
    rules: Sequence[Optional[str]],
    profile: Union[pd.DataFrame, "pd.Series", dict],
    covar: Optional[Sequence[str]] = None,
) -> np.ndarray:
    """Assign profiles to the first matching terminal rule (1-based).

    .. note:: **Deliberate divergence.** For overlapping (user-supplied)
       rules R's ``disstree.assign`` aborts with a vector-length assignment
       error; traminex keeps the documented first-match semantics. Rules
       produced by :func:`disstree_get_rules` partition the space, so the two
       agree on every fitted tree.

    Covariate names may contain dots (``N.Eastern``) or be backquoted as in R
    (```my var```). A bare non-syntactic name, which R's
    ``disstree.get.rules`` writes for such a covariate (``my var <= 5``), is
    read too; R's own ``disstree.assign`` cannot parse those rules
    (deliberate correction, recorded in the semantic ledger).
    """
    if isinstance(profile, pd.DataFrame):
        frame = profile.copy()
    elif isinstance(profile, pd.Series):
        frame = profile.to_frame().T
    else:
        frame = pd.DataFrame([profile])

    # R's NA rule (a root-only tree's) has no terms and evaluates to NA: it names no covariate and matches no profile.
    parsed_rules = [
        None if _is_na_rule(rule) else tuple(_parse_term(piece) for piece in _split_rule(str(rule)))
        for rule in rules
    ]
    inferred = []
    for parsed in parsed_rules:
        for term in parsed or ():
            if term[1] not in inferred:
                inferred.append(term[1])
    required = list(covar) if covar is not None else inferred
    missing_columns = [name for name in required if name not in frame.columns]
    if missing_columns:
        raise ValueError(f"covariates missing in profile: {missing_columns}")
    if not required and all(parsed is not None for parsed in parsed_rules):
        raise ValueError("list of covariate names covar is NULL")

    result = np.full(len(frame), np.nan, dtype=float)
    for row_number, (_, row) in enumerate(frame.iterrows()):
        for rule_number, parsed in enumerate(parsed_rules, start=1):
            if parsed is not None and all(_eval_term(term, row) for term in parsed):
                result[row_number] = float(rule_number)
                break
    return result
