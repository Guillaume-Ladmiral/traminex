#!/usr/bin/env python3
"""Generate the TraMineR/TraMineRextras API audit on Linux."""

import argparse
import ast
import csv
import datetime
import difflib
import importlib
import json
import re
import sys
from pathlib import Path

if __package__:
    from ._r_namespace import parse_namespace
else:
    from _r_namespace import parse_namespace

# Every input and the default output are anchored to this repository, so the
# script (even --help) works from any working directory (audit S7-18).
REPO = Path(__file__).resolve().parents[1]

# Helpers


# Load R NAMESPACE
traminer_ns = REPO / "SOURCES" / "TraMineR" / "NAMESPACE"
extras_ns = REPO / "SOURCES" / "TraMineRextras" / "NAMESPACE"
traminer_names = parse_namespace(traminer_ns)
extras_names = parse_namespace(extras_ns)
traminer_exports = {name for name, kind in traminer_names.items() if kind == "export"}
traminer_s3 = {name for name, kind in traminer_names.items() if kind == "s3"}
extras_exports = {name for name, kind in extras_names.items() if kind == "export"}
extras_s3 = {name for name, kind in extras_names.items() if kind == "s3"}

r_exports = traminer_exports | extras_exports
r_s3 = traminer_s3 | extras_s3
r_all = sorted(r_exports | r_s3)

# Load Python public API from canonical manifest (fallback: __all__)
manifest_path = REPO / "traminex" / "api_manifest.py"
init_path = REPO / "traminex" / "__init__.py"
manifest_src = manifest_path.read_text(encoding="utf-8") if manifest_path.exists() else ""
init_src = init_path.read_text(encoding="utf-8")
manifest_mod = ast.parse(manifest_src, filename=str(manifest_path)) if manifest_src else None
init_mod = ast.parse(init_src, filename=str(init_path))


def _extract_literal_str_list(node: ast.AST):
    if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
        out = []
        for elt in node.elts:
            if isinstance(elt, ast.Constant) and isinstance(elt.value, str):
                out.append(elt.value)
            else:
                return None
        return out
    return None


def _extract_all_names(node: ast.AST):
    direct = _extract_literal_str_list(node)
    if direct is not None:
        return direct
    if isinstance(node, ast.Call):
        # list([...]) / tuple([...]) wrappers
        if (
            isinstance(node.func, ast.Name)
            and node.func.id in {"list", "tuple"}
            and len(node.args) == 1
        ):
            inner = _extract_all_names(node.args[0])
            if inner is not None:
                return list(inner)
        # list(dict.fromkeys([...])) pattern (used in traminex __all__)
        if (
            isinstance(node.func, ast.Name)
            and node.func.id == "list"
            and len(node.args) == 1
            and isinstance(node.args[0], ast.Call)
        ):
            inner_call = node.args[0]
            if (
                isinstance(inner_call.func, ast.Attribute)
                and inner_call.func.attr == "fromkeys"
                and isinstance(inner_call.func.value, ast.Name)
                and inner_call.func.value.id == "dict"
                and len(inner_call.args) >= 1
            ):
                keys = _extract_all_names(inner_call.args[0])
                if keys is not None:
                    return list(dict.fromkeys(keys))
    return None


def _collect_named_assignments(module: ast.Module, var_name: str):
    names = []
    for node in module.body:
        value = None
        if isinstance(node, ast.Assign):
            if any(isinstance(t, ast.Name) and t.id == var_name for t in node.targets):
                value = node.value
        elif isinstance(node, ast.AnnAssign):
            if isinstance(node.target, ast.Name) and node.target.id == var_name:
                value = node.value
        if value is not None:
            parsed = _extract_all_names(value)
            if parsed is not None:
                names.extend([str(x) for x in parsed])
    return names


py_api = []
if manifest_mod is not None:
    py_api = _collect_named_assignments(manifest_mod, "PUBLIC_API")
if not py_api:
    py_api = _collect_named_assignments(init_mod, "__all__")
fallback_error = None
if not py_api:
    try:
        sys.path.insert(0, str(REPO))
        tx = importlib.import_module("traminex")
        names = getattr(tx, "__all__", [])
        if isinstance(names, (list, tuple, set)):
            py_api.extend([str(x) for x in names])
    except Exception as exc:
        fallback_error = f"{type(exc).__name__}: {exc}"
if not py_api:
    msg = "Failed to resolve traminex public API (__all__) from AST and import fallback."
    if fallback_error:
        msg += f" Import fallback error: {fallback_error}"
    raise RuntimeError(msg)
py_api_set = set(py_api)
py_api_list = sorted(py_api_set)

# Mapping heuristics


def sanitize_r_name(name: str) -> str:
    return name.replace(".", "_").replace("<-", "")


known_map = {
    "plot.stslist": ["seqplot", "sequence_index_plot"],
    "plot.stslist.freq": ["seqfplot"],
    "plot.stslist.meant": ["seqmtplot"],
    "plot.stslist.modst": ["seqmsplot"],
    "plot.stslist.rep": ["seqrplot"],
    "plot.stslist.statd": ["seqHtplot", "seqdHplot"],
    "plot.seqrf": ["seqrfplot"],
}

exact = []
exact_by_sanitized = []
suggested = []
missing = []

for rname in r_all:
    if rname in py_api_set:
        exact.append((rname, rname))
        continue
    r_san = sanitize_r_name(rname)
    if r_san in py_api_set:
        exact_by_sanitized.append((rname, r_san))
        continue
    suggestions = []
    if rname in known_map:
        for tgt in known_map[rname]:
            if tgt in py_api_set:
                suggestions.append((tgt, 1.0))
    if not suggestions:
        scores = [
            (p, difflib.SequenceMatcher(None, r_san.lower(), p.lower()).ratio())
            for p in py_api_list
        ]
        scores.sort(key=lambda t: t[1], reverse=True)
        for p, s in scores[:5]:
            if s >= 0.60:
                suggestions.append((p, s))
    if suggestions:
        suggested.append((rname, suggestions))
    else:
        missing.append(rname)


# Versions
def read_desc_version(desc_path: Path):
    if not desc_path.exists():
        return "unknown"
    txt = desc_path.read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"^Version:\s*([^\r\n]+)", txt, flags=re.M)
    return m.group(1).strip() if m else "unknown"


def read_traminex_version(pyproject_path: Path):
    try:
        txt = pyproject_path.read_text(encoding="utf-8")
        m = re.search(r'^version\s*=\s*["\']([^"\']+)["\']', txt, flags=re.M)
        return m.group(1) if m else "unknown"
    except Exception:
        return "unknown"


traminer_ver = read_desc_version(REPO / "SOURCES" / "TraMineR" / "DESCRIPTION")
extras_ver = read_desc_version(REPO / "SOURCES" / "TraMineRextras" / "DESCRIPTION")
traminex_ver = read_traminex_version(REPO / "pyproject.toml")

summary = {
    "TraMineR_exports": len(traminer_exports),
    "TraMineR_s3": len(traminer_s3),
    "TraMineR_total": len(traminer_exports | traminer_s3),
    "TraMineRextras_exports": len(extras_exports),
    "TraMineRextras_s3": len(extras_s3),
    "TraMineRextras_total": len(extras_exports | extras_s3),
    "R_total_unique": len(set(r_all)),
    "Python_public_api": len(py_api_set),
    "exact_matches": len(exact),
    "sanitized_name_matches": len(exact_by_sanitized),
    "with_suggestions": len(suggested),
    "missing": len(missing),
}

# Outputs
parser = argparse.ArgumentParser(description="Generate the traminex API audit.")
parser.add_argument("--output-dir", type=Path, default=REPO / "reports" / "audit")
parser.add_argument(
    "--generated-at",
    help="fixed timezone-aware ISO-8601 timestamp for byte-reproducible outputs",
)
args = parser.parse_args()
out_dir = args.output_dir
out_dir.mkdir(parents=True, exist_ok=True)
if args.generated_at:
    try:
        generated_at = datetime.datetime.fromisoformat(args.generated_at.replace("Z", "+00:00"))
    except ValueError as exc:
        parser.error(f"invalid --generated-at timestamp: {exc}")
    if generated_at.tzinfo is None:
        parser.error("--generated-at must include a timezone offset")
    generated_at = generated_at.astimezone(datetime.timezone.utc)
else:
    generated_at = datetime.datetime.now(datetime.timezone.utc)

# Markdown
now = generated_at.strftime("%Y-%m-%d %H:%M:%SZ")
md_lines = []
md_lines.append("# traminex vs TraMineR/TraMineRextras API audit")
md_lines.append("")
md_lines.append(f"Generated on: {now} (UTC)")
md_lines.append("")
md_lines.append("## Versions")
md_lines.append(f"- TraMineR: {traminer_ver}")
md_lines.append(f"- TraMineRextras: {extras_ver}")
md_lines.append(f"- traminex: {traminex_ver}")
md_lines.append("")
md_lines.append("## Summary")
for k in [
    "TraMineR_exports",
    "TraMineR_s3",
    "TraMineR_total",
    "TraMineRextras_exports",
    "TraMineRextras_s3",
    "TraMineRextras_total",
    "R_total_unique",
    "Python_public_api",
    "exact_matches",
    "sanitized_name_matches",
    "with_suggestions",
    "missing",
]:
    md_lines.append(f"- {k}: {summary[k]}")
md_lines.append("")
md_lines.append("## Notable equivalences (manually curated)")
for r, tgts in known_map.items():
    tgt_list = ", ".join(f"`{t}`" for t in tgts)
    md_lines.append(f"- `{r}` -> {tgt_list}")
md_lines.append("")
md_lines.append("## Exact matches")
for r, p in sorted(exact):
    md_lines.append(f"- `{r}` -> `{p}`")
if not exact:
    md_lines.append("(none)")
md_lines.append("")
md_lines.append("## Sanitized-name matches")
for r, p in sorted(exact_by_sanitized):
    md_lines.append(f"- `{r}` -> `{p}`")
if not exact_by_sanitized:
    md_lines.append("(none)")
md_lines.append("")
md_lines.append("## Suggested mappings")
for rname, suggs in sorted(suggested, key=lambda t: t[0].lower()):
    s = ", ".join([f"`{p}` ({score:.2f})" for p, score in suggs]) if suggs else "(none)"
    md_lines.append(f"- `{rname}` -> {s}")
if not suggested:
    md_lines.append("(none)")
md_lines.append("")
md_lines.append("## Missing (no suggestion >= 0.60)")
for r in sorted(missing):
    md_lines.append(f"- `{r}`")
if not missing:
    md_lines.append("(none)")
md_lines.append("")
md_lines.append("## Appendix: Python public API (manifest / __all__ fallback)\n")
md_lines.append(", ".join(f"`{n}`" for n in sorted(py_api_set)))
md_lines.append("")

(out_dir / "traminex_api_audit.md").write_text("\n".join(md_lines), encoding="utf-8", newline="\n")

# JSON
json_obj = {
    "generated_at_utc": generated_at.isoformat(),
    "versions": {
        "TraMineR": traminer_ver,
        "TraMineRextras": extras_ver,
        "traminex": traminex_ver,
    },
    "summary": summary,
    "exact_matches": [{"r": r, "py": p} for r, p in exact],
    "sanitized_name_matches": [{"r": r, "py": p} for r, p in exact_by_sanitized],
    "suggested_mappings": [
        {"r": r, "suggestions": [{"py": p, "score": round(float(s), 4)} for p, s in suggs]}
        for r, suggs in suggested
    ],
    "missing": missing,
    "python_public_api": sorted(py_api_set),
}
(out_dir / "traminex_api_audit.json").write_text(
    json.dumps(json_obj, indent=2, ensure_ascii=False) + "\n",
    encoding="utf-8",
    newline="\n",
)

# CSV
with (out_dir / "traminex_api_audit.csv").open("w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["r_name", "match_type", "py_name", "score"])
    for r, p in exact:
        w.writerow([r, "exact", p, "1.0"])
    for r, p in exact_by_sanitized:
        w.writerow([r, "sanitized", p, "1.0"])
    for r, suggs in suggested:
        for p, s in suggs:
            w.writerow([r, "suggested", p, f"{float(s):.4f}"])
    for r in missing:
        w.writerow([r, "missing", "", ""])

print("Audit complete: wrote MD, JSON, CSV to", out_dir)
print(
    "Summary: R_total_unique=",
    len(r_all),
    "Python_public_api=",
    len(py_api_set),
    "exact_matches=",
    len(exact),
)
