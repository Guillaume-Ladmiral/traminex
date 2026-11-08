#!/usr/bin/env python3
"""Synchronize stubs/traminex/__init__.pyi __all__ from traminex/api_manifest.py."""

from __future__ import annotations

import argparse
import ast
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


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


def _extract_assigned_list(path: Path, var_name: str) -> list[str]:
    mod = ast.parse(_read(path), filename=str(path))
    for node in mod.body:
        value = None
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
    raise RuntimeError(f"Could not parse {var_name} from {path}")


def _render_all_block(names: list[str]) -> str:
    lines = ["__all__ = ["]
    for name in names:
        lines.append(f'    "{name}",')
    lines.append("]")
    return "\n".join(lines)


def _replace_all_block(text: str, new_block: str) -> str:
    pattern = re.compile(r"(?ms)^__all__\s*=\s*\[.*?\]\s*")
    if not pattern.search(text):
        raise RuntimeError("Could not find __all__ block in stubs/traminex/__init__.pyi")
    return pattern.sub(new_block + "\n", text, count=1)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="Check only; do not write changes")
    ap.add_argument("--manifest", type=Path, default=REPO / "traminex" / "api_manifest.py")
    ap.add_argument("--stub-init", type=Path, default=REPO / "stubs" / "traminex" / "__init__.pyi")
    args = ap.parse_args()

    manifest_path = args.manifest
    stub_init = args.stub_init
    if not manifest_path.exists():
        raise SystemExit("Missing traminex/api_manifest.py")
    if not stub_init.exists():
        raise SystemExit("Missing stubs/traminex/__init__.pyi")

    names = _extract_assigned_list(manifest_path, "PUBLIC_API")
    current = _extract_assigned_list(stub_init, "__all__")
    if current == names:
        print(f"Stub __all__ matches API manifest exactly ({len(names)} ordered names).")
        return 0

    if args.check:
        print("Stub __all__ does not match API manifest.")
        return 1

    updated = _replace_all_block(_read(stub_init), _render_all_block(names))
    stub_init.write_text(updated, encoding="utf-8")
    print(f"Updated {stub_init} __all__ from API manifest.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
