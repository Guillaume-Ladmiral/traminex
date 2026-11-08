#!/usr/bin/env python3
"""Check bidirectional module, symbol, and signature parity for maintained stubs.

Scope and known limitation (cycle-9 audit N-10): the checker compares module
membership, symbol presence, parameter names/kinds/order, and literal default
values. It deliberately does NOT compare type annotations: stub annotations
refine runtime ones by design (most runtime signatures are partially or wholly
unannotated), so textual annotation equality would flag nearly every stub.
Annotation drift in ``stubs/`` is therefore not machine-policed; reviewers
must check annotations by eye when a signature changes.

Class constructors are compared too (audit S7-02, 2026-09-22): a stub class
declares its constructor either with an explicit ``__init__`` or with
``@dataclass`` (whose constructor is its annotated fields, in order). When the
runtime class defines its own ``__init__`` the stub must declare one, and the
two parameter contracts must agree. For a runtime dataclass every field must
be declared as an annotated stub attribute, and every annotated stub attribute
must exist at runtime. Instance attributes of non-dataclass classes are set in
``__init__`` and cannot be enumerated statically, so they are not compared.
"""

from __future__ import annotations

import argparse
import ast
import dataclasses
import importlib
import inspect
import sys
from collections import defaultdict
from pathlib import Path
from typing import Callable, cast

if __package__:
    from ._checkout_imports import prefer_checkout_imports
else:
    from _checkout_imports import prefer_checkout_imports

REPO = Path(__file__).resolve().parents[1]

ParameterContract = tuple[str, str, str]
FunctionNode = ast.FunctionDef | ast.AsyncFunctionDef


def _module_name(path: Path, root: Path, package: str) -> str:
    rel = path.relative_to(root).with_suffix("")
    parts = list(rel.parts)
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return package + ("" if not parts else "." + ".".join(parts))


def _declared_runtime_modules(runtime_root: Path, package: str) -> dict[str, Path]:
    """Return the maintained public-module set.

    Package modules and modules whose basename does not start with ``_`` are
    declared public. Leading-underscore implementation modules are intentionally
    outside the external stub policy.
    """
    modules: dict[str, Path] = {}
    for path in sorted(runtime_root.rglob("*.py")):
        module = _module_name(path, runtime_root, package)
        if path.name == "__init__.py" or not path.stem.startswith("_"):
            modules[module] = path
    return modules


def _stub_modules(stub_root: Path, package: str) -> dict[str, Path]:
    return {
        _module_name(path, stub_root, package): path for path in sorted(stub_root.rglob("*.pyi"))
    }


def _public_definitions(path: Path) -> tuple[dict[str, list[FunctionNode]], set[str]]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    functions: dict[str, list[FunctionNode]] = defaultdict(list)
    classes: set[str] = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and not node.name.startswith(
            "_"
        ):
            functions[node.name].append(node)
        elif isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
            classes.add(node.name)
    return dict(functions), classes


def _literal_all(path: Path) -> list[str] | None:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in tree.body:
        value: ast.AST | None = None
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "__all__" for target in node.targets
        ):
            value = node.value
        elif (
            isinstance(node, ast.AnnAssign)
            and isinstance(node.target, ast.Name)
            and node.target.id == "__all__"
        ):
            value = node.value
        if isinstance(value, (ast.List, ast.Tuple)) and all(
            isinstance(item, ast.Constant) and isinstance(item.value, str) for item in value.elts
        ):
            return [str(item.value) for item in value.elts if isinstance(item, ast.Constant)]
    return None


def _is_overload(definition: FunctionNode) -> bool:
    return any(
        isinstance(decorator, ast.Name)
        and decorator.id == "overload"
        or isinstance(decorator, ast.Attribute)
        and decorator.attr == "overload"
        for decorator in definition.decorator_list
    )


def _reexports(
    path: Path, module: str, public: set[str] | None, runtime_module: object | None = None
) -> list[tuple[str, str, str]]:
    """``(exported name, source module, original name)`` for each re-export of a stub.

    A re-export is a ``from ... import`` alias whose exported name is in the stub's
    ``__all__``, or, without a literal ``__all__``, the explicit ``from x import y as y``
    form or a renamed ``x as y`` with a public ``y``: such a stub describes ``runtime.y``
    as ``x`` (TX-REP-20261019-004), and a ``y`` the runtime module lacks is reported
    (batch-16 review M-4).
    Imports from outside the package are left out.
    """
    package = module.split(".")[0]
    base = module if path.name == "__init__.pyi" else module.rpartition(".")[0]
    found: list[tuple[str, str, str]] = []
    for node in ast.parse(path.read_text(encoding="utf-8")).body:
        if not isinstance(node, ast.ImportFrom):
            continue
        if node.level:
            parts = base.split(".")
            anchor = ".".join(parts[: len(parts) - node.level + 1])
            source = anchor + ("." + node.module if node.module else "")
        else:
            source = node.module or ""
        if source != package and not source.startswith(package + "."):
            continue
        for alias in node.names:
            exported = alias.asname or alias.name
            if public is not None and exported in public:
                found.append((exported, source, alias.name))
            elif public is None and alias.asname == alias.name:
                found.append((exported, source, alias.name))
            elif (
                public is None
                and alias.asname is not None
                and not alias.asname.startswith("_")
                and runtime_module is not None
            ):
                found.append((exported, source, alias.name))
    return found


def _reexport_errors(
    module: str, path: Path, runtime_module: object, public: set[str] | None
) -> tuple[list[str], int]:
    """A stub re-export must name the object the runtime module exposes (TX-AUD-20261010-009)."""
    errors: list[str] = []
    checked = 0
    for exported, source, original in _reexports(path, module, public, runtime_module):
        if not hasattr(runtime_module, exported):
            if public is None:
                # Without __all__ the name comparison reads definitions only, so a
                # re-export of a name the module lacks is reported here (TX-REP-20261023-006).
                errors.append(
                    f"{module}: '{exported}' re-exports {source}.{original}, but the runtime "
                    f"module has no '{exported}'"
                )
            continue  # with __all__, reported by the name comparison
        try:
            claimed = getattr(importlib.import_module(source), original)
        except (ImportError, AttributeError):
            errors.append(
                f"{module}: '{exported}' re-exports {source}.{original}, which does not exist"
            )
            continue
        checked += 1
        actual = getattr(runtime_module, exported)
        if claimed is not actual:
            where = (
                f"{getattr(actual, '__module__', '?')}.{getattr(actual, '__qualname__', exported)}"
            )
            errors.append(
                f"{module}: '{exported}' re-exports {source}.{original}, but the runtime "
                f"object is {where}"
            )
    return errors, checked


def _stub_parameters(arguments: ast.arguments) -> list[ParameterContract]:
    positional = [*arguments.posonlyargs, *arguments.args]
    first_default = len(positional) - len(arguments.defaults)
    out: list[ParameterContract] = []
    for index, arg in enumerate(positional):
        kind = "POSITIONAL_ONLY" if index < len(arguments.posonlyargs) else "POSITIONAL_OR_KEYWORD"
        status = "default" if index >= first_default else "required"
        out.append((kind, arg.arg, status))
    if arguments.vararg is not None:
        out.append(("VAR_POSITIONAL", arguments.vararg.arg, "variadic"))
    for arg, default in zip(arguments.kwonlyargs, arguments.kw_defaults):
        out.append(("KEYWORD_ONLY", arg.arg, "required" if default is None else "default"))
    if arguments.kwarg is not None:
        out.append(("VAR_KEYWORD", arguments.kwarg.arg, "variadic"))
    return out


_SIMPLE_DEFAULT_TYPES = (type(None), bool, int, float, str, bytes, tuple, frozenset)


def _stub_default_values(arguments: ast.arguments) -> dict[str, object]:
    """Map parameter names to the concrete default values a stub declares.

    Only plain literals participate: ``...`` (the typeshed "unspecified"
    marker) and non-literal expressions are skipped, so this catches a stub
    that states a *wrong* value (audit F-23: ``seqedist`` declared ``True`` /
    ``True`` against runtime ``"No"`` / ``"YujianBo"``) without forcing every
    default to be spelled out.
    """
    values: dict[str, object] = {}
    positional = [*arguments.posonlyargs, *arguments.args]
    tail = positional[len(positional) - len(arguments.defaults) :]
    pairs = list(zip(tail, arguments.defaults))
    pairs.extend(
        (arg, default)
        for arg, default in zip(arguments.kwonlyargs, arguments.kw_defaults)
        if default is not None
    )
    for arg, default in pairs:
        try:
            value = ast.literal_eval(default)
        except (ValueError, TypeError, SyntaxError, MemoryError, RecursionError):
            continue
        if value is Ellipsis:
            continue
        values[arg.arg] = value
    return values


def _runtime_default_values(obj: Callable[..., object]) -> dict[str, object]:
    values: dict[str, object] = {}
    for parameter in inspect.signature(obj).parameters.values():
        if parameter.default is parameter.empty:
            continue
        if isinstance(parameter.default, _SIMPLE_DEFAULT_TYPES):
            values[parameter.name] = parameter.default
    return values


def _defaults_agree(stub_value: object, runtime_value: object) -> bool:
    if isinstance(stub_value, float) and isinstance(runtime_value, float):
        if stub_value != stub_value and runtime_value != runtime_value:  # both NaN
            return True
    if isinstance(stub_value, bool) != isinstance(runtime_value, bool):
        return False
    return bool(stub_value == runtime_value)


def _default_value_errors(
    module: str,
    name: str,
    definitions: list[FunctionNode],
    runtime_obj: Callable[..., object],
) -> list[str]:
    implementations = [definition for definition in definitions if not _is_overload(definition)]
    selected = implementations or definitions
    runtime_values = _runtime_default_values(runtime_obj)
    variant_errors: list[list[str]] = []
    for definition in selected:
        stub_values = _stub_default_values(definition.args)
        mismatches = [
            f"{module}.{name} default mismatch for '{parameter}': "
            f"stub={stub_values[parameter]!r} runtime={runtime_values[parameter]!r}"
            for parameter in sorted(stub_values.keys() & runtime_values.keys())
            if not _defaults_agree(stub_values[parameter], runtime_values[parameter])
        ]
        if not mismatches:
            return []
        variant_errors.append(mismatches)
    return min(variant_errors, key=len) if variant_errors else []


def _runtime_parameters(obj: Callable[..., object]) -> list[ParameterContract]:
    out: list[ParameterContract] = []
    for parameter in inspect.signature(obj).parameters.values():
        if parameter.kind in (parameter.VAR_POSITIONAL, parameter.VAR_KEYWORD):
            status = "variadic"
        else:
            status = "required" if parameter.default is parameter.empty else "default"
        out.append((parameter.kind.name, parameter.name, status))
    return out


def _stub_contracts(definitions: list[FunctionNode]) -> list[list[ParameterContract]]:
    implementations = [definition for definition in definitions if not _is_overload(definition)]
    selected = implementations or definitions
    contracts = [_stub_parameters(definition.args) for definition in selected]
    return [
        contract for index, contract in enumerate(contracts) if contract not in contracts[:index]
    ]


def _stub_class_nodes(path: Path) -> dict[str, ast.ClassDef]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    return {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}


def _decorator_name(decorator: ast.expr) -> str | None:
    target = decorator.func if isinstance(decorator, ast.Call) else decorator
    if isinstance(target, ast.Name):
        return target.id
    if isinstance(target, ast.Attribute):
        return target.attr
    return None


def _dataclass_decorator(node: ast.ClassDef) -> ast.expr | None:
    for decorator in node.decorator_list:
        if _decorator_name(decorator) == "dataclass":
            return decorator
    return None


def _is_classvar(annotation: ast.expr) -> bool:
    target = annotation.value if isinstance(annotation, ast.Subscript) else annotation
    return _decorator_name(target) == "ClassVar"


def _stub_fields(node: ast.ClassDef, classes: dict[str, ast.ClassDef]) -> list[ast.AnnAssign]:
    """Annotated instance attributes of a stub class, in dataclass field order.

    Fields of base classes declared in the same stub module come first; a
    redeclared field keeps its original position, as ``dataclasses`` does.
    """
    merged: dict[str, ast.AnnAssign] = {}
    for base in node.bases:
        if isinstance(base, ast.Name) and base.id in classes and base.id != node.name:
            for inherited in _stub_fields(classes[base.id], classes):
                merged[cast(ast.Name, inherited.target).id] = inherited
    for statement in node.body:
        if (
            isinstance(statement, ast.AnnAssign)
            and isinstance(statement.target, ast.Name)
            and not _is_classvar(statement.annotation)
        ):
            merged[statement.target.id] = statement
    return list(merged.values())


def _stub_members(node: ast.ClassDef, classes: dict[str, ast.ClassDef]) -> set[str]:
    """Every name a stub class declares (attributes, methods, properties), with bases."""
    names = {cast(ast.Name, field.target).id for field in _stub_fields(node, classes)}
    for base in node.bases:
        if isinstance(base, ast.Name) and base.id in classes and base.id != node.name:
            names |= _stub_members(classes[base.id], classes)
    for statement in node.body:
        if isinstance(statement, (ast.FunctionDef, ast.AsyncFunctionDef)):
            names.add(statement.name)
        elif isinstance(statement, ast.Assign):
            names.update(target.id for target in statement.targets if isinstance(target, ast.Name))
    return names


def _stub_constructor(
    node: ast.ClassDef, classes: dict[str, ast.ClassDef]
) -> tuple[list[list[ParameterContract]], dict[str, object]] | None:
    """Return a stub class's constructor contracts and literal defaults, or None.

    An explicit ``__init__`` wins (``@dataclass`` keeps a hand-written one);
    otherwise a ``@dataclass`` stub's constructor is its fields; otherwise the
    first same-module stub base that declares a constructor supplies it.
    """
    inits = [
        statement
        for statement in node.body
        if isinstance(statement, (ast.FunctionDef, ast.AsyncFunctionDef))
        and statement.name == "__init__"
    ]
    if inits:
        contracts = [contract[1:] for contract in _stub_contracts(inits)]
        implementations = [init for init in inits if not _is_overload(init)] or inits
        return contracts, _stub_default_values(implementations[0].args)
    decorator = _dataclass_decorator(node)
    if decorator is not None:
        kw_only = isinstance(decorator, ast.Call) and any(
            keyword.arg == "kw_only"
            and isinstance(keyword.value, ast.Constant)
            and keyword.value.value is True
            for keyword in decorator.keywords
        )
        kind = "KEYWORD_ONLY" if kw_only else "POSITIONAL_OR_KEYWORD"
        contract: list[ParameterContract] = []
        defaults: dict[str, object] = {}
        for field in _stub_fields(node, classes):
            name = cast(ast.Name, field.target).id
            contract.append((kind, name, "required" if field.value is None else "default"))
            if field.value is None:
                continue
            try:
                value = ast.literal_eval(field.value)
            except (ValueError, TypeError, SyntaxError, MemoryError, RecursionError):
                continue
            if value is not Ellipsis:
                defaults[name] = value
        return [contract], defaults
    for base in node.bases:
        if isinstance(base, ast.Name) and base.id in classes and base.id != node.name:
            inherited = _stub_constructor(classes[base.id], classes)
            if inherited is not None:
                return inherited
    return None


def _format_contract(contract: list[ParameterContract]) -> str:
    parts: list[str] = []
    keyword_marker = False
    for kind, name, _status in contract:
        if kind == "VAR_POSITIONAL":
            parts.append(f"*{name}")
            keyword_marker = True
        elif kind == "VAR_KEYWORD":
            parts.append(f"**{name}")
        else:
            if kind == "KEYWORD_ONLY" and not keyword_marker:
                parts.append("*")
                keyword_marker = True
            parts.append(name)
    return "(" + ", ".join(parts) + ")"


def _class_errors(
    module: str,
    name: str,
    node: ast.ClassDef,
    classes: dict[str, ast.ClassDef],
    runtime_cls: type,
) -> tuple[list[str], bool]:
    """Compare one stub class with its runtime class; return (errors, compared)."""
    errors: list[str] = []
    qualified = f"{module}.{name}"
    declared = _stub_constructor(node, classes)
    runtime_defines_init = "__init__" in vars(runtime_cls)
    compared = False
    if declared is None:
        if runtime_defines_init:
            try:
                shown = _format_contract(_runtime_parameters(runtime_cls))
            except (TypeError, ValueError):
                shown = "(...)"
            errors.append(
                f"{qualified}: stub declares no constructor, but the runtime class "
                f"defines __init__{shown}"
            )
    else:
        contracts, stub_defaults = declared
        try:
            actual = _runtime_parameters(runtime_cls)
        except (TypeError, ValueError) as exc:
            errors.append(f"could not inspect {qualified} constructor: {exc}")
        else:
            compared = True
            if actual not in contracts:
                errors.append(
                    f"{qualified} constructor mismatch: stub={contracts} runtime={actual}"
                )
            runtime_defaults = _runtime_default_values(runtime_cls)
            errors.extend(
                f"{qualified}.__init__ default mismatch for '{parameter}': "
                f"stub={stub_defaults[parameter]!r} runtime={runtime_defaults[parameter]!r}"
                for parameter in sorted(stub_defaults.keys() & runtime_defaults.keys())
                if not _defaults_agree(stub_defaults[parameter], runtime_defaults[parameter])
            )
    if dataclasses.is_dataclass(runtime_cls):
        members = _stub_members(node, classes)
        # A private field outside the constructor (``field(init=False)``) is
        # implementation state, not public API, so its stub may omit it.
        runtime_fields = [
            field.name
            for field in dataclasses.fields(runtime_cls)
            if field.init or not field.name.startswith("_")
        ]
        errors.extend(
            f"{qualified}: stub does not declare field '{field}'"
            for field in runtime_fields
            if field not in members
        )
        errors.extend(
            f"{qualified}: stub declares attribute '{attribute}' that the runtime "
            "dataclass does not have"
            for attribute in (
                cast(ast.Name, field.target).id for field in _stub_fields(node, classes)
            )
            if attribute not in runtime_fields and not hasattr(runtime_cls, attribute)
        )
    return errors, compared


def check_package(
    runtime_root: Path, stub_root: Path, package: str
) -> tuple[list[str], dict[str, int]]:
    errors: list[str] = []
    runtime_modules = _declared_runtime_modules(runtime_root, package)
    stub_modules = _stub_modules(stub_root, package)

    for module in sorted(runtime_modules.keys() - stub_modules.keys()):
        errors.append(f"missing stub module {module}")
    for module in sorted(stub_modules.keys() - runtime_modules.keys()):
        errors.append(f"orphan stub module {module}")

    sys.path.insert(0, str(runtime_root.parent))
    # An editable install's meta-path finder outranks sys.path, so without this
    # the source files under --runtime-root were parsed while another tree's
    # modules were imported and compared (audit S7-12).
    prefer_checkout_imports()
    function_count = 0
    class_count = 0
    signature_count = 0
    constructor_count = 0
    reexport_count = 0
    try:
        for module in sorted(runtime_modules.keys() & stub_modules.keys()):
            runtime_functions, runtime_classes = _public_definitions(runtime_modules[module])
            stub_functions, stub_classes = _public_definitions(stub_modules[module])
            stub_class_nodes = _stub_class_nodes(stub_modules[module])

            try:
                runtime_module = importlib.import_module(module)
            except Exception as exc:  # noqa: BLE001
                errors.append(f"failed to import {module}: {exc}")
                continue
            origin = getattr(runtime_module, "__file__", None)
            if origin is None or not Path(origin).resolve().is_relative_to(runtime_root.resolve()):
                errors.append(
                    f"{module} was imported from {origin}, outside --runtime-root "
                    f"{runtime_root}; another installation shadows the tree being checked"
                )
                continue

            runtime_all = _literal_all(runtime_modules[module])
            runtime_names = (
                set(runtime_all)
                if runtime_all is not None
                else set(runtime_functions) | runtime_classes
            )
            stub_all = _literal_all(stub_modules[module])
            stub_names = (
                set(stub_all) if stub_all is not None else set(stub_functions) | stub_classes
            )
            reexport_errors, reexports_checked = _reexport_errors(
                module,
                stub_modules[module],
                runtime_module,
                set(stub_all) if stub_all is not None else None,
            )
            errors.extend(reexport_errors)
            reexport_count += reexports_checked

            for name in sorted(runtime_names):
                runtime_obj = getattr(runtime_module, name, None)
                if inspect.isclass(runtime_obj):
                    class_count += 1
                elif callable(runtime_obj):
                    function_count += 1

            for name in sorted(runtime_names - stub_names):
                runtime_obj = getattr(runtime_module, name, None)
                kind = (
                    "class"
                    if inspect.isclass(runtime_obj)
                    else "function" if callable(runtime_obj) else "symbol"
                )
                errors.append(f"{module}: stub missing public {kind} '{name}'")
            for name in sorted(stub_names - runtime_names):
                if not hasattr(runtime_module, name):
                    kind = (
                        "class"
                        if name in stub_classes
                        else "function" if name in stub_functions else "symbol"
                    )
                    errors.append(f"runtime module {module} missing stub {kind} '{name}'")

            for name in sorted(stub_functions):
                if not hasattr(runtime_module, name):
                    if name in runtime_names:
                        errors.append(f"runtime module {module} missing public function '{name}'")
                    continue
                runtime_obj = getattr(runtime_module, name)
                if not callable(runtime_obj):
                    errors.append(f"runtime attribute {module}.{name} is not callable")
                    continue
                expected = _stub_contracts(stub_functions[name])
                try:
                    actual = _runtime_parameters(runtime_obj)
                except (TypeError, ValueError) as exc:
                    errors.append(f"could not inspect {module}.{name}: {exc}")
                    continue
                signature_count += 1
                if actual not in expected:
                    errors.append(
                        f"{module}.{name} signature mismatch: stub={expected} runtime={actual}"
                    )
                errors.extend(
                    _default_value_errors(module, name, stub_functions[name], runtime_obj)
                )

            for name in sorted(stub_classes):
                if not hasattr(runtime_module, name):
                    if name in runtime_names:
                        errors.append(f"runtime module {module} missing public class '{name}'")
                elif not inspect.isclass(getattr(runtime_module, name)):
                    errors.append(f"runtime attribute {module}.{name} is not a class")
                else:
                    class_errors, compared = _class_errors(
                        module,
                        name,
                        stub_class_nodes[name],
                        stub_class_nodes,
                        getattr(runtime_module, name),
                    )
                    errors.extend(class_errors)
                    constructor_count += int(compared)
    finally:
        sys.path.pop(0)

    totals = {
        "runtime_modules": len(runtime_modules),
        "stub_modules": len(stub_modules),
        "functions": function_count,
        "classes": class_count,
        "signatures": signature_count,
        "constructors": constructor_count,
        "reexports": reexport_count,
    }
    return errors, totals


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runtime-root", type=Path, default=REPO / "traminex")
    parser.add_argument("--stub-root", type=Path, default=REPO / "stubs" / "traminex")
    parser.add_argument("--package", default="traminex")
    args = parser.parse_args(argv)

    if not args.runtime_root.exists():
        parser.error(f"runtime root not found: {args.runtime_root}")
    if not args.stub_root.exists():
        parser.error(f"stub root not found: {args.stub_root}")

    errors, totals = check_package(args.runtime_root, args.stub_root, args.package)
    if errors:
        print("Stub contract check FAILED:")
        for error in errors:
            print(f"- {error}")
    else:
        print("Stub contract check PASSED.")
    print(f"- Runtime modules: {totals['runtime_modules']}")
    print(f"- Stub modules: {totals['stub_modules']}")
    print(f"- Public functions: {totals['functions']}")
    print(f"- Public classes: {totals['classes']}")
    print(f"- Signatures checked: {totals['signatures']}")
    print(f"- Constructors checked: {totals['constructors']}")
    print(f"- Re-exports checked: {totals['reexports']}")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
