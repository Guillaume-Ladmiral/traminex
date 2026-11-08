# Type Stubs

This folder contains the maintained `.pyi` type stubs for `traminex`.

The external-stub policy covers every package module and every runtime module
whose basename does not start with an underscore. Within those 66 declared
modules, literal ``__all__`` lists define the public symbol set; modules without
``__all__`` expose their top-level public functions and classes. Re-exported
stub declarations are accepted only when the runtime module actually exposes
the same object: the checker resolves each ``from ... import name as alias``
whose alias is public (listed in the stub's literal ``__all__`` or, in a stub
without one, any alias not starting with an underscore) to its source module and
original name, and requires the runtime module to have that attribute and the
attribute to be that object. Leading-underscore implementation modules are checked by the
canonical whole-package mypy gate, but are intentionally not duplicated here.

Installed distributions use the inline annotations in ``traminex`` and ship a
PEP 561 ``py.typed`` marker. This maintained stub tree is a review and drift
gate; it is not packaged as a second, conflicting ``traminex`` tree.

## Check and synchronize

From the project root on Ubuntu:

```bash
make stub-check
```

To synchronize the root ``__all__`` block with the canonical API manifest and
then check signatures:

```bash
make stubs
```

The bidirectional check reports module, public function/class, and signature
totals; it detects missing/orphan modules, missing/extra symbols, and parameter
name/kind/required-default drift. It also compares class constructors: a stub
class declares its constructor with an explicit `__init__` or with
`@dataclass`, must do so whenever the runtime class defines `__init__`, and
must declare exactly the fields of a runtime dataclass, and it compares every
public re-export with the runtime object (a stub that re-exports `seqtransn` as
`seqtab` fails). `make stub-check` then type-checks the stub
tree itself with `MYPYPATH=stubs mypy -p traminex`, so a stub that re-exports a
name its module does not declare, or that leaves a class abstract, fails the
gate. The sync command only updates
`stubs/traminex/__init__.pyi`; ``--check`` is non-mutating and verifies the exact
ordered 226-name manifest. It does not replace reviewed annotations with
bulk-generated output.

The check prints its current totals (declared modules, public functions and
classes, and the signatures, constructors and re-exports it compared); they are not
repeated here because they change as the API grows (`docs/DOCS_POLICY.md`).
Separately, `make typecheck` runs bare mypy across every package source
module. An isolated installed-wheel smoke verifies that downstream mypy
discovers the inline annotations through the packaged `py.typed` marker.
