#!/usr/bin/env python3
"""Validate the syntax and minimum schema of the maintained GitHub workflows."""

from __future__ import annotations

import argparse
import re
import shlex
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKFLOWS = (
    ROOT / ".github" / "workflows" / "ci.yml",
    ROOT / ".github" / "workflows" / "wheels.yml",
    ROOT / ".github" / "workflows" / "render-docs.yml",
)


class _DuplicateKeyError(yaml.YAMLError):
    pass


class _UniqueKeyLoader(yaml.BaseLoader):
    pass


def _construct_unique_mapping(
    loader: _UniqueKeyLoader, node: yaml.MappingNode, deep: bool = False
) -> dict[str, Any]:
    mapping: dict[str, Any] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in mapping:
            line = key_node.start_mark.line + 1
            raise _DuplicateKeyError(f"line {line}: duplicate key {key!r}")
        mapping[key] = loader.construct_object(value_node, deep=deep)
    return mapping


_UniqueKeyLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
    _construct_unique_mapping,
)


def validate_workflow(path: Path) -> list[str]:
    """Return deterministic validation errors for one workflow file."""

    try:
        data = yaml.load(path.read_text(encoding="utf-8"), Loader=_UniqueKeyLoader)
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        return [f"{path}: {exc}"]

    errors: list[str] = []
    if not isinstance(data, dict):
        return [f"{path}: workflow root must be a mapping"]
    if not isinstance(data.get("name"), str) or not data["name"].strip():
        errors.append(f"{path}: name must be a non-empty scalar")
    if "on" not in data:
        errors.append(f"{path}: missing on trigger")
    jobs = data.get("jobs")
    if not isinstance(jobs, dict) or not jobs:
        errors.append(f"{path}: jobs must be a non-empty mapping")
        return errors
    for job_name, job in jobs.items():
        if not isinstance(job, dict):
            errors.append(f"{path}: jobs.{job_name} must be a mapping")
            continue
        if "runs-on" not in job:
            errors.append(f"{path}: jobs.{job_name}.runs-on is required")
        steps = job.get("steps")
        if not isinstance(steps, list) or not steps:
            errors.append(f"{path}: jobs.{job_name}.steps must be a non-empty sequence")
    return errors


SUPPORTED_PYTHONS = ("3.11", "3.12")
_SHELL_SEPARATOR = re.compile(r"&&|\|\||[;|\n]")
#: A backslash before a newline continues the command on the next line.
_CONTINUATION = re.compile(r"\\\n")
_MAKE_TARGET_LINE = re.compile(r"^([A-Za-z0-9][A-Za-z0-9_.-]*)\s*:(?!=)")


def _makefile_targets(makefile: Path) -> set[str]:
    """Return the explicit rule targets defined in ``makefile``."""
    try:
        text = makefile.read_text(encoding="utf-8")
    except (OSError, UnicodeError):
        return set()
    return {
        match.group(1) for line in text.splitlines() if (match := _MAKE_TARGET_LINE.match(line))
    }


def _joined(script: str) -> str:
    """Join backslash-continued lines, as the shell reads them (TX-AUD-20261019-021)."""
    return _CONTINUATION.sub("", script)


def _make_invocations(script: str) -> list[str]:
    """Return the targets named by ``make`` commands in a workflow ``run`` script."""
    targets: list[str] = []
    for command in _SHELL_SEPARATOR.split(_joined(script)):
        words = command.split()
        if not words or words[0] != "make":
            continue
        targets.extend(word for word in words[1:] if not word.startswith("-") and "=" not in word)
    return targets


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else [value]


def _newest_dependency_errors(path: Path, jobs: dict[str, Any]) -> list[str]:
    """ci.yml type-checks with the newest NumPy, pandas and SciPy the floors permit.

    The dependency floors have no upper bounds, so a release can meet a NumPy
    whose stubs the type-check gate has never seen (audit TX-AUD-20261010-010).
    mypy targets the interpreter it runs under, so the job runs on Python 3.11
    and 3.12: a 3.12-only run would pass standard-library calls that 3.11 lacks.
    """
    job = jobs.get("typecheck-newest")
    if not isinstance(job, dict):
        return [
            f"{path}: job typecheck-newest (make typecheck with the newest permitted "
            "numpy, pandas and scipy) is missing"
        ]
    steps = job.get("steps")
    scripts = [
        step["run"]
        for step in (steps if isinstance(steps, list) else [])
        if isinstance(step, dict) and isinstance(step.get("run"), str)
    ]
    text = "\n".join(scripts)
    errors: list[str] = []
    if "pip install --upgrade numpy pandas scipy" not in text:
        errors.append(f"{path}: jobs.typecheck-newest must upgrade numpy, pandas and scipy")
    if "typecheck" not in _make_invocations(text):
        errors.append(f"{path}: jobs.typecheck-newest must run make typecheck")
    strategy = job.get("strategy")
    matrix = strategy.get("matrix") if isinstance(strategy, dict) else None
    versions = matrix.get("python-version") if isinstance(matrix, dict) else None
    rows = {str(version) for version in _as_list(versions)}
    missing = [version for version in SUPPORTED_PYTHONS if version not in rows]
    if missing:
        errors.append(
            f"{path}: jobs.typecheck-newest must run on Python "
            f"{', '.join(SUPPORTED_PYTHONS)} (matrix python-version lacks {', '.join(missing)})"
        )
    setup = [
        step.get("with")
        for step in (steps if isinstance(steps, list) else [])
        if isinstance(step, dict) and str(step.get("uses", "")).startswith("actions/setup-python")
    ]
    if not any(
        isinstance(inputs, dict)
        and str(inputs.get("python-version", "")).replace(" ", "") == "${{matrix.python-version}}"
        for inputs in setup
    ):
        errors.append(
            f"{path}: jobs.typecheck-newest must set up python-version: "
            "${{ matrix.python-version }}"
        )
    return errors


_GATE_COMMAND = re.compile(r"(?:^|[\s;&|(])(?:make|pytest|\S*run_pybind11_tests\.sh)(?:\s|$)")
_TEST_COMMAND = re.compile(r"(?:^|[\s;&|(])(?:pytest|\S*run_pybind11_tests\.sh)(?:\s|$)")
#: After a gate command: any ``||`` list, or ``; true`` / ``; :`` / ``; exit 0``.
_SWALLOWED_FAILURE = re.compile(r"\|\||;\s*(?:true|:|exit\s+0)\b")
#: A single ``|``: without pipefail the pipeline's status is its last command's.
_PIPE = re.compile(r"(?<!\|)\|(?!\|)")
#: ``-k``/``-m`` also with the value attached (``-mslow``); ``--deselect``, ``--ignore``,
#: ``--ignore-glob``.
_TEST_SELECTION = re.compile(r"(?:^|\s)(?:-[km]|--(?:deselect|ignore|ignore-glob)(?:[\s=]|$))")
#: Options with which pytest runs no test (``--co``/``--collect-only`` are checked on the line).
_NO_RUN = re.compile(
    r"--(?:setup-only|setup-plan|fixtures|fixtures-per-test|markers|version|help)|-[Vh]"
)
#: Options that rewrite pytest's configuration (addopts, testpaths), also with the value attached.
_RECONFIGURE = re.compile(r"-[oc]|--(?:override-ini|config-file)(?:=|$)")
_REDIRECTION = re.compile(r"\d*[<>]")
_ERREXIT_OFF = re.compile(r"(?:^|[\s;&(])set\s+(?:\+e|\+o\s+errexit)\b")
#: The only condition a gate step may carry: the test job's matrix backend.
_BACKEND_CONDITION = re.compile(
    r"\s*(?:\$\{\{\s*)?matrix\.backend\s*==\s*'(?:python|pybind11)'\s*(?:\}\})?\s*"
)
#: Shells GitHub runs with ``-e``; ``bash`` (given explicitly) also sets pipefail.
_ERREXIT_SHELLS = {"bash", "sh"}
_DEFAULT_BRANCHES = {"main", "master"}


def _trigger_errors(path: Path, data: dict[str, Any]) -> list[str]:
    """ci.yml runs its gates on every push and pull request (audit TX-AUD-20261010-057)."""
    if path.name != "ci.yml":
        return []
    triggers = data.get("on")
    names = (
        set(triggers)
        if isinstance(triggers, dict)
        else set(_as_list(triggers)) if triggers is not None else set()
    )
    missing = sorted({"push", "pull_request"} - names)
    if missing:
        return [f"{path}: ci.yml must run on push and pull_request (missing {', '.join(missing)})"]
    errors: list[str] = []
    for event in ("push", "pull_request"):
        settings = triggers.get(event) if isinstance(triggers, dict) else None
        if not isinstance(settings, dict):
            continue
        branches = settings.get("branches")
        if branches is not None and not {str(b) for b in _as_list(branches)} & _DEFAULT_BRANCHES:
            errors.append(
                f"{path}: on.{event}.branches must keep the default branch (main or master)"
            )
        ignored = {str(b) for b in _as_list(settings.get("branches-ignore"))}
        if ignored & _DEFAULT_BRANCHES:
            errors.append(
                f"{path}: on.{event}.branches-ignore may not drop the default branch "
                f"({', '.join(sorted(ignored & _DEFAULT_BRANCHES))})"
            )
    return errors


def _tolerates_failure(settings: dict[str, Any]) -> bool:
    """``continue-on-error`` set to anything but false (an expression may be true)."""
    if "continue-on-error" not in settings:
        return False
    value = settings["continue-on-error"]
    return not (value is False or str(value).strip().lower() == "false")


def _default_shell(settings: Any) -> Any:
    defaults = settings.get("defaults") if isinstance(settings, dict) else None
    run = defaults.get("run") if isinstance(defaults, dict) else None
    return run.get("shell") if isinstance(run, dict) else None


def _test_arguments(command: str) -> list[str] | None:
    """The arguments of a ``pytest`` or ``run_pybind11_tests.sh`` command, else ``None``.

    The test command is the first word naming pytest or the script, whatever runs it
    (``python -m``, ``$PYTHON -m``, ``coverage run -m``, ``env VAR=value``; batch-8
    review I-3); a word ``install`` before it marks a package installation.
    """
    try:
        words = shlex.split(command)
    except ValueError:
        words = command.split()
    for index, word in enumerate(words):
        if word == "install":
            return None
        if word.rsplit("/", 1)[-1] in ("pytest", "run_pybind11_tests.sh"):
            return words[index + 1 :]
    return None


def _mentions(value: Any, text: str) -> bool:
    """Whether ``text`` occurs in a key or string anywhere in a parsed workflow."""
    if isinstance(value, str):
        return text in value
    if isinstance(value, dict):
        return any(_mentions(key, text) or _mentions(item, text) for key, item in value.items())
    if isinstance(value, list):
        return any(_mentions(item, text) for item in value)
    return False


def _test_path(arguments: list[str]) -> str | None:
    """The first argument of a test command that is not an option: a test path.

    A detached option value (``--option value``) reads the same way, so a gate writes
    option values attached (``--option=value``). Reading stops at a redirection and at a
    selection option, which the line check reports.
    """
    for word in arguments:
        if _REDIRECTION.match(word) or _TEST_SELECTION.match(word):
            return None
        if not word.startswith("-"):
            return word
    return None


def _neutraliser_errors(path: Path, data: dict[str, Any]) -> list[str]:
    """A gate step may not swallow its failure or skip running its tests.

    Audit TX-AUD-20261010-057 rejected ``|| true``, ``--collect-only`` and
    ``continue-on-error: true``; batch-5 review M-3 found that short list let
    ``|| echo``, ``; true``, ``set +e``, ``continue-on-error: ${{ true }}``,
    ``if: false`` and ``-k`` through. A gate step (one whose script runs make,
    pytest or run_pybind11_tests.sh) now runs under ``-e`` with no ``||`` or
    ``; true`` after the gate, a pipe only with pipefail, no test selection,
    no tolerated failure and no condition other than the matrix backend; and
    the job holding it carries no condition. Continued lines are read joined,
    and selection includes ``-m``, ``--ignore``, ``--ignore-glob`` and test
    paths (TX-AUD-20261019-021: a ``|| true`` on a continuation line and those
    selections passed). The test command is recognised whatever runs it, and
    PYTEST_ADDOPTS, options with which pytest runs no test and options that
    rewrite its configuration are refused (batch-8 review I-3). Path filters and
    job dependencies are not checked here.
    """
    errors: list[str] = []
    if _mentions(data, "PYTEST_ADDOPTS"):
        errors.append(
            f"{path}: PYTEST_ADDOPTS adds options to every pytest run, so a gate's tests "
            "cannot be read from its command"
        )
    jobs = data.get("jobs")
    for job_name, job in (jobs.items() if isinstance(jobs, dict) else []):
        if not isinstance(job, dict):
            continue
        shell_default = _default_shell(job) or _default_shell(data)
        steps = job.get("steps")
        job_has_gate = False
        for number, step in enumerate(steps if isinstance(steps, list) else [], start=1):
            if not isinstance(step, dict) or not isinstance(step.get("run"), str):
                continue
            where = f"{path}: jobs.{job_name}.steps[{number}]"
            script = _joined(step["run"])
            shell = step.get("shell", shell_default)
            pipefail = shell == "bash" or "pipefail" in script
            gate = False
            for line in script.splitlines():
                match = _GATE_COMMAND.search(line)
                if not match:
                    continue
                gate = True
                tail = line[match.start() :]
                if _SWALLOWED_FAILURE.search(tail):
                    errors.append(f"{where}: gate command cannot fail ({line.strip()!r})")
                if _PIPE.search(tail) and not pipefail:
                    errors.append(
                        f"{where}: gate command feeds a pipe without set -o pipefail, so the "
                        f"pipe's last command decides the status ({line.strip()!r})"
                    )
                if "--collect-only" in line or re.search(r"(?:^|\s)--co(?:\s|$)", line):
                    errors.append(
                        f"{where}: gate command collects tests without running them ({line.strip()!r})"
                    )
                test = _TEST_COMMAND.search(line)
                # Only the test command's own arguments: a later ``python -m`` is not -m.
                if test and _TEST_SELECTION.search(_SHELL_SEPARATOR.split(line[test.end() :])[0]):
                    errors.append(
                        f"{where}: gate command selects a subset of the tests ({line.strip()!r})"
                    )
                for command in _SHELL_SEPARATOR.split(line):
                    arguments = _test_arguments(command)
                    if arguments is None:
                        continue
                    if any(_NO_RUN.fullmatch(word) for word in arguments):
                        errors.append(f"{where}: gate command runs no test ({line.strip()!r})")
                    if any(_RECONFIGURE.match(word) for word in arguments):
                        errors.append(
                            f"{where}: gate command overrides pytest's configuration, which "
                            f"can select its tests ({line.strip()!r})"
                        )
                    selected = _test_path(arguments)
                    if selected is not None:
                        errors.append(
                            f"{where}: gate command selects a subset of the tests (test path "
                            f"{selected!r}; an option's value is written --option=value) "
                            f"({line.strip()!r})"
                        )
            if not gate:
                continue
            job_has_gate = True
            if _ERREXIT_OFF.search(script):
                errors.append(f"{where}: a gate step may not run set +e")
            if shell is not None and str(shell).strip() not in _ERREXIT_SHELLS:
                errors.append(
                    f"{where}: a gate step must run under bash or sh, which stop on the first "
                    f"failure (shell {shell!r})"
                )
            if _tolerates_failure(step) or _tolerates_failure(job):
                errors.append(f"{where}: a gate step may not set continue-on-error")
            if "if" in step and not (
                isinstance(step["if"], str) and _BACKEND_CONDITION.fullmatch(step["if"])
            ):
                errors.append(
                    f"{where}: a gate step's only allowed condition is matrix.backend "
                    f"(got {step['if']!r})"
                )
        if job_has_gate and "if" in job:
            errors.append(
                f"{path}: jobs.{job_name}: a job with gate steps may not carry a condition "
                f"(got {job['if']!r})"
            )
    return errors


def _validate_r5_structure(path: Path, makefile: Path) -> list[str]:
    """Check runners, Python versions and make targets on the parsed jobs.

    The text rules in :func:`validate_r5_policy` accepted a value that appeared
    anywhere in the file, so a ``windows-2022`` job, a Python 3.13 matrix entry
    or a make target that does not exist passed as long as the expected words
    occurred somewhere else (audit S7-09). These rules read the job structure.
    """
    try:
        data = yaml.load(path.read_text(encoding="utf-8"), Loader=_UniqueKeyLoader)
    except (OSError, UnicodeError, yaml.YAMLError):
        return []  # validate_workflow reports unreadable or invalid YAML
    if not isinstance(data, dict) or not isinstance(data.get("jobs"), dict):
        return []

    errors: list[str] = []
    versions: set[str] = set()
    make_targets: set[str] = set()
    for job_name, job in data["jobs"].items():
        if not isinstance(job, dict):
            continue
        runners = _as_list(job.get("runs-on"))
        if not all(isinstance(runner, str) and runner.startswith("ubuntu-") for runner in runners):
            errors.append(
                f"{path}: jobs.{job_name}.runs-on must be an Ubuntu runner, "
                f"got {job.get('runs-on')!r}"
            )
        strategy = job.get("strategy")
        matrix = strategy.get("matrix") if isinstance(strategy, dict) else None
        if isinstance(matrix, dict) and "python-version" in matrix:
            versions.update(str(value) for value in _as_list(matrix["python-version"]))
        steps = job.get("steps")
        for step in steps if isinstance(steps, list) else []:
            if not isinstance(step, dict):
                continue
            uses = step.get("uses")
            settings = step.get("with")
            if (
                isinstance(uses, str)
                and uses.startswith("actions/setup-python")
                and isinstance(settings, dict)
            ):
                version = settings.get("python-version")
                if isinstance(version, str) and "${{" not in version:
                    versions.add(version)
            script = step.get("run")
            if isinstance(script, str):
                make_targets.update(_make_invocations(script))

    for version in sorted(versions - set(SUPPORTED_PYTHONS)):
        errors.append(f"{path}: unsupported Python {version} is configured")
    if path.name in {"ci.yml", "wheels.yml"}:
        for version in SUPPORTED_PYTHONS:
            if version not in versions:
                errors.append(f"{path}: supported Python {version} is not configured")
    for target in sorted(make_targets - _makefile_targets(makefile)):
        errors.append(f"{path}: make target {target!r} is not defined in the Makefile")
    if path.name == "ci.yml" and "release-check" not in make_targets:
        errors.append(f"{path}: consistency job must run make release-check")
    if path.name == "ci.yml":
        errors.extend(_newest_dependency_errors(path, data["jobs"]))
    errors.extend(_trigger_errors(path, data))
    errors.extend(_neutraliser_errors(path, data))
    return errors


def validate_r5_policy(path: Path, makefile: Path = ROOT / "Makefile") -> list[str]:
    """Return active Ubuntu/version/backend policy errors for one workflow."""

    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        return [f"{path}: {exc}"]
    lowered = text.lower()
    errors: list[str] = []
    forbidden = ("windows-latest", "macos-latest", "nanobind", "follow-imports")
    for value in forbidden:
        if value in lowered:
            errors.append(f"{path}: forbidden active workflow value {value!r}")
    for version in ("3.9", "3.10"):
        if re.search(rf"(?<!\d){re.escape(version)}(?!\d)", text):
            errors.append(f"{path}: unsupported Python {version} is active")
    if "ubuntu-latest" not in lowered:
        errors.append(f"{path}: every maintained job must use Ubuntu")

    errors.extend(_validate_r5_structure(path, makefile))
    if path.name == "wheels.yml":
        if "TRAMINEX_BUILD_PYBIND11=ON" not in text:
            errors.append(f"{path}: explicit pybind11 build is missing")
        if "accel_openmp_enabled" not in text:
            errors.append(f"{path}: OpenMP wheel assertion is missing")
    # A built wheel installed with --no-deps has none of its runtime
    # dependencies, so the step that then imports traminex.distances fails on
    # numpy (audit S7-10).
    for number, line in enumerate(text.splitlines(), start=1):
        if re.search(r"\bpip\s+install\b", line) and "--no-deps" in line and ".whl" in line:
            errors.append(
                f"{path}:{number}: installs a built wheel with --no-deps, so importing "
                "traminex afterwards runs without its runtime dependencies "
                "(numpy, pandas, scipy)"
            )
    if path.name == "render-docs.yml":
        for permission in ("pages: write", "id-token: write"):
            if permission not in lowered:
                errors.append(f"{path}: Pages permission {permission!r} is missing")
    return errors


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", type=Path, help="workflow YAML paths")
    args = parser.parse_args(argv)
    paths = tuple(args.paths) or DEFAULT_WORKFLOWS
    errors = [
        error for path in paths for error in (*validate_workflow(path), *validate_r5_policy(path))
    ]
    if errors:
        for error in errors:
            print(error)
        return 1
    print(f"Workflow validation PASSED ({len(paths)} files).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
