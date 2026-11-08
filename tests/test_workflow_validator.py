from __future__ import annotations

from pathlib import Path

from scripts.validate_workflows import DEFAULT_WORKFLOWS, validate_r5_policy, validate_workflow


def test_workflow_validator_accepts_minimal_workflow(tmp_path: Path):
    workflow = tmp_path / "valid.yml"
    workflow.write_text(
        "name: Example\non: workflow_dispatch\njobs:\n  test:\n"
        "    runs-on: ubuntu-latest\n    steps:\n      - run: 'true'\n",
        encoding="utf-8",
    )

    assert validate_workflow(workflow) == []


def test_workflow_validator_reports_duplicate_keys(tmp_path: Path):
    workflow = tmp_path / "invalid.yml"
    workflow.write_text(
        "name: Example\non: workflow_dispatch\njobs:\n  test:\n    steps: []\n    steps: []\n",
        encoding="utf-8",
    )

    errors = validate_workflow(workflow)

    assert any("duplicate key 'steps'" in error for error in errors)


def test_workflow_validator_reports_missing_runner(tmp_path: Path):
    workflow = tmp_path / "invalid.yml"
    workflow.write_text(
        "name: Example\non: workflow_dispatch\njobs:\n  test:\n"
        "    steps:\n      - run: 'true'\n",
        encoding="utf-8",
    )

    errors = validate_workflow(workflow)

    assert any("jobs.test.runs-on" in error for error in errors)


def test_repository_workflows_satisfy_r5_ubuntu_policy():
    errors = [error for path in DEFAULT_WORKFLOWS for error in validate_r5_policy(path)]

    assert errors == []


_WHEELS_TEMPLATE = """name: Wheels
on: workflow_dispatch
jobs:
  build-wheels:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.11", "3.12"]
    steps:
      - uses: actions/setup-python@v5
        with:
          python-version: ${{{{ matrix.python-version }}}}
      - run: python -m pip install build cmake ninja pybind11 scikit-build-core
      - run: echo "-DTRAMINEX_BUILD_PYBIND11=ON"
      - run: |
          {install}
          python -c "from traminex import _pb_kernels, distances; _pb_kernels.accel_openmp_enabled()"
"""


def test_wheel_verification_must_install_the_runtime_dependencies(tmp_path: Path):
    """``wheels.yml`` verified its wheel in an environment without numpy (S7-10).

    The step installed the wheel with ``--no-deps`` after installing only build
    tools (build, cmake, ninja, pybind11, scikit-build-core, none of which pulls
    numpy), then imported ``traminex.distances``, which imports numpy, pandas
    and scipy at module level (audit S7-10, 2026-09-22; reproduced in a fresh
    environment).
    """
    workflow = tmp_path / "wheels.yml"
    workflow.write_text(
        _WHEELS_TEMPLATE.format(
            install="python -m pip install --force-reinstall --no-deps dist/*.whl"
        ),
        encoding="utf-8",
    )

    errors = validate_r5_policy(workflow)

    assert any("--no-deps" in error for error in errors), errors

    workflow.write_text(
        _WHEELS_TEMPLATE.format(install="python -m pip install --force-reinstall dist/*.whl"),
        encoding="utf-8",
    )
    assert not [error for error in validate_r5_policy(workflow) if "--no-deps" in error]


_MAKEFILE = "release-check: lint\n\ttrue\n\nlint:\n\ttrue\n\ntypecheck:\n\ttrue\n"


def _policy_errors(tmp_path: Path, name: str, text: str) -> list[str]:
    workflow = tmp_path / name
    workflow.write_text(text, encoding="utf-8")
    makefile = tmp_path / "Makefile"
    makefile.write_text(_MAKEFILE, encoding="utf-8")
    return validate_r5_policy(workflow, makefile=makefile)


def test_policy_reads_runners_versions_and_make_targets_from_the_jobs(tmp_path: Path):
    """The R5 policy used to accept any file that merely contained its words (S7-09).

    This workflow, from the audit probe, has a ``windows-2022`` job, a Python
    3.13 matrix and a make target that does not exist; it passed because
    ``ubuntu-latest``, ``3.11``, ``3.12`` and ``make release-check`` all
    appeared somewhere in the text (audit S7-09, 2026-09-22).
    """
    errors = _policy_errors(
        tmp_path,
        "ci.yml",
        "name: CI\non: [push]\njobs:\n"
        "  consistency:\n    runs-on: ubuntu-latest\n    steps:\n"
        "      - run: make release-check && echo 3.11 3.12\n"
        "  test:\n    runs-on: windows-2022\n    strategy:\n      matrix:\n"
        '        python-version: ["3.13"]\n    steps:\n'
        "      - run: make this-target-does-not-exist\n",
    )

    assert any("jobs.test.runs-on must be an Ubuntu runner" in e for e in errors), errors
    assert any("unsupported Python 3.13 is configured" in e for e in errors), errors
    assert any("supported Python 3.11 is not configured" in e for e in errors), errors
    assert any("supported Python 3.12 is not configured" in e for e in errors), errors
    assert any(
        "make target 'this-target-does-not-exist' is not defined" in e for e in errors
    ), errors
    assert not any("must run make release-check" in e for e in errors), errors


def test_policy_requires_ci_to_invoke_release_check_as_a_command(tmp_path: Path):
    errors = _policy_errors(
        tmp_path,
        "ci.yml",
        "name: CI\non: [push]\njobs:\n  consistency:\n    runs-on: ubuntu-latest\n"
        '    strategy:\n      matrix:\n        python-version: ["3.11", "3.12"]\n'
        "    steps:\n      - run: echo make release-check\n",
    )

    assert any("must run make release-check" in e for e in errors), errors


def test_policy_accepts_setup_python_versions_and_make_variables(tmp_path: Path):
    errors = _policy_errors(
        tmp_path,
        "ci.yml",
        "name: CI\non: [push, pull_request]\njobs:\n  consistency:\n    runs-on: ubuntu-latest\n"
        "    steps:\n      - uses: actions/setup-python@v5\n        with:\n"
        '          python-version: "3.11"\n'
        "      - uses: actions/setup-python@v5\n        with:\n"
        '          python-version: "3.12"\n'
        "      - run: make release-check PYTHON=python\n"
        # ci.yml must also type-check with the newest dependencies (TX-AUD-20261010-010).
        "  typecheck-newest:\n    runs-on: ubuntu-latest\n"
        '    strategy:\n      matrix:\n        python-version: ["3.11", "3.12"]\n'
        "    steps:\n      - uses: actions/setup-python@v5\n"
        "        with:\n          python-version: ${{ matrix.python-version }}\n"
        "      - run: python -m pip install --upgrade numpy pandas scipy\n"
        "      - run: make typecheck PYTHON=python\n",
    )

    assert errors == []
