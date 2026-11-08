"""The workflow validator rejects a CI workflow that can no longer fail (TX-AUD-20261010-057).

A ci.yml left with only the ``workflow_dispatch`` trigger and both test commands
turned into ``--collect-only || true`` passed validation: the checks read the
job structure and make targets but not whether the gates run on pushes or can
fail. ci.yml must now trigger on push and pull_request, and no gate step (a
make target, pytest, run_pybind11_tests.sh) may end in ``|| true``, collect
tests without running them, or carry ``continue-on-error: true``.

That first rule set was a short deny-list (batch-5 review M-3): ``|| echo``,
``; true``, ``set +e``, ``continue-on-error: ${{ true }}``, ``if: false``,
``-k no_such_test`` and a push filter naming a branch that never exists all
passed. Now no ``||`` may follow a gate command, nor ``; true``, ``; :`` or
``; exit 0``; a pipe after one needs ``set -o pipefail`` or ``shell: bash``;
a gate step may not run ``set +e`` or under another shell; ``continue-on-error``
must be absent or false; a gate step's ``if:`` may only pick the matrix
backend, and a job holding gate steps may have none; pytest may not select
tests with ``-k`` or ``--deselect``; and the push and pull_request filters of
ci.yml must keep the default branch (main or master).

That check still read each line of a script on its own and knew only ``-k`` and
``--deselect`` (TX-AUD-20261019-021): continued lines are now read joined, as the
shell reads them, and ``-m``, ``--ignore``, ``--ignore-glob`` and test paths are
selection too; a gate writes an option's value attached (``--option=value``).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import validate_workflows  # noqa: E402

WORKFLOWS = REPO / ".github" / "workflows"
CI = (WORKFLOWS / "ci.yml").read_text(encoding="utf-8")
PYTEST = "run: python -m pytest -q --cov=traminex --cov-branch --cov-report=term-missing"


def _errors(path: Path) -> list[str]:
    return [
        *validate_workflows.validate_workflow(path),
        *validate_workflows.validate_r5_policy(path),
    ]


def _mutated(tmp_path: Path, old: str, new: str) -> Path:
    assert CI.count(old) == 1
    path = tmp_path / "ci.yml"
    path.write_text(CI.replace(old, new), encoding="utf-8")
    return path


@pytest.mark.parametrize("name", ["ci.yml", "wheels.yml", "render-docs.yml"])
def test_the_tracked_workflows_pass(name: str) -> None:
    assert _errors(WORKFLOWS / name) == []


@pytest.mark.parametrize(
    "old, new, expected",
    [
        (
            "on:\n  push:\n    branches: [main, master]\n  pull_request:\n"
            "    branches: [main, master]\n  workflow_dispatch:\n",
            "on:\n  workflow_dispatch:\n",
            "must run on push and pull_request",
        ),
        (PYTEST, PYTEST + " || true", "cannot fail"),
        (PYTEST, "run: python -m pytest --collect-only", "collects tests without running them"),
        (
            "      - name: Run release consistency gates\n",
            "      - name: Run release consistency gates\n        continue-on-error: true\n",
            "continue-on-error",
        ),
    ],
    ids=["dispatch-only", "or-true", "collect-only", "continue-on-error"],
)
def test_a_neutralised_ci_fails(tmp_path, old, new, expected) -> None:
    errors = _errors(_mutated(tmp_path, old, new))
    assert any(expected in error for error in errors), errors


RELEASE_STEP = "      - name: Run release consistency gates\n"
RELEASE_RUN = "        run: make release-check PYTHON=python\n"
PUSH = "  push:\n    branches: [main, master]\n"


@pytest.mark.parametrize(
    "old, new, expected",
    [
        (PYTEST, PYTEST + " || echo failed", "cannot fail"),
        (PYTEST, PYTEST + "; true", "cannot fail"),
        (PYTEST, PYTEST + " | tee pytest.log", "pipe"),
        (
            RELEASE_RUN,
            "        run: |\n          set +e\n          make release-check PYTHON=python\n",
            "set +e",
        ),
        (RELEASE_STEP, RELEASE_STEP + "        shell: bash {0}\n", "shell"),
        (
            RELEASE_STEP,
            RELEASE_STEP + "        continue-on-error: ${{ true }}\n",
            "continue-on-error",
        ),
        (RELEASE_STEP, RELEASE_STEP + "        if: false\n", "condition"),
        (
            "  consistency:\n    runs-on: ubuntu-latest\n",
            "  consistency:\n    if: false\n    runs-on: ubuntu-latest\n",
            "condition",
        ),
        (PYTEST, PYTEST + " -k no_such_test_name", "selects"),
        (PUSH, "  push:\n    branches: [never-used]\n", "default branch"),
        (PUSH, "  push:\n    branches-ignore: [main, master]\n", "default branch"),
    ],
    ids=[
        "or-echo",
        "semicolon-true",
        "pipe",
        "set-plus-e",
        "custom-shell",
        "continue-on-error-expression",
        "step-if-false",
        "job-if-false",
        "keyword-filter",
        "push-other-branch",
        "push-ignores-default-branch",
    ],
)
def test_a_weakened_gate_fails(tmp_path, old, new, expected) -> None:
    errors = _errors(_mutated(tmp_path, old, new))
    assert any(expected in error for error in errors), errors


@pytest.mark.parametrize(
    "old, new",
    [
        (RELEASE_STEP, RELEASE_STEP + "        continue-on-error: false\n"),
        (RELEASE_STEP, RELEASE_STEP + "        shell: bash\n"),
        (
            PYTEST,
            "run: |\n          set -o pipefail\n"
            "          python -m pytest -q --cov=traminex --cov-branch | tee pytest.log",
        ),
        (PUSH, "  push:\n    branches: [master]\n"),
    ],
    ids=["continue-on-error-false", "shell-bash", "pipe-with-pipefail", "push-master"],
)
def test_a_gate_that_still_fails_passes(tmp_path, old, new) -> None:
    assert _errors(_mutated(tmp_path, old, new)) == []


#: The audit's neutralised copies A and C (TX-AUD-20261019-021) passed validation: tests
#: selected by marker, ignore and path, and ``|| true`` on a continuation line, which the
#: validator read as a line of its own. Copy B's ``|| true`` after a continued make line was
#: reported only as an undefined make target ``\``.
AUDIT_COPY_A = (
    "run: python -m pytest -q -m 'not slow and not r' --ignore=tests tests/test_version.py"
)
AUDIT_COPY_B = (
    "        run: |\n          make release-check PYTHON=python \\\n            || true\n"
)
AUDIT_COPY_C = (
    "run: |\n          python -m pytest -q --cov=traminex \\\n            --cov-branch || true"
)
PYBIND11 = 'run: TRAMINEX_PYTHON="$(command -v python)" scripts/run_pybind11_tests.sh -q'


@pytest.mark.parametrize(
    "old, new, expected",
    [
        (PYTEST, AUDIT_COPY_A, "selects"),
        (RELEASE_RUN, AUDIT_COPY_B, "cannot fail"),
        (PYTEST, AUDIT_COPY_C, "cannot fail"),
        (PYTEST, "run: python -m pytest -q -m 'not slow'", "selects"),
        (PYTEST, "run: python -m pytest -q -mslow", "selects"),
        (PYTEST, "run: python -m pytest -q --ignore=tests/test_version.py", "selects"),
        (PYTEST, "run: python -m pytest -q --ignore-glob='*goldens*'", "selects"),
        (PYTEST, "run: python -m pytest -q tests/test_version.py", "selects"),
        (PYTEST, "run: python -m pytest -q tests/test_version.py::test_version", "selects"),
        (PYTEST, "run: python -m pytest -q --co", "collects tests without running them"),
        (PYBIND11, PYBIND11 + " tests/test_version.py", "selects"),
    ],
    ids=[
        "audit-copy-a",
        "audit-copy-b",
        "audit-copy-c",
        "marker",
        "marker-attached",
        "ignore",
        "ignore-glob",
        "path",
        "node-id",
        "co",
        "pybind11-path",
    ],
)
def test_selection_and_continued_neutralisers_fail(tmp_path, old, new, expected) -> None:
    errors = _errors(_mutated(tmp_path, old, new))
    assert any(expected in error for error in errors), errors


@pytest.mark.parametrize(
    "old, new",
    [
        (
            PYTEST,
            "run: |\n          python -m pytest -q --cov=traminex \\\n"
            "            --cov-branch --cov-report=term-missing",
        ),
        (
            RELEASE_RUN,
            "        run: |\n          make release-check \\\n            PYTHON=python\n",
        ),
        (PYTEST, PYTEST + " 2>&1"),
        (PYTEST, PYTEST + " && python -m coverage report"),
    ],
    ids=["continued-pytest", "continued-make", "redirection", "later-python-m"],
)
def test_continued_and_option_only_gates_pass(tmp_path, old, new) -> None:
    assert _errors(_mutated(tmp_path, old, new)) == []


#: Batch-8 review I-3: the test command was recognised only after a literal ``python``, so
#: these copies still passed; PYTEST_ADDOPTS, options that run no test and options that
#: rewrite pytest's configuration passed as well.
@pytest.mark.parametrize(
    "old, new, expected",
    [
        (PYTEST, "run: $PYTHON -m pytest -q tests/test_version.py", "selects"),
        (PYTEST, "run: ${{ env.PY }} -m pytest -q tests/test_version.py", "selects"),
        (PYTEST, "run: coverage run -m pytest -q tests/test_version.py", "selects"),
        (PYTEST, "run: env MPLBACKEND=Agg pytest -q tests/test_version.py", "selects"),
        (PYTEST, "run: PYTEST_ADDOPTS='-k version' python -m pytest -q", "PYTEST_ADDOPTS"),
        (
            "          TRAMINEX_BINDING: python\n",
            "          TRAMINEX_BINDING: python\n          PYTEST_ADDOPTS: -k version\n",
            "PYTEST_ADDOPTS",
        ),
        (PYTEST, "run: python -m pytest -q --setup-plan", "runs no test"),
        (PYTEST, "run: python -m pytest -q --fixtures", "runs no test"),
        (PYTEST, "run: python -m pytest --version", "runs no test"),
        (PYTEST, "run: python -m pytest -q -oaddopts=-kversion", "configuration"),
        (PYTEST, "run: python -m pytest -q --override-ini=testpaths=tests/x", "configuration"),
        (PYTEST, "run: python -m pytest -q -ctox.ini", "configuration"),
    ],
    ids=[
        "variable-interpreter",
        "expression-interpreter",
        "coverage-run",
        "env-prefix",
        "addopts-inline",
        "addopts-step-env",
        "setup-plan",
        "fixtures",
        "version",
        "override-ini-attached",
        "override-ini-long",
        "config-file-attached",
    ],
)
def test_review_bypasses_fail(tmp_path, old, new, expected) -> None:
    errors = _errors(_mutated(tmp_path, old, new))
    assert any(expected in error for error in errors), errors


@pytest.mark.parametrize(
    "old, new",
    [
        (PYTEST, "run: $PYTHON -m pytest -q --cov=traminex --cov-branch --cov-report=term-missing"),
        (
            "      - name: Install test profile\n        run: |\n",
            "      - name: Install test profile\n        run: |\n"
            "          python -m pip install pytest hypothesis\n",
        ),
    ],
    ids=["variable-interpreter-whole-suite", "install-names-pytest"],
)
def test_whole_suite_commands_still_pass(tmp_path, old, new) -> None:
    assert _errors(_mutated(tmp_path, old, new)) == []
