# Ubuntu-native development and release gates.

.PHONY: audit dev-environment-check docs docs-check docs-render examples examples-lab format-check \
	generated-check lint \
	install-profiles release-check sdist-manifest stub-check stubs test-pybind11 test-python \
	typecheck workflow-check

PYTHON ?= $(HOME)/miniforge3/envs/datasci/bin/python

# `install-profiles` passes these to a verifier that requires an existing file, so
# a bare name that is not on PATH fails argparse immediately (cycle-4b F-017).
# Resolve from PATH when available and otherwise fall back to this machine's known
# location, consistent with PYTHON above. `dot` ships in the same environment as
# PYTHON, so its fallback is derived rather than repeated.
PYTHON311 ?= $(shell command -v python3.11 2>/dev/null || \
	echo $(HOME)/.local/share/uv/python/cpython-3.11-linux-x86_64-gnu/bin/python3.11)
UV ?= $(shell command -v uv 2>/dev/null || echo $(HOME)/.local/bin/uv)
DOT ?= $(shell command -v dot 2>/dev/null || echo $(dir $(PYTHON))dot)
SDIST ?=
PURE_WHEEL ?=
COMPILED_WHEEL ?=
# Under the git-ignored build/ (TX-AUD-20261010-088); the verifier creates the directory.
INSTALL_PROFILE_EVIDENCE ?= build/install-profiles.json

examples:
	$(PYTHON) scripts/open_examples.py

examples-lab:
	$(PYTHON) scripts/open_examples.py --lab

docs:
	@echo "Open docs/USAGE_new_apis.md in your viewer."
	@echo " - $(CURDIR)/docs/USAGE_new_apis.md"

docs-render:
	cd docs && quarto render

docs-check:
	$(PYTHON) scripts/check_docs_consistency.py

generated-check:
	$(PYTHON) scripts/check_generated_artifacts.py

# Names any missing development module and the install command (TX-AUD-20261010-046).
dev-environment-check:
	$(PYTHON) scripts/check_dev_environment.py

format-check:
	$(PYTHON) scripts/check_dev_environment.py
	$(PYTHON) -m black --check traminex tests scripts benchmarks examples docs \
		verify_audit_improvements.py

# Same paths as format-check, so ruff sees every tracked Python file (S7-22).
lint:
	$(PYTHON) -m ruff check traminex tests scripts benchmarks examples docs \
		verify_audit_improvements.py

typecheck:
	$(PYTHON) -m mypy

workflow-check:
	$(PYTHON) scripts/validate_workflows.py

# The last line type-checks the stub tree itself (audit S7-01): MYPYPATH puts
# stubs/ ahead of the runtime package, so `-p traminex` checks the .pyi files.
stub-check:
	$(PYTHON) scripts/sync_stub_public_api.py --check
	$(PYTHON) scripts/check_stub_signatures.py
	MYPYPATH=stubs $(PYTHON) -m mypy --cache-dir .mypy_cache/stubs -p traminex

stubs:
	$(PYTHON) scripts/sync_stub_public_api.py
	$(PYTHON) scripts/check_stub_signatures.py

# Rewrites pyproject.toml's sdist.include from `git ls-files` and the policy in the
# script; tests/test_packaging_profiles.py fails while the list is out of date.
sdist-manifest:
	$(PYTHON) scripts/sync_sdist_manifest.py

test-python:
	TRAMINEX_BINDING=python $(PYTHON) -m pytest -q

test-pybind11:
	TRAMINEX_PYTHON="$(PYTHON)" scripts/run_pybind11_tests.sh -q

install-profiles:
	@test -n "$(SDIST)" -a -n "$(PURE_WHEEL)" -a -n "$(COMPILED_WHEEL)"
	$(PYTHON) scripts/verify_install_profiles.py \
		--sdist "$(SDIST)" --pure-wheel "$(PURE_WHEEL)" \
		--compiled-wheel "$(COMPILED_WHEEL)" \
		--python311 "$(PYTHON311)" --python312 "$(PYTHON)" \
		--uv "$(UV)" --dot "$(DOT)" \
		--evidence "$(INSTALL_PROFILE_EVIDENCE)"

release-check: dev-environment-check format-check lint typecheck workflow-check docs-check stub-check \
	generated-check
	TRAMINEX_MAKE_PYTHON="$(PYTHON)" $(PYTHON) -m pytest tests/test_supported_versions.py \
		tests/test_version_single_source.py tests/test_public_api_dupes.py \
		tests/test_semantic_parity_ledger.py -q
	$(PYTHON) verify_audit_improvements.py

# Every artifact docs/DOCS_POLICY.md says `make audit` refreshes (TX-AUD-20261010-086);
# `make generated-check` compares them with fresh copies.
audit:
	$(PYTHON) scripts/audit_traminex.py
	$(PYTHON) scripts/generate_parity_table.py
	$(PYTHON) scripts/generate_api_mappings.py
	$(PYTHON) scripts/generate_semantic_parity_ledger.py
