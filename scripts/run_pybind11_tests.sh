#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
cd -- "$ROOT"

# make test-pybind11 passes its PYTHON (TX-AUD-20261011-078); run directly, set
# TRAMINEX_PYTHON or the python3 on PATH is used.
PY=${TRAMINEX_PYTHON:-python3}
source "$ROOT/scripts/_native_link.sh"
clear_dead_native_links "$ROOT/traminex"
PYBIND11_DIR=$($PY -m pybind11 --cmakedir)
BUILD_DIR=$(mktemp -d "${NATIVE_BUILD_PREFIX}XXXXXX")
# A later run reads this to tell a killed run's link from one in use (TX-AUD-20261010-088).
echo "$$" > "$BUILD_DIR/$NATIVE_PID_FILE"
PB_LINK=

cleanup() {
  if [[ -n "$PB_LINK" ]]; then
    rm -f -- "$PB_LINK"
  fi
  rm -rf -- "$BUILD_DIR"
}
trap cleanup EXIT

cmake -S . -B "$BUILD_DIR" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DTRAMINEX_BUILD_PYBIND11=ON \
  -DPython_EXECUTABLE="$PY" \
  -Dpybind11_DIR="$PYBIND11_DIR"
cmake --build "$BUILD_DIR" --parallel

PB_SO=$(find "$BUILD_DIR" -maxdepth 1 -type f -name '_pb_kernels*.so' -print -quit)
test -n "$PB_SO"

PB_CANDIDATE="$ROOT/traminex/$(basename "$PB_SO")"
refuse_existing_module "$PB_CANDIDATE" || exit 1

ln -s -- "$PB_SO" "$PB_CANDIDATE"
PB_LINK="$PB_CANDIDATE"

TRAMINEX_BINDING=pybind11 "$PY" - <<'PY'
from pathlib import Path

from scripts._checkout_imports import prefer_checkout_imports

prefer_checkout_imports()

import traminex
from traminex import _pb_kernels, distances

assert Path(traminex.__file__).resolve().parent == Path.cwd() / "traminex"
assert distances._impl_name() == "pybind11", distances._impl_name()
assert hasattr(_pb_kernels, "accel_openmp_enabled")
assert _pb_kernels.accel_openmp_enabled() is True
print("verified backend:", distances._impl_name())
print("OpenMP enabled:", _pb_kernels.accel_openmp_enabled())
PY

TRAMINEX_BINDING=pybind11 "$PY" - "$@" <<'PY'
import sys

from scripts._checkout_imports import prefer_checkout_imports

prefer_checkout_imports()

import pytest

from scripts._native_skip_gate import NativeSkipGate

# A compiled-kernel skip fails this gate (docs/DEV_BUILD.md; TX-AUD-20261011-075).
raise SystemExit(pytest.main(sys.argv[1:], plugins=[NativeSkipGate()]))
PY
