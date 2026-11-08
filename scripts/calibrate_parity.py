#!/usr/bin/env python
import argparse
import json
import os
from pathlib import Path

import numpy as np

from traminex import seqdef
from traminex.distances.costs import trate_substitution_matrix
from traminex.distances.python_kernels import pairwise as py_pairwise

# Defaults are anchored to this repository, not the working directory (S7-18).
REPO = Path(__file__).resolve().parents[1]
DEFAULT_REFS = REPO / "tests" / "references"
DEFAULT_OUT = DEFAULT_REFS / "calibrated_params.json"

# A candidate is recorded only when every value is within these absolute tolerances of
# the R reference (those of tests/references/README.md, without np.allclose's relative
# part, which let a TRATE cost near 2 pass about 2e-5 off; batch-17 review M-1).
# Otherwise nothing is recorded (TX-AUD-20261019-087). The written file is what
# seqdist(..., trate_profile="calibrated") and dhd_profile="calibrated" read.
TRATE_ATOL = 1e-6
DHD_ATOL = 2e-3


def _reproduces(candidate, ref, atol):
    return candidate.shape == ref.shape and bool(np.allclose(candidate, ref, rtol=0.0, atol=atol))


def calibrate_trate(x, ref_path):
    """The TRATE options that reproduce the R reference, or None when none does."""
    ref = np.loadtxt(ref_path, delimiter=",")
    best = None
    for mode in ["symmetric", "row", "col"]:
        for normalize in ["none", "mean_offdiag"]:
            for smoothing in [0.0, 1e-9, 1e-6, 1e-4]:
                sm = trate_substitution_matrix(
                    x.data, smoothing=smoothing, mode=mode, normalize=normalize
                )
                if not _reproduces(sm, ref, TRATE_ATOL):
                    continue
                err = float(np.mean((sm - ref) ** 2))
                cand = (err, {"smoothing": smoothing, "mode": mode, "normalize": normalize})
                if best is None or cand[0] < best[0]:
                    best = cand
    return best[1] if best else None


def calibrate_dhd(x, ref_path):
    """The frequency-variant options that reproduce the R reference, or None.

    R's default DHD uses transition costs, which the public seqdist(method="DHD")
    implements; the low-level frequency variants are a different distance and
    reproduce it only by coincidence.
    """
    ref = np.loadtxt(ref_path, delimiter=",")
    best = None
    for variant in ["mean", "sum"]:
        for smooth in [None, 3, 5]:
            D = py_pairwise(x.data[:2], method="DHD", variant=variant, smooth_window=smooth)
            if not _reproduces(D, ref, DHD_ATOL):
                continue
            err = float(np.mean((D - ref) ** 2))
            cand = (err, {"variant": variant, "smooth_window": smooth})
            if best is None or cand[0] < best[0]:
                best = cand
    return best[1] if best else None


def main():
    p = argparse.ArgumentParser(
        description=(
            "Calibrate TRATE and DHD to reference CSVs. Writes calibrated_params.json with a"
            " section only for parameters that reproduce the R reference (every value within"
            " an absolute 1e-6 for TRATE, 2e-3 for DHD)."
        )
    )
    p.add_argument("--refs", default=str(DEFAULT_REFS), help="Path to references directory")
    p.add_argument("--out", default=str(DEFAULT_OUT))
    args = p.parse_args()

    input_path = os.path.join(args.refs, "input_wide.csv")
    trate_ref = os.path.join(args.refs, "expected_trate_subm.csv")
    dhd_ref = os.path.join(args.refs, "expected_dhd.csv")
    if not os.path.exists(input_path):
        print("No input_wide.csv; nothing to calibrate")
        return
    import pandas as pd

    df = pd.read_csv(input_path)
    x = seqdef(df)
    out = {}
    if os.path.exists(trate_ref):
        trate = calibrate_trate(x, trate_ref)
        if trate is None:
            print(f"TRATE: no candidate reproduces the R reference within {TRATE_ATOL:g}.")
        else:
            out["trate"] = trate
            print("Calibrated TRATE:", trate)
    if os.path.exists(dhd_ref):
        dhd = calibrate_dhd(x, dhd_ref)
        if dhd is None:
            print(
                f"DHD: no frequency variant reproduces the R reference within {DHD_ATOL:g};"
                ' R\'s default DHD is the public seqdist(method="DHD")'
                " (tests/test_parity_refs.py compares it)."
            )
        else:
            out["dhd"] = dhd
            print("Calibrated DHD:", dhd)
    if not (os.path.exists(trate_ref) or os.path.exists(dhd_ref)):
        print("No reference files found to calibrate")
        return
    # The file is written whenever a reference was read, even empty, so an earlier
    # file's sections never outlive the calibration that replaced them (batch-17 review M-2).
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    if out:
        print("Wrote", args.out)
    else:
        print(
            f"No parameters reproduce the R references; wrote an empty calibration to {args.out},"
            " which a calibrated profile refuses."
        )


if __name__ == "__main__":
    main()
