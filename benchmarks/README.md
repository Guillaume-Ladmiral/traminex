# Reproducible benchmarks

`run_benchmarks.py` is the maintained benchmark entry point. It replaces the
former `benchmark_accelerators.py`, `benchmark_p2_accelerators.py`,
`bindings_compare.py`, and `python_baseline.py` scripts.

Run a tiny agreement smoke from the repository environment:

```bash
PY=~/miniforge3/envs/datasci/bin/python
$PY -m benchmarks.run_benchmarks --method OM --backend python \
  --backend numba --N 8 --L 6 --K 3 --warmups 1 --repeats 2
```

The supported backends are `python`, `numba`, `pybind11`, and `auto`. `auto`
prefers pybind11, then Numba, then pure Python for the selected method. A
backend explicitly requested but not installed or not valid for that method is
an error; it is never relabeled as a successful fallback.

Use repeated `--backend` options to compare implementations. The common method
set is `DHD`, `DUR`, `HAM`, `JSD_TRANS`, `LCP`, `LCS`, `OM`, and `RLCP`. Numba
does not accelerate `HAM`; pybind11 covers neither `DUR` nor `JSD_TRANS`. The DHD
comparison intentionally uses the legacy mean-frequency variant, which is the one
contract shared by Python, Numba, and pybind11. This does not change the default
TraMineR-faithful DHD behavior. `DUR` is benchmarked with `metric="chi2"`,
because that is the only `DUR` metric with a Numba kernel.

`DUR`-chi2 and `JSD_TRANS` were added deliberately: they are the two library
Numba routes whose kernels reassociate logs and divisions, so they are the only
ones where agreement is not byte-identical by construction. Without them the
agreement-before-timing contract exercised no non-trivial case (cycle-4b F-016).
Measured on the reference fixture, `DUR`-chi2 is in fact byte-identical and
`JSD_TRANS` differs by ~1e-16, comfortably inside the tolerance below.

Agreement is evaluated before warmup and timing. Promised pybind11 kernels must
be byte-identical to Python with `numpy.array_equal`; Numba uses the established
`rtol=1e-10`, `atol=1e-12` comparison. No benchmark applies a performance
threshold. The strict per-site byte-identity comparison lives in
`tests/test_numba_backend_agreement.py`, which drives both backends through the
supported `TRAMINEX_NUMBA` switch.

Write a reproducibility report with:

```bash
$PY -m benchmarks.run_benchmarks --method LCS --backend auto \
  --N 200 --L 30 --K 5 --seed 42 --warmups 1 --repeats 5 \
  --json /tmp/traminex-lcs-benchmark.json
```

The JSON records the Python version, the imported traminex (`metadata.package`:
its `__file__`, `__version__` and `source_sha256`, a SHA-256 of its `.py`
sources), the kernel options timed (`config.kernel_options`, for example
`{"variant": "mean"}` for DHD), resolved backend, platform, CPU, OpenMP and
thread information, deterministic input identity, every timing sample, summary
statistics, and the separate agreement result. R5 audit timing
reports belong only under the ignored reliability evidence directory.

For a compiled smoke, first make `_pb_kernels` importable with the maintained
pybind11/OpenMP build, then require it explicitly:

```bash
TRAMINEX_BINDING=pybind11 $PY -m benchmarks.run_benchmarks \
  --method OM --backend python --backend pybind11 \
  --N 8 --L 6 --K 3 --warmups 1 --repeats 2
```

The CLI exits nonzero if a requested backend is missing, if the method/backend
pair is unsupported, or if agreement fails. Timing samples are descriptive
evidence; they are intentionally not a release threshold.
