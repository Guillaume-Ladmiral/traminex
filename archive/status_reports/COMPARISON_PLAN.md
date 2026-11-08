# M2 binding comparison plan (nanobind vs pybind11 vs Cython)

Goal
- Implement and benchmark C++ kernels for OM and HAM, compare nanobind vs pybind11 (and a minimal Cython) for our workloads.
- Default to nanobind builds if results show a clear advantage; keep a pybind11 path switchable.

Scope (M2)
- C++17 reference implementations for:
  - OM (edit distance, uniform or matrix substitution)
  - HAM (Hamming distance)
- Three bindings:
  - nanobind module: `traminex._nb_kernels`
  - pybind11 module: `traminex._pb_kernels`
  - Cython fallback: `traminex._cy_kernels` (optional, useful for packaging edge cases)

Build system
- scikit-build-core + CMake + Ninja
- Wheel matrix for Python 3.9–3.12 on Windows/macOS/Linux (x86_64, arm64 macOS)

Harness outline
- benchmarks/bindings_compare.py
  - Generates random integer-coded sequences (N,L,K); repeats a few trials
  - Times pairwise OM/HAM for each backend (nb, pb, cy) and for the Python baseline
  - Emits CSV/JSON summary with medians and speedups vs Python
- CLI example:
  ```
  .venv/Scripts/python benchmarks/bindings_compare.py --N 300 --L 40 --K 8 --trials 3 --out results.json
  ```
- Optional: pin thread env (e.g., OMP_NUM_THREADS=1) to reduce variance

Selection rule
- Choose the default binding by median runtime on representative sizes:
  - N in {200, 400, 800}, L in {30, 60}, K in {5, 10}
  - If nanobind median is <= 0.9× pybind11 across most cells, set nanobind as default
  - Keep an environment switch TRAMINEX_BINDING=pybind11|nanobind|cython|python

Implementation notes
- Shared core C++ sources, separate thin binding layers
- Keep function signatures aligned across backends:
  - pairwise_om(seq: int32[N,L], indel: float, subm: scalar or float32[K,K]) -> float64[N,N]
  - pairwise_ham(seq: int32[N,L]) -> float64[N,N]
- Validate on small deterministic cases vs Python kernels

Deliverables
- C++ sources + minimal bindings
- Bench harness (+ CI job to run a small matrix on PRs)
- Short report with results and recommendation (nanobind default unless results contradict)

Risks & mitigations
- Windows packaging: verify MSVC settings early; test local wheel install
- Substitution matrix broadcasting: pre-validate shapes to avoid Python<->C++ overhead
- Memory usage for pairwise: chunked computation option for large N (later)
