# Name-mapping artifacts

This directory contains discovery artifacts that compare R namespace spellings
with the Python public API:

- `mapping_traminex_vs_TraMineR.csv`
- `mapping_traminex_vs_TraMineRextras.csv`
- `traminer_API_list.csv`
- `traminer_API_list_with_traminex.csv`

Exact, sanitized, and fuzzy scores are name-mapping suggestions. They do not
establish compatible values, options, containers, errors, or behavior. The
hand-reviewed semantic source is
[`docs/API_PARITY_DECISIONS.csv`](../../docs/API_PARITY_DECISIONS.csv); its
generated views are
[`docs/API_PARITY_LEDGER.md`](../../docs/API_PARITY_LEDGER.md) and
[`reports/audit/traminex_semantic_parity.json`](../audit/traminex_semantic_parity.json).

## Prerequisites

Run from the repository root on Ubuntu with the project interpreter. The two
vendored `NAMESPACE` files and `traminex/api_manifest.py` must be present.
`annotate_traminer_list.py` additionally requires the tracked
`docs/API_PARITY_TABLE.csv` and `reports/parity/traminer_API_list.csv`.

```bash
PY=~/miniforge3/envs/datasci/bin/python
```

## Regeneration

Build the namespace/API audit first. `--output-dir` changes its destination;
`--generated-at` accepts a timezone-aware ISO-8601 value when a reproducible
audit timestamp is needed.

```bash
$PY scripts/audit_traminex.py
```

Generate the human-readable name table from that audit JSON. Both paths can be
overridden with `--audit-json` and `--output`.

```bash
$PY scripts/generate_parity_table.py
cp reports/audit/traminex_api_audit.csv docs/API_PARITY_TABLE.csv
```

Split the audit CSV into the two package mapping CSVs. The CLI also accepts
`--parity-csv`, `--extras-namespace`, and `--out-dir`.

```bash
$PY scripts/generate_api_mappings.py
```

Finally, annotate the fixed TraMineR inventory. This script has no CLI options;
it reads and writes the tracked paths named in the prerequisites.

```bash
$PY scripts/annotate_traminer_list.py
```

`scripts/check_generated_artifacts.py` regenerates and checks the audit
Markdown/JSON/CSV, `docs/API_PARITY_TABLE.md`,
`docs/API_PARITY_TABLE.csv`, the two package mapping CSVs, and the generated
semantic Markdown/JSON:

```bash
$PY scripts/check_generated_artifacts.py
```

The annotated `traminer_API_list_with_traminex.csv` is not produced by that
freshness checker; regenerate it explicitly after its two fixed inputs change.
