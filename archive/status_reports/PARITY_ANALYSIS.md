# Parity Analysis

This file is a high-level companion to the machine-generated parity artifacts.

## Source of Truth

- Audit report (human-readable): `reports/audit/traminex_api_audit.md`
- Audit report (machine-readable): `reports/audit/traminex_api_audit.json`
- Function table for docs: `docs/API_PARITY_TABLE.md`
- Narrative parity overview: `docs/API_PARITY.md`

## Current Snapshot

Use the generated audit report for exact counts and mappings. Do not rely on manually maintained totals in this file.

## Regeneration

```powershell
powershell -ExecutionPolicy Bypass -File scripts/audit_traminex.ps1
.\scripts\project_python.ps1 scripts/generate_parity_table.py
.\scripts\project_python.ps1 scripts/generate_api_mappings.py
.\scripts\project_python.ps1 scripts/annotate_traminer_list.py
```

## Notes

- Historical parity narratives and milestone-era analyses are archived under `archive/status_reports/`.
- If parity numbers conflict between documents, trust `reports/audit/traminex_api_audit.json`.
