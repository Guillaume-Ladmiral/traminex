#!/usr/bin/env python3
import argparse
import csv
import sys
from pathlib import Path

if __package__:
    from ._r_namespace import parse_namespace
else:
    from _r_namespace import parse_namespace


def load_parity_csv(parity_path: Path):
    rows = []
    with parity_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        required = {"r_name", "py_name", "match_type", "score"}
        if not required.issubset(reader.fieldnames or []):
            raise RuntimeError(f"Unexpected columns in {parity_path}: {reader.fieldnames}")
        for row in reader:
            rows.append({k: row.get(k, "") for k in ("r_name", "py_name", "match_type", "score")})
    return rows


def extract_extras_names(namespace_path: Path) -> set[str]:
    if not namespace_path.exists():
        raise FileNotFoundError(f"TraMineRextras NAMESPACE not found: {namespace_path}")
    # The shared parser strips R comments, so a commented-out export is not
    # read as a name (audit S7-17: the old regex produced "##seqeordplot").
    return set(parse_namespace(namespace_path))


def write_csv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["r_name", "py_name", "match_type", "score"],
            lineterminator="\n",
        )
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split parity CSV into TraMineR and TraMineRextras mapping reports."
    )
    parser.add_argument(
        "--parity-csv",
        default=None,
        help="Source parity CSV. Defaults to reports/audit/traminex_api_audit.csv or docs/API_PARITY_TABLE.csv.",
    )
    parser.add_argument(
        "--extras-namespace",
        default=None,
        help="Path to SOURCES/TraMineRextras/NAMESPACE.",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Output directory for generated mapping CSVs.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    if args.parity_csv is None:
        parity_candidates = [
            repo_root / "reports" / "audit" / "traminex_api_audit.csv",
            repo_root / "docs" / "API_PARITY_TABLE.csv",
        ]
        parity_path = next((p for p in parity_candidates if p.exists()), None)
        if parity_path is None:
            raise FileNotFoundError(
                "No parity CSV found. Expected one of: "
                + ", ".join(str(p) for p in parity_candidates)
            )
    else:
        parity_path = Path(args.parity_csv)
    extras_namespace = (
        Path(args.extras_namespace)
        if args.extras_namespace is not None
        else repo_root / "SOURCES" / "TraMineRextras" / "NAMESPACE"
    )
    out_dir = Path(args.out_dir) if args.out_dir is not None else repo_root / "reports" / "parity"
    out_traminer = out_dir / "mapping_traminex_vs_TraMineR.csv"
    out_extras = out_dir / "mapping_traminex_vs_TraMineRextras.csv"

    rows = load_parity_csv(parity_path)
    extras_names = extract_extras_names(extras_namespace)

    rows_extras = [r for r in rows if r.get("r_name", "") in extras_names]
    rows_traminer = [r for r in rows if r.get("r_name", "") not in extras_names]

    write_csv(out_traminer, rows_traminer)
    write_csv(out_extras, rows_extras)

    print("Wrote:")
    print(out_traminer)
    print(out_extras)
    print(f"TraMineR rows: {len(rows_traminer)} | TraMineRextras rows: {len(rows_extras)}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
