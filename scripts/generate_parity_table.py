import argparse
import json
from pathlib import Path

# Defaults are anchored to this repository, not the working directory (S7-18).
REPO = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate docs/API_PARITY_TABLE.md from audit JSON."
    )
    parser.add_argument(
        "--audit-json",
        default=str(REPO / "reports" / "audit" / "traminex_api_audit.json"),
        help="Path to the audit JSON input.",
    )
    parser.add_argument(
        "--output",
        default=str(REPO / "docs" / "API_PARITY_TABLE.md"),
        help="Path to the Markdown output file.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    audit_path = Path(args.audit_json)
    out_path = Path(args.output)
    if not audit_path.exists():
        raise SystemExit(f"Audit JSON not found: {audit_path}")
    data = json.loads(audit_path.read_text(encoding="utf-8"))

    exact = data.get("exact_matches", [])
    sanitized = data.get("sanitized_name_matches", [])
    suggested = data.get("suggested_mappings", [])
    missing = data.get("missing", [])

    def md_table(rows, headers):
        # rows: list of tuple/list of strings
        out = []
        out.append("| " + " | ".join(headers) + " |")
        out.append("| " + " | ".join(["---"] * len(headers)) + " |")
        for r in rows:
            out.append("| " + " | ".join(r) + " |")
        return "\n".join(out)

    lines = []
    lines.append("# API Parity Table (Function-by-Function)")
    lines.append("")
    summary = data.get("summary", {})
    if summary:
        lines.append("Summary: " + ", ".join([f"{k}={v}" for k, v in summary.items()]))
        lines.append("")
        # Two artifacts published different R denominators with no reconciliation
        # anywhere (cycle-6 F-008). State it beside the figure it explains.
        total_unique = summary.get("R_total_unique")
        if total_unique is not None:
            lines.append(
                f"`R_total_unique={total_unique}` counts the exported names of both packages "
                "plus every S3 method their `NAMESPACE` files register, de-duplicated. "
                "`docs/API_PARITY_LEDGER.md` reports a larger figure because it additionally "
                "carries one row per verified unexported R internal that traminex promotes to "
                "a public name (`scripts/_r_namespace.py::PROMOTED_INTERNALS`); that "
                "difference is the whole of the gap between the two numbers."
            )
            lines.append("")

    # Exact matches
    lines.append("## Exact name matches")
    if exact:
        rows = sorted([(e["r"], e["py"]) for e in exact], key=lambda t: t[0].lower())
        lines.append(md_table(rows, ["R", "Python"]))
    else:
        lines.append("(none)")
    lines.append("")

    # Sanitized-name matches
    lines.append("## Sanitized-name matches (R name normalized)")
    if sanitized:
        rows = sorted([(e["r"], e["py"]) for e in sanitized], key=lambda t: t[0].lower())
        lines.append(md_table(rows, ["R", "Python"]))
    else:
        lines.append("(none)")
    lines.append("")

    # Suggested mappings (top up to 3 suggestions)
    lines.append("## Suggested mappings (top suggestions by similarity)")
    if suggested:
        sugg_rows = []
        for s in suggested:
            rname = s.get("r", "")
            suggs = s.get("suggestions", [])
            if suggs:
                top = ", ".join([f"{it['py']} ({it['score']:.2f})" for it in suggs[:3]])
            else:
                top = ""
            sugg_rows.append((rname, top))
        sugg_rows.sort(key=lambda t: t[0].lower())
        lines.append(md_table(sugg_rows, ["R", "Suggestions"]))
    else:
        lines.append("(none)")
    lines.append("")

    # Missing
    lines.append("## Missing (no suggestion >= 0.60)")
    if missing:
        miss_rows = sorted([(m, "") for m in missing], key=lambda t: t[0].lower())
        lines.append(md_table(miss_rows, ["R", ""]))
    else:
        lines.append("(none)")
    lines.append("")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
