#!/usr/bin/env python3
import csv
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parity = root / "docs" / "API_PARITY_TABLE.csv"
parity_dir = root / "reports" / "parity"
traminer_list = parity_dir / "traminer_API_list.csv"
outp = parity_dir / "traminer_API_list_with_traminex.csv"

parity_dir.mkdir(parents=True, exist_ok=True)

pref = {"exact": 0, "sanitized": 1, "suggested": 2, "missing": 3}
mp = {}
with parity.open("r", encoding="utf-8", newline="") as f:
    rdr = csv.DictReader(f)
    for row in rdr:
        r = (row.get("r_name", "") or "").strip()
        py = (row.get("py_name", "") or "").strip()
        mt = (row.get("match_type", "") or "").strip()
        sc = (row.get("score", "") or "").strip()
        if not r:
            continue
        rank = pref.get(mt, 99)
        cur = mp.get(r)
        if cur is None or rank < cur[0]:
            mp[r] = (rank, py, mt, sc)

rows = []
with traminer_list.open("r", encoding="utf-8", newline="") as f:
    rdr = csv.DictReader(f)
    flds = list(rdr.fieldnames or [])
    if "traminer_API" in flds:
        rcol = "traminer_API"
    elif "traminer_api" in flds:
        rcol = "traminer_api"
    elif "r_name" in flds:
        rcol = "r_name"
    else:
        rcol = flds[0] if flds else "name"
    for row in rdr:
        rname = (row.get(rcol, "") or "").strip()
        best = mp.get(rname)
        if best:
            _, py, mt, sc = best
        else:
            py, mt, sc = "", "missing", ""
        row["traminex_API"] = py
        row["match_type"] = mt
        row["score"] = sc
        rows.append(row)

with outp.open("w", encoding="utf-8", newline="") as f:
    fieldnames = rows[0].keys() if rows else ["traminer_API", "traminex_API", "match_type", "score"]
    # .gitattributes declares *.csv eol=lf; csv's default is CRLF (audit S7-19).
    w = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow(r)

print(outp)
