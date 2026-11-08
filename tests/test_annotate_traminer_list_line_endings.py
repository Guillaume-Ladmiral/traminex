"""``annotate_traminer_list.py`` must write LF line endings (audit S7-19, 2026-09-22).

What was wrong: the script wrote ``reports/parity/traminer_API_list_with_traminex.csv``
with ``csv.DictWriter``'s default ``lineterminator="\\r\\n"``, while
``.gitattributes`` declares ``*.csv text eol=lf``; regenerating the tracked file
therefore rewrote every line ending. The sibling generator
``generate_api_mappings.py`` already passes ``lineterminator="\\n"``.

The script resolves its inputs and output from its own location, so it runs on
a copy of exactly the files it reads.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_annotated_parity_list_is_written_with_lf_line_endings(tmp_path: Path) -> None:
    for relative in (
        "scripts/annotate_traminer_list.py",
        "docs/API_PARITY_TABLE.csv",
        "reports/parity/traminer_API_list.csv",
    ):
        target = tmp_path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / relative, target)

    result = subprocess.run(
        [sys.executable, str(tmp_path / "scripts" / "annotate_traminer_list.py")],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    written = (tmp_path / "reports" / "parity" / "traminer_API_list_with_traminex.csv").read_bytes()
    assert written, "the script wrote an empty file"
    assert b"\r\n" not in written
    # Regenerating from the tracked inputs reproduces the tracked output.
    tracked = ROOT / "reports" / "parity" / "traminer_API_list_with_traminex.csv"
    assert written == tracked.read_bytes()
