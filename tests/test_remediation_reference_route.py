"""Exercise optional reference consumers with an independent pinned R oracle.

TX-AUD-20261005-020/-021/-027: incomplete references must skip, complete
comparisons must execute, and calibration output names must all be usable.
Expected CSVs are produced by TraMineR, never by the Python implementation.
"""

import json
import os
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
FREQUENCIES = "expected_seqstatd_freq.csv"
DISTANCES = "expected_dist_om.csv"

R_ORACLE = r"""
suppressPackageStartupMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
dest <- commandArgs(TRUE)[1]
dir.create(dest, recursive=TRUE, showWarnings=FALSE)
df <- data.frame(T1=c("A","A","B"), T2=c("A","B","A"), T3=c("B","B","A"))
write.csv(df, file.path(dest,"input_wide.csv"), row.names=FALSE)
alphabet <- unique(as.vector(t(as.matrix(df))))
x <- seqdef(df, alphabet=alphabet)
write.table(t(seqstatd(x)$Frequencies), file.path(dest,"expected_seqstatd_freq.csv"),
            sep=",", row.names=FALSE, col.names=FALSE)
sm <- seqsubm(x, method="CONSTANT", cval=2)
d <- seqdist(x[1:2,], method="OM", indel=1, sm=sm)
write.table(as.matrix(d), file.path(dest,"expected_dist_om.csv"),
            sep=",", row.names=FALSE, col.names=FALSE)
trate <- seqsubm(x, method="TRATE")
write.table(trate, file.path(dest,"expected_trate_subm.csv"),
            sep=",", row.names=FALSE, col.names=FALSE)
dhd <- seqdist(x[1:2,], method="DHD")
write.table(as.matrix(dhd), file.path(dest,"expected_dhd.csv"),
            sep=",", row.names=FALSE, col.names=FALSE)
cat("oracle version", as.character(packageVersion("TraMineR")), "\n")
"""


@pytest.fixture(scope="module")
def r_references(tmp_path_factory):
    work = tmp_path_factory.mktemp("reference-oracle")
    script = work / "oracle.R"
    script.write_text(R_ORACLE, encoding="utf-8")
    refs = work / "references"
    run = subprocess.run(
        ["Rscript", str(script), str(refs)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 0, run.stdout + run.stderr
    assert "oracle version 2.2.14" in run.stdout
    return refs


def run_optional_check(work, oracle, filenames):
    """Run the actual maintained test, with exactly the requested CSVs present."""
    test_dir = work / "tests"
    refs = test_dir / "references"
    refs.mkdir(parents=True)
    shutil.copyfile(ROOT / "tests/test_parity_refs.py", test_dir / "test_parity_refs.py")
    for filename in filenames:
        shutil.copyfile(oracle / filename, refs / filename)
    report = work / "results.xml"
    run = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "-rA",
            "tests/test_parity_refs.py",
            f"--junitxml={report}",
        ],
        cwd=work,
        env={**os.environ, "PYTHONPATH": str(ROOT)},
        capture_output=True,
        text=True,
        check=False,
    )
    cases = ET.parse(report).findall(".//testcase")
    assert len(cases) == 1, run.stdout + run.stderr
    return run, cases[0]


@pytest.mark.parametrize(
    "filenames",
    [[], ["input_wide.csv"], [FREQUENCIES]],
    ids=["no-input", "input-only", "expectation-only"],
)
def test_incomplete_optional_references_report_skip(tmp_path, r_references, filenames):
    run, case = run_optional_check(tmp_path, r_references, filenames)
    assert run.returncode == 0, run.stdout + run.stderr
    assert case.find("skipped") is not None, run.stdout + run.stderr


@pytest.mark.parametrize("expectations", [[FREQUENCIES], [DISTANCES], [FREQUENCIES, DISTANCES]])
def test_complete_independent_r_comparison_executes(tmp_path, r_references, expectations):
    run, case = run_optional_check(tmp_path, r_references, ["input_wide.csv", *expectations])
    assert run.returncode == 0, run.stdout + run.stderr
    assert case.find("skipped") is None, run.stdout + run.stderr


def test_partial_reference_coverage_is_reported(tmp_path, r_references):
    run, case = run_optional_check(tmp_path, r_references, ["input_wide.csv", FREQUENCIES])
    assert run.returncode == 0, run.stdout + run.stderr
    assert case.find("skipped") is None
    output = run.stdout + run.stderr
    assert "partial" in output.lower(), output
    assert DISTANCES in output, output


def test_wrong_distance_fixture_fails_a_real_comparison(tmp_path, r_references):
    # OM(AAB, ABB) = 2 in the R output; deliberate corruption must fail.
    broken = tmp_path / "broken-oracle"
    broken.mkdir()
    shutil.copyfile(r_references / "input_wide.csv", broken / "input_wide.csv")
    (broken / DISTANCES).write_text("0,99\n99,0\n", encoding="utf-8")
    run, case = run_optional_check(tmp_path / "consumer", broken, ["input_wide.csv", DISTANCES])
    assert run.returncode == 1, run.stdout + run.stderr
    assert case.find("failure") is not None, run.stdout + run.stderr


@pytest.mark.parametrize("output_kind", ["basename", "nested", "absolute"])
def test_calibration_creates_requested_output_path(tmp_path, r_references, output_kind):
    refs = tmp_path / "refs"
    refs.mkdir()
    for filename in ("input_wide.csv", "expected_trate_subm.csv"):
        shutil.copyfile(r_references / filename, refs / filename)
    destination = {
        "basename": Path("calibrated.json"),
        "nested": Path("results") / "nested" / "calibrated.json",
        "absolute": tmp_path / "absolute" / "calibrated.json",
    }[output_kind]
    run = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/calibrate_parity.py"),
            "--refs",
            str(refs),
            "--out",
            str(destination),
        ],
        cwd=tmp_path,
        env={**os.environ, "PYTHONPATH": str(ROOT)},
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 0, run.stdout + run.stderr
    artifact = destination if destination.is_absolute() else tmp_path / destination
    config = json.loads(artifact.read_text(encoding="utf-8"))
    assert set(config) == {"trate"}
    assert set(config["trate"]) == {"smoothing", "mode", "normalize"}
    assert config["trate"]["mode"] in {"symmetric", "row", "col"}
    assert config["trate"]["normalize"] in {"none", "mean_offdiag"}
    assert config["trate"]["smoothing"] >= 0


# TX-AUD-20261019-087: a calibration records only parameters that reproduce the R
# reference within the tolerance tests/test_parity_refs.py uses (TRATE 1e-6, DHD 2e-3).


def _calibration_sequences(refs):
    import pandas as pd

    from traminex import seqdef

    return seqdef(pd.read_csv(refs / "input_wide.csv"))


def _run_calibration(tmp_path, refs, destination):
    return subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/calibrate_parity.py"),
            "--refs",
            str(refs),
            "--out",
            str(destination),
        ],
        cwd=tmp_path,
        env={**os.environ, "PYTHONPATH": str(ROOT)},
        capture_output=True,
        text=True,
        check=False,
    )


def test_dhd_calibration_records_nothing_r_does_not_reproduce(r_references):
    # R's default DHD uses transition costs; no frequency variant reproduces it.
    from scripts import calibrate_parity

    x = _calibration_sequences(r_references)

    assert calibrate_parity.calibrate_dhd(x, r_references / "expected_dhd.csv") is None


def test_trate_calibration_keeps_a_candidate_that_reproduces_r(r_references):
    from scripts import calibrate_parity
    from traminex.distances.costs import trate_substitution_matrix

    x = _calibration_sequences(r_references)
    params = calibrate_parity.calibrate_trate(x, r_references / "expected_trate_subm.csv")

    reference = np.loadtxt(r_references / "expected_trate_subm.csv", delimiter=",")
    assert params is not None
    assert np.allclose(trate_substitution_matrix(x.data, **params), reference, atol=1e-6)


def test_trate_calibration_records_nothing_when_no_candidate_reproduces_r(tmp_path, r_references):
    from scripts import calibrate_parity

    reference = np.loadtxt(r_references / "expected_trate_subm.csv", delimiter=",")
    shifted = tmp_path / "expected_trate_subm.csv"
    np.savetxt(shifted, reference + 0.01, delimiter=",")
    x = _calibration_sequences(r_references)

    assert calibrate_parity.calibrate_trate(x, shifted) is None


def test_trate_calibration_holds_an_absolute_tolerance(tmp_path, r_references):
    # Batch-17 review M-1: np.allclose's default relative tolerance let a TRATE cost near 2
    # pass up to about 2e-5 off, where the stated tolerance is an absolute 1e-6.
    from scripts import calibrate_parity

    reference = np.loadtxt(r_references / "expected_trate_subm.csv", delimiter=",")
    shifted_values = reference + 5e-6 * ~np.eye(len(reference), dtype=bool)
    assert np.allclose(shifted_values, reference, atol=1e-6)
    assert not np.allclose(shifted_values, reference, rtol=0.0, atol=1e-6)
    shifted = tmp_path / "expected_trate_subm.csv"
    np.savetxt(shifted, shifted_values, delimiter=",")
    x = _calibration_sequences(r_references)

    assert calibrate_parity.calibrate_trate(x, shifted) is None


def test_calibration_writes_no_dhd_section_and_names_the_public_dhd(tmp_path, r_references):
    refs = tmp_path / "refs"
    refs.mkdir()
    for filename in ("input_wide.csv", "expected_trate_subm.csv", "expected_dhd.csv"):
        shutil.copyfile(r_references / filename, refs / filename)
    destination = tmp_path / "calibrated.json"

    run = _run_calibration(tmp_path, refs, destination)

    assert run.returncode == 0, run.stdout + run.stderr
    assert set(json.loads(destination.read_text(encoding="utf-8"))) == {"trate"}
    assert 'seqdist(method="DHD")' in run.stdout


def test_calibration_replaces_an_earlier_file_when_nothing_reproduces_r(tmp_path, r_references):
    # Batch-17 review M-2: an earlier file kept its sections, here the DHD variant the old
    # calibration recorded, and a calibrated profile went on reading them.
    import traminex as tx

    refs = tmp_path / "refs"
    refs.mkdir()
    for filename in ("input_wide.csv", "expected_dhd.csv"):
        shutil.copyfile(r_references / filename, refs / filename)
    destination = tmp_path / "calibrated.json"
    destination.write_text('{"dhd": {"variant": "sum", "smooth_window": 5}}', encoding="utf-8")

    run = _run_calibration(tmp_path, refs, destination)

    assert run.returncode == 0, run.stdout + run.stderr
    assert json.loads(destination.read_text(encoding="utf-8")) == {}
    assert "empty calibration" in run.stdout
    sequences = _calibration_sequences(refs)
    with pytest.raises(ValueError):
        tx.seqdist(
            sequences, method="DHD", dhd_profile="calibrated", calibrated_params_path=destination
        )


def test_full_independent_r_reference_recipe_executes(tmp_path, r_references):
    run, case = run_optional_check(
        tmp_path,
        r_references,
        ["input_wide.csv", FREQUENCIES, DISTANCES, "expected_trate_subm.csv", "expected_dhd.csv"],
    )
    assert run.returncode == 0, run.stdout + run.stderr
    assert case.find("skipped") is None, run.stdout + run.stderr
    assert "partial reference coverage" not in run.stdout + run.stderr
