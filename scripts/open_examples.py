#!/usr/bin/env python3
"""
Open example notebooks for traminex.

Usage:
  python scripts/open_examples.py [--lab]

Tries to launch Jupyter (Notebook by default). If --lab is given, tries jupyter lab.
If Jupyter is not installed, prints paths to the notebooks.
"""

import argparse
import os
import subprocess
import sys
from shutil import which

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

NOTEBOOKS = [
    os.path.join(ROOT, "examples", "disstree_example.ipynb"),
    os.path.join(ROOT, "examples", "hspell_overlap_example.ipynb"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lab", action="store_true", help="Use JupyterLab instead of classic notebook")
    args = ap.parse_args()

    cmd = "jupyter lab" if args.lab else "jupyter notebook"
    exe = cmd.split()[0]
    if which(exe) is None:
        print("Jupyter does not seem to be installed.")
        print("You can install it with: pip install notebook OR pip install jupyterlab")
        print("Example notebooks are located at:")
        for nb in NOTEBOOKS:
            print(" -", nb)
        sys.exit(0)
    try:
        subprocess.Popen(cmd.split() + NOTEBOOKS, cwd=ROOT)
        print("Launched:", cmd, *NOTEBOOKS)
    except Exception as e:
        print("Failed to launch Jupyter:", e, file=sys.stderr)
        print("Example notebooks are located at:")
        for nb in NOTEBOOKS:
            print(" -", nb)
        sys.exit(1)


if __name__ == "__main__":
    main()
