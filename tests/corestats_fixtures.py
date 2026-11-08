"""Shared A/B/C(+weights) fixtures for the RP3 core-stats R-goldens (T3-T11).

Rebuilds the exact fixtures of ``scripts/r_goldens/make_goldens_corestats.R`` so
the Python side feeds TraMineR the identical numbers:

- ``A``  void-free, no missing (alphabet A,B,C,D)
- ``B``  interior missing ``*`` (void-free; R uses ``right=NA``)
- ``C``  ragged: trailing ``*`` -> void ``%``, interior ``*`` -> nr (R ``right="DEL"``)

The weighted variants attach ``WEIGHTS``. State codes align 0-based with R's
``seqasnum`` because the alphabet is given explicitly.
"""

from __future__ import annotations

import numpy as np

from traminex.format import seqdef

WEIGHTS = np.array([1.5, 2.0, 0.5, 3.0], dtype=float)

ALPH_A = ["A", "B", "C", "D"]
SEQS_A = [
    ["A", "A", "B", "C", "C", "C"],
    ["A", "B", "B", "B", "C", "A"],
    ["D", "D", "D", "D", "D", "D"],
    ["A", "B", "A", "B", "A", "B"],
]
ALPH_B = ["A", "B", "C"]
SEQS_B = [
    ["A", "A", "*", "B", "B", "C"],
    ["*", "A", "A", "A", "*", "*"],
    ["A", "*", "A", "B", "B", "B"],
    ["C", "C", "C", "C", "C", "C"],
]
# Fixture C: trailing "*" written directly as the void "%" (matches R right="DEL"),
# interior "*" kept as nr.
ALPH_C = ["A", "B", "C"]
SEQS_C = [
    ["A", "A", "B", "C", "C", "C"],
    ["A", "B", "B", "%", "%", "%"],
    ["A", "*", "A", "B", "%", "%"],
    ["C", "C", "%", "%", "%", "%"],
]


def fixture_a(weighted: bool = False):
    w = WEIGHTS if weighted else None
    return seqdef(np.array(SEQS_A, dtype=object), alphabet=ALPH_A, weights=w)


def fixture_b(weighted: bool = False):
    w = WEIGHTS if weighted else None
    return seqdef(np.array(SEQS_B, dtype=object), alphabet=ALPH_B, nr="*", weights=w)


def fixture_c(weighted: bool = False):
    w = WEIGHTS if weighted else None
    return seqdef(np.array(SEQS_C, dtype=object), alphabet=ALPH_C, nr="*", weights=w)
