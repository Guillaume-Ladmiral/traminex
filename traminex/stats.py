"""Facade module for sequence statistics.

Implementation is split across:
- stats_basic.py: core counts/distributions
- stats_complexity.py: turbulence and complexity indices
- stats_utils.py: helper and indicator utilities
- stats_p1.py: seqmeant (also routed at top level by traminex/__init__.py)
- stats_p2.py: seqici, seqpm, seqsep (also routed at top level)
"""

from __future__ import annotations

from .stats_basic import entropy, seqistatd, seqlength, seqstatd, seqtab, seqtransn
from .stats_complexity import (
    seqdss,
    seqdur,
    seqibad,
    seqidegrad,
    seqinsecurity,
    seqintegr,
    seqipos,
    seqivardur,
    seqivolatility,
    seqlogp,
    seqmaintokens,
    seqprecarity,
    seqprecstart,
    seqST,
    seqtransn_count,
    turbulence_exact,
    turbulence_proxy,
)
from .stats_p1 import seqmeant
from .stats_p2 import seqici, seqpm, seqsep
from .stats_utils import (
    group_p,
    rowmode,
    seqasnum,
    seqhasmiss,
    seqindic,
    seqindic_dyn,
    seqnum,
    seqtabstocc,
    seqtrate,
    sorti,
    sortv,
)

__all__ = [
    "seqlength",
    "seqistatd",
    "seqstatd",
    "seqtab",
    "seqtransn",
    "entropy",
    "seqtransn_count",
    "turbulence_proxy",
    "seqST",
    "seqdss",
    "seqdur",
    "seqivolatility",
    "turbulence_exact",
    "seqmeant",
    "seqivardur",
    "seqlogp",
    "seqmaintokens",
    "seqpm",
    "seqprecarity",
    "seqprecstart",
    "seqidegrad",
    "seqinsecurity",
    "seqibad",
    "seqici",
    "seqsep",
    "seqintegr",
    "seqipos",
    "group_p",
    "sorti",
    "sortv",
    "rowmode",
    "seqtabstocc",
    "seqasnum",
    "seqindic",
    "seqindic_dyn",
    "seqnum",
    "seqtrate",
    "seqhasmiss",
]
