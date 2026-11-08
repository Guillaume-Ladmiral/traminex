"""Event-sequence framework.

This package keeps the public import surface of the former ``traminex.events``
module unchanged (``from traminex.events import …`` and ``from traminex import
seqe…`` both keep working). It is split into:

- :mod:`traminex.events._object` -- the **R-faithful** (RP5a) eseq/seqelist
  object, accessors and transition machinery (golden-pinned to TraMineR 2.2-14);
- :mod:`traminex.events._engine` -- the internal pure-Python port of TraMineR's
  vendored C++ event engine (RP5b);
- :mod:`traminex.events._mining` -- the **R-faithful** (RP5b) event-mining surface
  built on the engine (``seqeconstraint`` / ``seqefsub`` / ``seqeapplysub`` /
  ``seqerules`` / ``seqecmpgroup`` / ``seqedist`` / ``seqentrans``);
- :mod:`traminex.events._helpers` -- traminex-only convenience helpers with no R
  equivalent.
"""

from __future__ import annotations

from ._helpers import (
    EventSequence,
    events_to_state_sequences,
    events_transition_degrees,
    events_transition_matrix,
    frequent_event_pairs,
    frequent_event_pairs_timegap,
    frequent_events,
    frequent_ngrams,
)
from ._mining import (
    SeqEConstraint,
    SubSeqElist,
    SubSeqElistChisq,
    is_subseqelist,
    seqeapplysub,
    seqecmpgroup,
    seqeconstraint,
    seqedist,
    seqefsub,
    seqentrans,
    seqerules,
)
from ._object import (
    Eseq,
    SeqElist,
    is_eseq,
    is_seqelist,
    is_stslist,
    seqe2stm,
    seqeage,
    seqecontain,
    seqecreate,
    seqeid,
    seqelength,
    seqetm,
    seqeweight,
)

__all__ = [
    # R-faithful object / accessors / transition machinery (RP5a)
    "Eseq",
    "SeqElist",
    "seqecreate",
    "seqeid",
    "seqeweight",
    "seqelength",
    "seqecontain",
    "is_eseq",
    "is_seqelist",
    "is_stslist",
    "seqetm",
    "seqe2stm",
    "seqeage",
    "seqentrans",
    # R-faithful event-mining engine (RP5b)
    "SeqEConstraint",
    "seqeconstraint",
    "seqefsub",
    "seqeapplysub",
    "seqerules",
    "seqecmpgroup",
    "seqedist",
    "SubSeqElist",
    "SubSeqElistChisq",
    "is_subseqelist",
    # Non-R convenience helpers
    "EventSequence",
    "events_transition_matrix",
    "events_transition_degrees",
    "events_to_state_sequences",
    "frequent_events",
    "frequent_ngrams",
    "frequent_event_pairs",
    "frequent_event_pairs_timegap",
]
