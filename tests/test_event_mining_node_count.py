"""Event mining's stop rule counts the nodes of its own tree (TX-AUD-20261011-021).

The miner stops when a pass adds no node to the prefix tree (R's
``nodeCount - lastNodeCount == 0``). traminex mirrored R's static C++ counter
with a class attribute, so another mining call creating or freeing nodes at the
same time changed the count and could stop a pass early: a legal interleaving
returned 6 subsequences instead of 7. Thread B below holds and frees nodes of
another tree exactly between thread A's second pass and its stop-rule check.
"""

from __future__ import annotations

import threading

import pandas as pd

import traminex as tx
from traminex.events import _engine

EVENTS = pd.DataFrame(
    {
        "id": [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4],
        "time": [0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1],
        "event": ["A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B"],
    }
)
SERIAL = ["(A)", "(A)-(B)", "(B)", "(A)-(B)-(C)", "(A)-(C)", "(B)-(C)", "(C)"]


def _eseq():
    return tx.seqecreate(
        id=EVENTS["id"].tolist(), timestamp=EVENTS["time"].tolist(), event=EVENTS["event"].tolist()
    )


def _foreign_nodes(count: int) -> list:
    """Nodes of another prefix tree, as a concurrent mining call holds them."""
    tree = _engine._PrefixTree()
    counter = getattr(tree, "counter", None)  # the per-tree counter of the repair
    if counter is None:
        return [_engine._TreeEventNode(1) for _ in range(count)]
    return [_engine._TreeEventNode(1, counter) for _ in range(count)]


def test_serial_mining_result() -> None:
    assert list(tx.seqefsub(_eseq(), pmin_support=0.5).subseq) == SERIAL


def test_concurrent_mining_does_not_change_the_stop_rule(monkeypatch) -> None:
    eseq = _eseq()
    held, release, released = threading.Event(), threading.Event(), threading.Event()
    original = _engine._PrefixTree.simplify_tree
    passes = {"A": 0}

    def thread_b() -> None:
        nodes = _foreign_nodes(3)  # the k=2 pass adds 3 nodes serially
        held.set()
        release.wait(timeout=30)
        for node in nodes:
            _engine._free_subtree(node)
        released.set()

    def interleaved(self, min_sup):
        out = original(self, min_sup)
        if threading.current_thread().name == "A":
            passes["A"] += 1
            if passes["A"] == 2:
                release.set()
                released.wait(timeout=30)
        return out

    monkeypatch.setattr(_engine._PrefixTree, "simplify_tree", interleaved)
    result: dict[str, list[str]] = {}

    def thread_a() -> None:
        held.wait(timeout=30)
        result["subseq"] = list(tx.seqefsub(eseq, pmin_support=0.5).subseq)

    b = threading.Thread(target=thread_b, name="B")
    a = threading.Thread(target=thread_a, name="A")
    b.start()
    a.start()
    a.join(timeout=60)
    b.join(timeout=60)

    assert result["subseq"] == SERIAL
