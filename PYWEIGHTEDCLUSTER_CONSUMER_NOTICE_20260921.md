# Notice: two provider-tree states that stop `pyweightedcluster` from importing

**Written:** 2026-09-21 (Asia/Tokyo)

**Author:** `pyweightedcluster` session, working at that project's commit
`149ae71826bc126aaa82763fa79a4e6477ba03b5`

**Audience:** `traminex` maintainers

**Status:** informational. Nothing in this repository was modified to produce
this notice, and no change is requested to traminex's own behaviour. One
sibling action was taken and is disclosed in §4; a later merge of the 2.2-14
branch and its reversal are disclosed in §0, which also records a defect that
**does** require a change before that branch is merged.

> **Status note (2026-10-24, added by the traminex repair).** This notice is a
> dated record of 2026-09-21; its text is kept. The re-pin digest it points to
> has been superseded by later traminex commits, and re-pinning the consumer is
> deferred by the owner (TX-AUD-20261010-099). Its line numbers in traminex
> files refer to the tree it read on 2026-09-21. Its "79 Python sources" are
> those of the consumer's pinned commitment, which still holds 79
> (TX-AUD-20261011-030).

## 0. CORRECTION -- read this first (added later on 2026-09-21)

> **RESOLVED 2026-09-22 (added by the traminex audit-and-fix pass).** The
> defect below is fixed by `e85beae2`, and `oracle-2214-remediation-20260920`
> is now merged into `master` (`6a5b299b`). The consumer's own call,
> `seqdist(X, method="OM", subm=s, with_missing=True)` for `CONSTANT`,
> `INDELSLOG` and `TRATE`, was re-run with `probe.py` under `6e76440b` and under
> the fix: 36/36 matrices identical, 0 errors (320/320 bit-identical in an
> unrounded variant). **The 2.2-14 work has now landed on `master`**, which is
> the signal §2 asks for; the final `master` revision, its `.py` inventory and
> digest for re-pinning are in `reports/audit/full_codebase_audit_20260922.md`.
> The text of this section below is the original, unchanged.

**§3 below is wrong in two places, and its conclusion must not be relied on.
Do not merge `oracle-2214-remediation-20260920` until the defect below is fixed.**

1. §3's `distances/costs.py` row says the consumer "never calls `seqcost`".
   False: `pyweightedcluster.seqnullcqi` forwards its `seqdist_args` to public
   `traminex.seqdist`, which calls `seqcost` with `with_missing=True` for `sm`
   values such as `CONSTANT` and `INDELSLOG`.
2. §3's finding that "no numerical change in this delta reaches the consumer's
   delegated surface" is refuted. At the `pyweightedcluster` maintainer's
   instruction the branch was merged into `master` (`c84b8d6a`). Before
   re-pinning, the call the consumer makes was executed under the old revision
   and the merge:

   | `subm` in `seqdist(X, method="OM", subm=..., with_missing=True)` | `6e76440b` | merge `c84b8d6a` |
   |---|---|---|
   | `CONSTANT` | ok | ok, identical |
   | `TRATE` | ok | ok, identical |
   | `INDELSLOG` | ok | `ValueError: OM indel vector length must equal the active state count` |

   36 matrices over 12 seeded inputs without missing values: 24 identical,
   0 differing in value, 12 errors -- every `INDELSLOG` case. **Live R 2.2.14
   returns a result for the same sequences with `with.missing=TRUE` and
   `FALSE`**, so the merge diverges from the oracle it adopts. The branch's own
   suite passed (3943 passed, 0 failed), because nothing in it exercises this
   call. The old revision was loaded with `python -S` from the
   `audit-cycle10-20260919` worktree, so the editable install could not shadow
   it, and `traminex.__file__` was checked on both sides.

   *Probable* cause, not yet confirmed in traminex's source: R's `seqdist`
   resets `with.missing` to FALSE when neither the data nor a sequence
   reference carries a missing value (`seqdist.R:111-116`, identical in both
   releases) before computing costs. traminex's `seqdist` does not, and
   2.2-13's `seqcost` downgrade was masking that. With the downgrade removed,
   `INDELSLOG`'s indel *vector* gains the missing state and no longer matches
   the active state count the OM kernel validates.

   **Superseded root cause.** The "probable" cause above is partly wrong:
   traminex's `seqdist` *does* reset `with_missing`, but only on its
   sequence-object path; a plain array skips the reset. Confirmed, with R's
   expected matrix and acceptance criteria, in
   `ISSUE_SEQDIST_WITH_MISSING_RAW_ARRAY_20260921.md`.

**The merge was undone** (`git reset --keep ce67776c`). `master` is back at
`ce67776c`, and `oracle-2214-remediation-20260920` still holds every commit.
`pyweightedcluster` was never re-pinned; its pin remains `6e76440b`
(`d6eafe73...`), and its strict suite passes against the restored `master`
(2369 passed, 0 failed, 0 skipped). Before a future merge, re-run the comparison
above and require zero errors. Reproduction evidence is retained in the consumer
repository at `audit/traminex-repin-20260921/merge_regression/`.

## Summary

`pyweightedcluster` verifies the installed traminex provider on **every**
import, against a packaged commitment of 79 Python sources for a reviewed
revision. The check is fail-closed and has no environment-variable bypass.

Between 2026-09-20 and 2026-09-21 the provider tree entered two different
states that each stopped the consumer from importing at all. Neither was a
defect in traminex's code, and neither is visible from inside traminex —
both are only observable downstream. They are recorded here so the build
and working-tree habits that produce them can be considered.

A third section reports, as a courtesy, the consumer-side review of the
unmerged `oracle-2214-remediation-20260920` work.

## 1. A build symlink in the package tree bricks the consumer

**Observed:** 2026-09-20 18:28 JST.

A build placed a symbolic link into the package tree:

```
traminex/_pb_kernels.cpython-312-x86_64-linux-gnu.so
  -> /tmp/traminex-r1-pybind11.mpaveQ/_pb_kernels.cpython-312-x86_64-linux-gnu.so
```

The consumer's verifier walks the whole package tree and raises on **any**
symlink entry it encounters, before deciding whether the entry is one of the
79 files it commits to (`_traminex_compatibility.py:258`). Every
`import weightedcluster` then failed with:

```
TraminexCompatibilityError: traminex package contains a symlink:
_pb_kernels.cpython-312-x86_64-linux-gnu.so
```

**Consumer impact, measured.** `pytest -q` went from its recorded baseline of
3 failures to **15 failures / 1899 passed / 444 xfailed**. The twelve extra
failures were all tests that spawn a fresh interpreter and re-import the
package; the error appears 21 times in that run's log. After the tree was
restored the same suite returned **3 failed / 1911 passed / 444 xfailed**,
confirming all twelve were environmental.

**Note the asymmetry.** This file is gitignored here (`.gitignore:19`), is a
compiled artifact rather than a source, and is **not** one of the 79 files the
commitment covers. traminex itself is unaffected either way, because
`traminex/distances/__init__.py:56-60` imports `_pb_kernels` inside a
`try/except ImportError` and falls back to the Python backend. So a file that
neither project considers part of its contract nevertheless stopped the
consumer completely.

There is also a documented tension worth flagging: the consumer's own reviewed
manifest procedure says to "exclude caches, compiled artifacts, and
distribution metadata", yet its symlink refusal still applies to exactly such
an artifact. That is the consumer's issue to reconcile, not traminex's, and it
is being recorded there too.

**Suggested mitigations** (traminex's call entirely):

- have the pybind11 build emit into a directory outside `traminex/`, or
- copy the built extension into place instead of linking it, or
- keep the link but accept that consumers cannot import while it exists.

## 2. Leaving the working tree on an unmerged branch also bricks the consumer

**Observed:** 2026-09-21, after the symlink was gone.

The working tree was checked out on `oracle-2214-remediation-20260920`
(`acc75529`). The commitment names `6e76440b`, so the consumer failed with a
digest mismatch instead:

```
expected d6eafe7300d36493d0f814692fe25e51fd5170036eea827305200e12832e5443
found    d73768785139e154e1a56a00622e7341b01069ae4ad73fa74138b7e1b3654519
```

Computed from git objects, without touching this working tree:

| Revision | Tree digest | Matches the consumer's pin |
|---|---|---|
| `master` (`6e76440b`) | `d6eafe73…` | **yes, exactly** |
| `acc75529` (branch tip) | `d7376878…` | no |

The file **set** is identical (79 in both); only contents differ.

The consumer's pin is therefore correct and current with respect to `master`.
The consumer deliberately did **not** re-pin to `acc75529`, because committing
a downstream project to an unmerged branch tip would leave a dangling
commitment if that branch is rebased, squashed on merge, or abandoned.

**Request:** when the 2.2-14 work lands on `master`, please say so. The
consumer will then re-pin through its reviewed procedure against the merged
revision.

## 3. Consumer-side review of the 2.2-14 work (courtesy; its conclusion is refuted -- see §0)

The delta `6e76440b..acc75529` was reviewed against the surface
`pyweightedcluster` actually delegates to — `distances.pairwise`,
`distances.pairwise_between_sets`, public `seqdist`, `diss.*`, `dissmfacw`,
`disstree.disstree`, `seqtree.seqtree`, `format.*`, `seqformat`.

5 commits, 25 changed `.py` files, file set unchanged.

| Area | Change | Reaches the consumer |
|---|---|---|
| `diss.py`, `disstree.py`, `seqtree.py`, `distances/_normalization.py` | docstrings 2.2-13 → 2.2-14 only | no |
| `distances/costs.py` | real: `seqcost` no longer downgrades `with_missing`; warning removed | **no** *(refuted -- see §0)* — the consumer never calls `seqcost`; its own `with_missing` is `seqpropclust`'s argument feeding `seqdss`/`seqdur`/`seqistatd` |
| `distances/__init__.py` | real: `om_trate_pairwise` automatic indel `0.5*mean(offdiag)` → `0.5*max(sm)`; new guard for an undersized reference `subm` | **no** — the consumer uses `pairwise`/`seqdist`, not `om_trate_pairwise` |
| `distances/_validation.py`, `distances/python_kernels.py` | earlier `DUR`/`eps` refusal; suppressed a spurious NumPy warning | message/warning only |

**Finding (refuted -- see §0): no numerical change in this delta reaches the
consumer's delegated surface.** The removed `seqcost` warning has no consumer contract either —
zero references to its text anywhere in `pyweightedcluster`.

This is a static review, not a verification. It does not replace re-running the
consumer's suite after the eventual re-pin, and it should not be cited as
evidence that the re-pin is safe.

## 4. Sibling action taken, disclosed

At the `pyweightedcluster` maintainer's explicit instruction, and only after
confirming the tracked tree was clean and nothing was stashed, this session ran
`git checkout master` in this repository to restore a pinned-compatible state.

No file was created, modified or deleted here except this notice.
`oracle-2214-remediation-20260920` retains every commit, including `acc75529`.

## 5. A note on the oracle version

`pyweightedcluster` still pins TraMineR **2.2-13** as its companion oracle,
while this project has moved to **2.2-14** and the machine's R library has
2.2.14 installed. That divergence is the consumer's to resolve and is being
acted on separately; it is mentioned only so the two projects' oracle
statements are not read as contradictory in the interim.
