# Owner decisions about traminex

The product decisions the owner made about traminex, gathered in one place. The texts below are copied word
for word from their sources, which are named. Family-wide decisions (layout, standby, records) are in the
parent repository's `docs/DECISIONS.md` (`PROJECT_SEQANALYSIS_and_SPATSCAN`).

## From ECOSYSTEM.md

Owner decision: 2026-10-12 — **traminex corrects defects of its R sources instead of porting them.** Where TraMineR 2.2-14 or TraMineRextras 0.6.9 behaves wrongly (a result that contradicts the method's definition or R's own documentation, an R error or quirk on valid input), traminex implements the correct behaviour. Each such correction is a deliberate change in the sense of "Source of truth" below: it is recorded as a `DIVERGENCE` row of `docs/API_PARITY_DECISIONS.csv` that states R's behaviour, the correction and the reason, and it is tested. R remains the source of truth for all other behaviour. weightedcluster and flexscanx stay on standby until traminex reaches a stable, release-ready version; this decision authorizes no sibling edit, synchronization or migration, and does not resume R6.

Owner decision: 2026-10-19 — **traminex's default alphabet is R's order.** Every traminex constructor that derives a state alphabet without `alphabet=` (`seqdef`, `seqformat`, `seqdef_from_long`, `TSE_to_STS`, `HSPELL_to_STS`, `events_to_state_sequences`, `seqMD`'s combined alphabet) orders the real states as TraMineR's `seqdef` does through `seqstatl`: by numeric value when every value (a `missing=` code included) is an R number, R's ICU root collation otherwise (ICU stays a hard requirement, decision of 2026-10-12). An explicit `alphabet=` keeps the caller's order, as in R. This reverses the earlier first-appearance default, so per-state results, colours by position, tie rules and state-rank defaults follow R without an explicit alphabet.

## From the repair run's report

Section 3 of the `REMEDIATION_REPORT.md` of repair run `repair-20261012T025723Z` (tracked until 2026-10-26;
it stays in git history, last present at commit d2ff400c, and a copy is in the family supervision
repository's `docs/history/run-reports/`). The entries it names are in that run's `REPAIR_STATE.md`; the
finding identifiers are in `docs/FINDINGS_REGISTER.csv`.

From the prompts' section 2 and REPAIR_STATE.md (the entries named):

- 2026-10-12: R defects are corrected, not ported (ECOSYSTEM.md); a correction is tested against the method's
  definition and recorded in the R name's ledger row as DIVERGENCE; everywhere else R is the reference, proved by an R
  golden. R's checkcost diagnostics are reproduced on generated cost matrices too.
- OQ-1 (2026-10-12): seqdef keeps ICU as a hard requirement when two or more text states have no alphabet=, without a
  fallback (no code change; the code already did so).
- OQ-2 and OQ-3 (2026-10-18): seqdef keeps refusing left/gaps/right = "NEUTRAL"; HSPELL_to_STS keeps blanking the
  periods from endObs on (endObs is the first unobserved period). Both were the code's behaviour; the wording was
  corrected (70072fcb).
- D1-D4 (2026-10-19, the second prompt): D1 the run continues on the same branch and records; D2 the default alphabet
  follows R's order; D3 seqpolyads keeps the row-wise pair-weight order its Rd states, as a correction of R's
  column-wise code; D4 seqrecode(otherwise=) applies to real states only. Release strategy (CI jobs, release gates,
  release metadata) is pending with the owner: none was added beyond a finding's own repair.
- R9-g, R10-q, R10-e (2026-10-22): seqdist and seqMD keep accepting the lower-case and alias spellings R refuses, as a
  recorded convenience (option B); gmean of OMloc, OMspell and TWED keeps R's bound while maxdist and YujianBo keep the
  corrected one (option A); the automatic label order follows the loaded ICU, as R's does, with alphabet= the
  documented way to an R-identical order under any ICU version, and the test output shows the loaded ICU (option E).
- 2026-10-23 (TX-AUD-20261019-033): seqrecode's otherwise leaves missing and void cells as markers.
- 2026-10-24 (batch 17): the DATA/ copies of ex1, ex2 and famform are replaced by the vendored bytes and tested; the
  flat package layout stays and is documented (the src/ layout ESCALATED); graphify-out/ stays as a dated snapshot;
  the release-strategy items (CI R pin, CI jobs, build identity) are kept for that discussion. The owner handles the
  drive a reviewer's forbidden whole-filesystem search caused to be mounted (REPAIR_STATE.md, 2026-10-24 incident).
- 2026-10-24 (batch 18): a size limit on long functions now, their restructuring later (ESCALATED); seqecontain keeps
  accepting a single event sequence (recorded); the 2026-10-08 audit's report and evidence index are tracked.
- 2026-10-25: run the completion of the second prompt.

## From the family's supervision repository

- 2026-11-02 (its `docs/DECISIONS.md`): "Yes upgrade to 4.6.2 and lme4 2.1.0." The oracle's R is 4.6.2.
