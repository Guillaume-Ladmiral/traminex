# API Parity Table (Function-by-Function)

Summary: TraMineR_exports=125, TraMineR_s3=65, TraMineR_total=190, TraMineRextras_exports=42, TraMineRextras_s3=8, TraMineRextras_total=50, R_total_unique=240, Python_public_api=226, exact_matches=129, sanitized_name_matches=24, with_suggestions=67, missing=20

`R_total_unique=240` counts the exported names of both packages plus every S3 method their `NAMESPACE` files register, de-duplicated. `docs/API_PARITY_LEDGER.md` reports a larger figure because it additionally carries one row per verified unexported R internal that traminex promotes to a public name (`scripts/_r_namespace.py::PROMOTED_INTERNALS`); that difference is the whole of the gap between the two numbers.

## Exact name matches
| R | Python |
| --- | --- |
| alphabet | alphabet |
| cpal | cpal |
| createdatadiscrete | createdatadiscrete |
| dissassoc | dissassoc |
| disscenter | disscenter |
| dissCompare | dissCompare |
| dissdomassoc | dissdomassoc |
| dissindic | dissindic |
| dissmergegroups | dissmergegroups |
| dissmfacw | dissmfacw |
| dissrep | dissrep |
| dissrf | dissrf |
| disstree | disstree |
| disstree2dot | disstree2dot |
| disstree2dotp | disstree2dotp |
| disstreedisplay | disstreedisplay |
| disstreeleaf | disstreeleaf |
| dissvar | dissvar |
| FCE_to_TSE | FCE_to_TSE |
| gower_matrix | gower_matrix |
| HSPELL_to_STS | HSPELL_to_STS |
| rowmode | rowmode |
| seqalign | seqalign |
| seqauto | seqauto |
| seqBIC | seqBIC |
| seqcomp | seqcomp |
| seqCompare | seqCompare |
| seqconc | seqconc |
| seqcost | seqcost |
| seqcta | seqcta |
| seqdecomp | seqdecomp |
| seqdef | seqdef |
| seqdHplot | seqdHplot |
| seqdiff | seqdiff |
| seqdim | seqdim |
| seqdist | seqdist |
| seqdistmc | seqdistmc |
| seqdomassoc | seqdomassoc |
| seqdss | seqdss |
| seqdur | seqdur |
| seqe2stm | seqe2stm |
| seqeapplysub | seqeapplysub |
| seqecmpgroup | seqecmpgroup |
| seqeconstraint | seqeconstraint |
| seqecontain | seqecontain |
| seqecreate | seqecreate |
| seqedist | seqedist |
| seqedplot | seqedplot |
| seqefsub | seqefsub |
| seqeid | seqeid |
| seqelength | seqelength |
| seqemlt | seqemlt |
| seqentrans | seqentrans |
| seqerulesdisc | seqerulesdisc |
| seqetm | seqetm |
| seqeweight | seqeweight |
| seqfind | seqfind |
| seqformat | seqformat |
| seqfplot | seqfplot |
| seqfpos | seqfpos |
| seqfposend | seqfposend |
| seqgen | seqgen |
| seqgranularity | seqgranularity |
| seqhasmiss | seqhasmiss |
| seqHtplot | seqHtplot |
| seqibad | seqibad |
| seqici | seqici |
| seqidegrad | seqidegrad |
| seqient | seqient |
| seqimplic | seqimplic |
| seqindic | seqindic |
| seqinsecurity | seqinsecurity |
| seqintegr | seqintegr |
| seqipos | seqipos |
| seqistatd | seqistatd |
| seqivardur | seqivardur |
| seqivolatility | seqivolatility |
| seqlength | seqlength |
| seqLLCP | seqLLCP |
| seqLLCS | seqLLCS |
| seqlogp | seqlogp |
| seqLRT | seqLRT |
| seqmaintokens | seqmaintokens |
| seqMD | seqMD |
| seqmeant | seqmeant |
| seqmodst | seqmodst |
| seqmpos | seqmpos |
| seqmsplot | seqmsplot |
| seqmtplot | seqmtplot |
| seqnum | seqnum |
| seqpcfilter | seqpcfilter |
| seqpcplot | seqpcplot |
| seqplot | seqplot |
| seqplotMD | seqplotMD |
| seqpm | seqpm |
| seqpolyads | seqpolyads |
| seqprecarity | seqprecarity |
| seqprecstart | seqprecstart |
| seqrecode | seqrecode |
| seqrep | seqrep |
| seqrf | seqrf |
| seqrfplot | seqrfplot |
| seqrplot | seqrplot |
| seqsamm | seqsamm |
| seqsammeha | seqsammeha |
| seqsammseq | seqsammseq |
| seqsep | seqsep |
| seqsha | seqsha |
| seqsplot | seqsplot |
| seqST | seqST |
| seqstart | seqstart |
| seqstatd | seqstatd |
| seqstatf | seqstatf |
| seqstatl | seqstatl |
| seqsubm | seqsubm |
| seqsubsn | seqsubsn |
| seqsurv | seqsurv |
| seqtab | seqtab |
| seqtabstocc | seqtabstocc |
| seqtransn | seqtransn |
| seqtrate | seqtrate |
| seqtree | seqtree |
| seqtree2dot | seqtree2dot |
| seqtreedisplay | seqtreedisplay |
| sorti | sorti |
| sortv | sortv |
| stlab | stlab |
| toPersonPeriod | toPersonPeriod |
| TSE_to_STS | TSE_to_STS |

## Sanitized-name matches (R name normalized)
| R | Python |
| --- | --- |
| alphabet<- | alphabet |
| cpal<- | cpal |
| disstree.assign | disstree_assign |
| disstree.get.rules | disstree_get_rules |
| dissvar.grp | dissvar_grp |
| group.p | group_p |
| is.eseq | is_eseq |
| is.seqelist | is_seqelist |
| is.stslist | is_stslist |
| plot.eseq | plot_eseq |
| plot.SAMM | plot_SAMM |
| plot.seqelist | plot_seqelist |
| plot.seqimplic | plot_seqimplic |
| plot.subseqelist | plot_subseqelist |
| plot.subseqelistchisq | plot_subseqelistchisq |
| seqelength<- | seqelength |
| seqeweight<- | seqeweight |
| seqgen.missing | seqgen_missing |
| seqindic.dyn | seqindic_dyn |
| seqplot.rf | seqplot_rf |
| seqplot.tentrop | seqplot_tentrop |
| seqplot.tentrop.m | seqplot_tentrop_m |
| seqrep.grp | seqrep_grp |
| stlab<- | stlab |

## Suggested mappings (top suggestions by similarity)
| R | Suggestions |
| --- | --- |
| [.seqelist | is_seqelist (0.86), plot_seqelist (0.78), seqedist (0.78) |
| [.stslist | STSList (0.88), is_stslist (0.84), is_subseqelist (0.61) |
| [.stslist.freq | STSList (0.67), is_stslist (0.67) |
| [.stslist.modst | is_stslist (0.64), STSList (0.64) |
| [.stslist.statd | is_stslist (0.64), STSList (0.64) |
| [.subseqelist | is_subseqelist (0.89), plot_subseqelist (0.83), is_seqelist (0.75) |
| as.character.seqelist | is_seqelist (0.62) |
| hist.dissassoc | dissassoc (0.78), dissdomassoc (0.69), dissassoc_test (0.64) |
| levels.eseq | is_eseq (0.67), plot_eseq (0.60) |
| levels.seqelist | is_seqelist (0.77), plot_seqelist (0.71), is_subseqelist (0.69) |
| plot.emlt | plot_seqelist (0.73), plot_seqimplic (0.70), plot_SAMM (0.67) |
| plot.seqalign | plot_seqelist (0.77), SeqAlign (0.76), seqalign (0.76) |
| plot.seqdiff | plot_eseq (0.76), seqdiff (0.74), plot_seqelist (0.72) |
| plot.seqpcplot | seqpcplot (0.78), plot_seqelist (0.74), plot_seqimplic (0.71) |
| plot.seqrf | seqrfplot (1.00) |
| plot.stslist | seqplot (1.00), sequence_index_plot (1.00) |
| plot.stslist.freq | seqfplot (1.00) |
| plot.stslist.meant | seqmtplot (1.00) |
| plot.stslist.modst | seqmsplot (1.00) |
| plot.stslist.rep | seqrplot (1.00) |
| plot.stslist.statd | seqHtplot (1.00), seqdHplot (1.00) |
| plot.stslist.surv | plot_seqelist (0.67), plot_subseqelist (0.67), plot_subseqelistchisq (0.63) |
| print.diss.rep | dissrep (0.67) |
| print.dissassoc | dissassoc (0.75), dissdomassoc (0.67), dissassoc_test (0.62) |
| print.dissrf | dissrf (0.67) |
| print.disstree | disstree (0.73), disstree2dot (0.62), disstree_fit (0.62) |
| print.eseq | plot_eseq (0.74), is_eseq (0.71) |
| print.eseqsum | plot_eseq (0.64), is_eseq (0.60), seqsubm (0.60) |
| print.sdomassoc | dissdomassoc (0.74), seqdomassoc (0.69) |
| print.seqalign | SeqAlign (0.73), seqalign (0.73), seqalignplot (0.62) |
| print.seqdiff | seqdiff (0.70) |
| print.seqeconstraint | seqeconstraint (0.82), seqecontain (0.71), seqentrans (0.60) |
| print.seqelist | plot_seqelist (0.81), is_seqelist (0.80), plot_subseqelist (0.73) |
| print.seqidegrad | seqidegrad (0.77) |
| print.seqimplic | plot_seqimplic (0.83), seqimplic (0.75) |
| print.seqipos | seqipos (0.70), plot_seqelist (0.62), seqfpos (0.60) |
| print.seqivardur | seqivardur (0.77) |
| print.seqrf | seqrf (0.62), plot_eseq (0.60) |
| print.seqtree | seqtree (0.70) |
| print.stslist | is_stslist (0.78), STSList (0.70), plot_subseqelist (0.62) |
| print.stslist.freq | is_stslist (0.64) |
| print.stslist.meant | is_stslist (0.62) |
| print.stslist.modst | is_stslist (0.62) |
| print.stslist.rep | is_stslist (0.67) |
| print.stslist.statd | is_stslist (0.62) |
| print.subseqelist | plot_subseqelist (0.85), is_subseqelist (0.84), plot_subseqelistchisq (0.74) |
| rbind.stslist | is_stslist (0.78), STSList (0.70) |
| read.tda.mdist | createdatadiscrete (0.62) |
| recodef | seqrecode (0.75), seqdef (0.62) |
| seqdplot | seqdHplot (0.94), seqedplot (0.94), seqplot (0.93) |
| seqfcheck | seqconc (0.62) |
| seqintegration | seqintegr (0.78), seqentrans (0.67), seqidegrad (0.67) |
| seqIplot | seqplot (0.93), seqfplot (0.88), seqrplot (0.88) |
| seqiplot | seqplot (0.93), seqfplot (0.88), seqrplot (0.88) |
| seqlegend | seqgen (0.80), state_legend (0.76), SeqAlign (0.71) |
| seqlength.align | seqlength (0.75), seqelength (0.72), SeqAlign (0.70) |
| seqprecorr | seqrecode (0.74), seqprecarity (0.73), seqprecstart (0.73) |
| str.eseq | is_eseq (0.80), plot_eseq (0.71), disstree (0.62) |
| str.seqelist | is_seqelist (0.87), plot_seqelist (0.80), is_subseqelist (0.77) |
| summary.ddomassoc | seqdomassoc (0.64), dissdomassoc (0.62) |
| summary.diss.rep | dissrep (0.61) |
| summary.dissrf | dissrf (0.60) |
| summary.eseq | is_eseq (0.63) |
| summary.seqelist | is_seqelist (0.74), is_subseqelist (0.67), plot_subseqelist (0.62) |
| summary.seqpcplot | seqpcplot (0.69) |
| summary.stslist | is_stslist (0.72), STSList (0.64) |
| summary.stslist.rep | is_stslist (0.62) |

## Missing (no suggestion >= 0.60)
| R |  |
| --- | --- |
| as.character.eseq |  |
| convert.g |  |
| hist.TraMineRPermut |  |
| pamward |  |
| plot.dynin |  |
| print.dissmultifactor |  |
| print.dynin |  |
| print.emlt |  |
| print.seqprec |  |
| print.TraMineRPermut |  |
| summary.seqrf |  |
| summary.TraMineRPermut |  |
| TraMineR.check.depr.args |  |
| TraMineRInternalLayout |  |
| TraMineRInternalLegend |  |
| TraMineRInternalNodeInit |  |
| TraMineRInternalSeqeage |  |
| TraMineRInternalSeqgbar |  |
| TraMineRInternalSplitInit |  |
| TraMineRInternalWeightedInertiaDist |  |
