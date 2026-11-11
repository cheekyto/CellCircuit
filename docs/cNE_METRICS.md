# cNE metrics

This page lists every measure the three cNE functions compute (`cNEDetection`, `cNEAnalysis` and `cNECompare`). It also marks which measures go into the final compiled table.

- **✓** marks a column of the final compiled table (`Compile_MEAproject`, `<group>_Compiled.xlsx`).
- The other measures are still computed and kept in each recording's full results. They are left out of the compiled table only to keep it readable.
- Exact definitions are in [FUNCTION_REFERENCE.md](FUNCTION_REFERENCE.md) and in the headers of the functions.

**Status (November 2026): proposal.** The pipeline (`Analyze_spikeTrain`) does not run the cNE functions yet, and the compiled table has no cNE sheets yet.

---

## At a glance: the final compiled table

There are two new sheets. Each starts with the same columns as the other sheets: Organoid, File, ttx_qc, Condition and dose.

### Sheet "cNEs": one row per condition, each with its own cNEs (22 columns and a note)

| Group | Columns |
|---|---|
| Ensembles | nCNE, medianMembers, unitsInCNE, activationRate, nAtChance, memberParticipation |
| Network bursts | largeWithCNE, shareCore, shareQuiet, nLeader |
| Order in the bursts | orderConsistency, orderMatch, orderLag, orderJitter |
| Small events | smallWithCNE, smallRecruited, smallAgreement, smallConsistency |
| Secondary events | secondaryWithCNE, secondaryRecruited, secondaryAgreement, secondaryConsistency |
| Note | why cells are blank (from noBursting and notDescribed) |

### Sheet "cNEs vs baseline": one row per drug condition (12 columns)

This sheet also has one no-drug row per recording, comparing the 2nd half of the baseline with the 1st half. It shows the change that time alone brings, as a yardstick for the drugs.

| Group | Columns |
|---|---|
| Ensembles detected | nCNE, nCNETest, nFound, nLost, nNew, foundShare |
| Baseline ensembles followed under the drug | rateRatio, memberRateRatio |
| Roles | rolesKept, rolesLost, rolesGained |
| Order | referenceAgreement |

---

## 1. Per condition, each with its own cNEs (`cNEAnalysis`, `report.summary`)

When a measure is a median, it is taken over the ensembles of the condition.

| Metric | What it tells you | Unit | Final |
|---|---|---|---|
| nCNE | number of ensembles | # | ✓ |
| nUnits | units analyzed for ensembles | # | |
| medianMembers | typical ensemble size: median number of member units | # units | ✓ |
| unitsInCNE | share of the analyzed units that belong to at least one ensemble | % | ✓ |
| activationRate | how often the ensembles fire together: median of their activations per minute | per min | ✓ |
| excessRate | activation rate minus the rate expected by chance (median) | per min | |
| nAtChance | ensembles whose activations are mostly what chance gives (quality flag: chance gives at least 50% of their activations) | # | ✓ |
| memberParticipation | share of the members that fire when the ensemble is active (median). It falls with ensemble size, so it is best for following one ensemble across conditions | % | ✓ |
| nLarge, nSmall, nSecondary, nSmallLeftOut | number of events of each type (small events that overlap a burst are left out) | # | |
| largeWithCNE | share of the network bursts with at least one active ensemble | % | ✓ |
| timeOnset … timeQuiet | share of the recording in each of the 6 periods: onset (0.5 s before a burst), core, secondary, tail, small, quiet | % | |
| shareCore | share of the activations inside the burst core (median) | % | ✓ |
| shareQuiet | share of the activations in quiet time, outside all events (median) | % | ✓ |
| shareOnset, shareSecondary, shareTail, shareSmall | share of the activations in the other periods (median) | % | |
| nLeader | ensembles that fire just before the bursts start (in the onset of at least 25% of the bursts) | # | ✓ |
| nCore, nTail, nReverberation, nSmallEvent, nBaseline | ensembles with each of the other roles | # | |
| leaderLead | how early the leaders fire before the burst starts | s | |
| orderConsistency | is the ensemble order the same in every burst? (0 = chance, 1 = the same order every time) | 0–1 | ✓ |
| orderP | p value of the order consistency | p | |
| orderMatch | match of each burst with the usual order (Kendall tau, median over the bursts) | −1 to 1 | ✓ |
| orderSpan, orderSpanFraction | time from the first to the last ensemble of the usual order | s, fraction of the burst | |
| orderLag | typical time between two ensembles in a burst (later minus earlier ensemble of the usual order; median over the pairs) | s | ✓ |
| orderJitter | how much that time varies from burst to burst (median absolute deviation; median over the pairs) | s | ✓ |
| smallWithCNE | share of the small events with at least one active ensemble | % | ✓ |
| smallRecruited | when a small event recruits ensembles, the share of the ensembles it recruits (mean over the small events with an active ensemble; added November 2026) | % | ✓ |
| smallPosition | which ensembles the small events recruit (0 = the leading ones, 1 = the late ones) | 0–1 | |
| smallAgreement | do the small events follow the burst order? (0.5 = chance, 1 = the same order) | 0–1 | ✓ |
| smallConsistency | is the order of the small events the same in every one? (their own order; 0 = chance, 1 = always; added November 2026) | 0–1 | ✓ |
| smallAgreementP, smallConsistencyP | p values | p | |
| secondaryWithCNE | share of the secondary events (reverberations) with at least one active ensemble | % | ✓ |
| secondaryRecruited | when a secondary event recruits ensembles, the share of the ensembles it recruits (mean over the secondary events with an active ensemble; added November 2026) | % | ✓ |
| secondaryPosition | which ensembles the secondary events recruit (0 = the leading ones, 1 = the late ones) | 0–1 | |
| secondaryAgreement | do the secondary events follow the burst order? (0.5 = chance, 1 = the same order) | 0–1 | ✓ |
| secondaryConsistency | is the order of the secondary events the same in every one? (their own order; 0 = chance, 1 = always) | 0–1 | ✓ |
| secondaryAgreementP, secondaryConsistencyP | p values | p | |
| referenceAgreement, referenceAgreementP, referenceMatch | bursts against the order of another condition. Only with a reference order; the comparison sheet reports referenceAgreement | 0–1, p, tau | |
| noBursting, notDescribed | why cells are blank (see below) | yes/no | ✓ as the note |

---

## 2. Drug vs baseline (`cNECompare`, `comparison.summary`)

"Baseline ensembles followed" means the ensembles detected in the baseline, followed under the drug with the same members. Ratios are drug / baseline; medians are over the baseline ensembles.

| Metric | What it tells you | Unit | Final |
|---|---|---|---|
| nCNE | ensembles of the baseline | # | ✓ |
| nCNETest | ensembles detected under the drug on their own | # | ✓ |
| nFound | baseline ensembles found again under the drug (weight similarity at least 0.5) | # | ✓ |
| nLost | baseline ensembles not found again | # | ✓ |
| nNew | drug ensembles that match no baseline ensemble | # | ✓ |
| foundShare | share of the baseline ensembles found again | % | ✓ |
| medianSimilarity | how similar the closest drug ensemble is (median) | cosine | |
| activeShare | share of the baseline ensembles that fire at all under the drug | % | |
| rateRatio | activation rate, drug / baseline (median) | ratio | ✓ |
| excessRatio, nAtChance | rateRatio after removing what chance gives; baseline ensembles at chance (left out of excessRatio) | ratio, # | |
| memberRateRatio | firing rate of the members, drug / baseline (median). If the activations fall more than this, the members fire less *together*, not only less | ratio | ✓ |
| participationChange | change of the member participation, drug − baseline (median) | points | |
| largeCoverChange | change of the share of the bursts with the ensemble active, drug − baseline (median) | points | |
| rolesKept, rolesLost, rolesGained | number of roles (ensemble and role) kept, lost and gained under the drug | # | ✓ |
| referenceAgreement | do the drug bursts follow the baseline order? (0.5 = chance, 1 = the same order) | 0–1 | ✓ |
| referenceAgreementP, referenceMatch | p value of the agreement; match of the drug order with the baseline order (Kendall tau; needs 3 ensembles) | p, tau | |
| nLargeReference, nLargeTest, and the flags of both conditions | bursts in each condition; noBursting and notDescribed of each | #, yes/no | |

---

## 3. Kept only in the full results (never in the compiled table)

- **Per ensemble** (`cNEAnalysis`, `report.cNEs`): members and their cluster ids, stability, rate, chance rate, chance share and flag, excess rate, participation inside and outside the bursts, share and rate in each period, share of each event type where the ensemble is active (largeCover, onsetCover, tailCover, smallCover, secondaryCover), onset lead, time in the usual order, rank, position, roles.
- **Per ensemble, drug vs baseline** (`cNECompare`, `comparison.cNEs`): each rate, member rate, participation and coverage in both conditions and their ratios, silent members, roles in both conditions, the matched drug ensemble and its similarity.
- **Per event** (`cNEAnalysis`, `report.events`): time of each ensemble in each burst, small event and secondary event; active bins; match of each burst with the usual order; share of the bursts where each ensemble is active at each time from the burst start.
- **Detection** (`cNEDetection`): weights, member thresholds, eigenvalues and their bounds, activity traces, active bins, chance events, ICA convergence, resamples, stability.

---

## Why a cell can be blank

- **noBursting**: no network burst in the window. Only the ensembles, their rates and their member participation are reported.
- **notDescribed**: fewer than 50% of the bursts have an active ensemble. The roles, the order and the small- and secondary-event comparisons are left blank.
- **Order measures** (orderConsistency, orderMatch, orderLag, orderJitter) need at least 2 ensembles and 5 bursts.
  - With 2 ensembles, the order is which ensemble leads.
  - orderMatch needs 3 active ensembles in a burst, so it stays blank with 2 ensembles.
  - orderLag and orderJitter need a pair of ensembles active together in at least 3 bursts.
- **Small and secondary events**: the agreement needs the order measures; the agreement and the consistency need at least 3 events with 2 active ensembles.
  - Small events rarely recruit 2 ensembles. In the November 2026 test recordings, smallConsistency could be computed in only 6 of 120 conditions.
  - In the same recordings, smallWithCNE was computed in 36 of 120 conditions and smallRecruited in 23 (in the other 13, no small event recruited an ensemble). The secondary-event measures were computed in only 3.
  - smallRecruited and secondaryRecruited count only the events that recruit an ensemble. They are blank when no event does (smallWithCNE or secondaryWithCNE is then 0%).
- **referenceAgreement**: needs the order measures in the baseline, and at least 3 drug bursts with 2 active ensembles.
- **Fewer than 5 good units** (MIN_UNITS of the compilation): every cNE column is blank.
