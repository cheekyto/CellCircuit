# CellCircuit: function reference

Detailed documentation of the stages and functions of the pipeline for Multielectrode array (MEA) recordings (see the [README](../README.md) for an overview, and `run_pipeline` to run it).



The pipeline has a pre-specified order:

0\. Kilosort + BombCell (spike sorting)

1\. ShapeSpikeTrainData

2\. QualityControl\_General

3\. QualityControl\_HyperSynch

4\. QualityControl\_Manual

5\. SelectElectrodes (optional)

6\. IntegrateDataSets (Kilosort and/or 3Brain datasets; it generates the file "XXX\_curated\_QCmanual.mat", or "XXX\_integrated\_QCmanual.mat" when both datasets are combined)

7\. Analyze\_spikeTrain

8\. Compile\_MEAproject



## Revision of October 2026

The whole pipeline was revised, phase by phase (the development history is kept in the lab's repository; this public repository starts from the revised version). The main points are:

- **Inputs:** every stage function accepts an optional structure `inputs` that overwrites any variable of its INPUTS / PARAMETERS sections (see `override_inputs.m`), e.g. `Analyze_spikeTrain(struct('folders', "...", 'files_data', "...", 'condition_analysis', -1))`. Without it, the values written in the code are used, as before. `inputs.system = 0` skips the `addpath` at the beginning.
- **Tests:** the folder `tests` contains a regression harness that runs each stage with the old and the new code on a reference recording and compares every saved variable, plus function-level tests (see `tests/README.md`).
- **Speed** (reference recording, ~800 clusters before QC): `QualityControl_HyperSynch` ~850 h -> ~30 min; `Analyze_spikeTrain` ~50 min -> ~2 min per recording condition. Most speed-ups give exactly the same results.
- **Changes in the results** (compare with results obtained before the revision with care):
    - Fano factor = CV^2 of the inter-spike / inter-burst intervals (`BurstFeatures`, `SpikeFeatures`).
    - `QualityControl_HyperSynch`: shuffled spike trains generated with an equivalent, faster random process; compact output file (see below).
    - Network bursts (`NetworkBurstDetection` and the functions it calls):
        - the participating units are counted only in solid events, for both the tonic and the phasic rules;
        - a level-0 primary event that coincides with a secondary event is kept only as a secondary event;
        - the edges of each event are set at the interpolated crossing of 20% of the peak height;
        - the activity filter (0.03 Hz) is applied to every step of the detection;
        - the tail of the population peaks is a ratio of activity per unit of time;
        - the secondary-per-primary metrics are normalized to the level-1 primary events;
        - the temporal structure (`NetworkBurstStructure`) uses the whole spike train of each neuron, between the neighbouring events;
        - each family of primary events has its own participating units (`familyParticipants`);
        - the clustering is reproducible (fixed random seed).
    - `Compile_MEAproject`: the recording conditions that were not analyzed yet are skipped (with a warning).
- **New helper functions:** `override_inputs`, `count_below`, `windowCounts`, `shuffleCache`, `shuffleHits`, `mergeSpikeTrains`, `familyParticipants`.

## Upgrade of October 2026

After the revision:

- **Variability of the spiking and bursting activity** (`SpikeFeatures`, `BurstFeatures`, `NetworkBurstFeatures`, new helper `countFanoFactor`). The previous `fanoFactor` (= CV^2) was replaced by:
    - the CV of the intervals (ISI, IBI);
    - the interval Fano factor in two forms, for the ISI, the IBI, the ISI within bursts and the ISI outside bursts: `FF_cv2` = CV^2 (dimensionless, the same value as the previous `fanoFactor`) and `FF_varmean` = variance / mean of the intervals (in ms for the ISI, in sec for the IBI);
    - the count Fano factor `FFcount` = variance / mean of the number of events (spikes, bursts, network bursts) counted in consecutive windows of 60 sec (`COUNT_WINDOW`), with the mean count per window `meanCount` (FFcount is biased toward 1 when there are less than ~5-10 events per window). For a stationary renewal process FFcount ~ CV^2 for long windows; FFcount > CV^2 suggests slow changes of the rate (non-stationarity). Reference: Rajdl, Lansky & Kostal 2020, Front. Comput. Neurosci. 14:569049;
    - Excel columns (`Compile_MEAproject`): Firing `CV_ISI, FFisi_CV2, FFisi_VarMean, FFcount_Spikes, MeanCount_Spikes` (they replace `Fano_Factor`); Bursting `CV_IBI, FFibi_CV2, FFibi_VarMean, FFcount_Burst, MeanCount_Burst` (they replace `RateFF_Burst`), `CV_ISI_Burst, FFisi_CV2_Burst, FFisi_VarMean_Burst` (they replace `AmplitudeFF_Burst`), and the same for `_nonBurst` and `_Tonic` (they replace `RateFF_Spiking_nonBurst` and `RateFF_Spiking_Tonic`); network bursts `FFibi_CV2_NETburst, FFibi_VarMean_NETburst, FFcount_NETburst, MeanCount_NETburst` (the CV is `RateCV_NETburst`).
- **Clustering of the network bursts:** `Analyze_spikeTrain` prints `quality_clustering`, and the network-burst tables have the column `Clustering_Quality` (2 = reliable; 1 = it seems to work; 0 = not significant, single family; NaN = not run). Since November 2026 it is the stability of the least stable family (0 to 1; 0 = no family; NaN = not run): see the families of network burst events below.
- **Secondary events:** the minimum number of participating units is 6 during the detection (`Optimal_secondaryPeak`, it was 5), as in the final quality control. In the last round of the detection (`NetworkBurstDetection`), after the quality control of each secondary event, the secondary events are accepted only if at least `GAP_WITH_PEAK` = 40% (75% until November 2026, then 50%; more than it, until November 2026) of the gaps after the level-1 primary events have at least `MIN_PEAK2_GAP` = 1 secondary event; otherwise they are all rejected, and the level-0 primary events that coincided with them are kept as primary events (the rule was in `NetworkBurstPlayers`, where it never ran). On the reference recording 98% of the gaps have secondary events: the results do not change.
- **Master function and single parameter file** (`run_pipeline`, `pipeline_parameters`, `load_parameters`, `apply_parameters`, `parameter_set`): the recordings are listed in a table (`recordings_template.csv`; `make_recordings_table` builds the columns `folder` and `file_name` from a folder, one row per recording found as `<name>_QCmanual.mat` in it or its subfolders), and `run_pipeline` runs the pipeline end to end (Kilosort) or stage by stage (also the 3Brain stages). Every parameter of every function (267) is set in `pipeline_parameters.m`; while `run_pipeline` runs, each function replaces the values written in its code by those of the file (a stage run by hand still uses the values in its code). The set used is saved for each recording (`<file_name>_parameters.mat`, with the history of the runs and the version of the code), and a run can reuse it (`'parameters', "saved"`), use a different file, or change a few values (`'changes'`). Since November 2026, a run with the saved sets can take some of their parts from `pipeline_parameters.m` instead (`'renew'`: a whole function, e.g., one revised after the last run, whose old parameters would stop the run, or single values); the rest of each saved set, such as the settings of a recording, is kept. `QualityControl_Manual` runs without windows (the clusters in the column `to_discard` are discarded). Only the progress is printed; the messages of the stages go to `<file_name>_pipeline_log.txt`. With the parameter set active, every stage gives exactly the same results as before on the reference recording.
- **Network bursts, option 4 (November 2026):** `NetworkBurstDetection(..., 4, ...)` (`NETWORK_BURST_DETECT_METHOD = 4` in `Analyze_spikeTrain`) runs the procedure of option 3 (screen of bin sizes, participating units, secondary events, families), but the population signal is the **participation**, the fraction of units with at least one spike in each bin, and the thresholds come from **surrogates** instead of standard deviations (`PopulationActivityDetection` with `opt_norm = 2`; new helper `participationSurrogate`). Each surrogate shifts the spike train of every unit circularly by its own random offset (N_SURROGATE = 1000; it keeps each unit's spikes and intervals and breaks the timing between units). Large (level-1) events: peak above the 95th percentile of the maxima of the surrogates (max statistic, ALPHA_LARGE = 0.05: a 5% chance of a false large event anywhere in the recording); each event gets its p-value, (1 + number of surrogate maxima >= peak)/(N_SURROGATE + 1) (`IndivBurst.pSurrogate`). Small (level-0) events: peak above the 99th percentile of all the surrogate bins pooled (PERCENTILE_SMALL; Bocchio et al. 2020, Nat Commun 11:4559). Secondary events: peak above the 95th percentile of the maxima of surrogates made within each gap (ALPHA_SECONDARY). No prominence is used: each stretch of the smoothed signal above the threshold is a candidate, with its peak at its highest point; its edges are where the participation falls to 50% of the height of the peak above the median of the signal (event window; half-height above the typical level of the recording) and to 20% (reported). Candidates whose 50% edges overlap are merged into one event spanning all of them. Tail fragments are dropped: a stretch with no rise of its own (its highest point is its first bin) that starts right where the signal had just left the threshold (a dip of at most 2 bins, not higher than where it left) or at the first bin of the window. The participation of the units is tested on the core of each event (start to peak, plus the same time after the peak), so that the persistent firing after large events, which is part of their duration, does not dilute it. The size check and everything downstream as in option 3. Its own bin settings: BIN_MIN = 0.02 s, BIN_MAX = 1.0 s, BIN_DEFAULT = 0.10 s.
- **Network bursts, option 5 (November 2026):** `NetworkBurstDetection(..., 5, ...)` runs option 4 first and keeps only its large (level-1) events, then runs option 3 as usual. A large event of option 4 that matches no event of option 3 (large, small or secondary; same event when the peaks are within `ADD_MATCH_TOL` = 0.5 sec or one peak is inside the other event) is added to option 3 as a small (level-0) primary event, after all the checks of option 3. It is characterized by the same code as the events of option 3, on the signal of the final detection of option 3 (same bin, units and smoothing): the peak is the highest local maximum of that signal within the event of option 4, the edges come from `QuickNetBurstLimits` and `characterizeEvent` (20% and 50% levels), and the size check (`PopulationActivitySize`), the tail (`tailPopulationPeaks`) and the spikes of each unit (`eventSpikes`) are those of option 3. Option 5 has no copy of options 3 and 4: it calls option 4 as `NetworkBurstDetection(..., 4, ...)`, runs the code of option 3 with its parameters, and the steps of option 3 that it reuses were moved to shared functions (`characterizeEvent`, from `PopulationActivityDetection`; `eventSpikes`, from `NetworkBurstPlayers`; same results), so a change to option 3 or 4 also applies to option 5. Not added: events of option 4 when option 3 finds no network bursting, events whose peak is inside a large event of option 3 or its tail, and events whose edges overlap an event of option 3. The added events are not tested by the recruitment check (they passed the surrogate test of option 4 as large events), but their relative recruitment is printed. Each event has `origin` = 3 (option 3) or 4 (added), and `NetworkBurstDetection.option5` counts what happened to the large events of option 4. The small and secondary events of option 4 are never used.
- **Network bursts, options 1-3 (November 2026):** candidate events whose edges (20% of the prominence) overlap are merged into one event spanning all of them, with the peak of the highest one; before, one of them was dropped. All the options report two durations of each event: at 20% and at 50% (half-width) of its height (`IndivBurst.start20/last20/start50/last50`; `duration20`, `duration50` in `NetworkBurstFeatures`, and their columns in `Compile_NetBurstResults`). Options 1-3 keep the 20% window as the event window; option 4 keeps the 50% window. Small (level-0) events of options 1-3: the tail of each large event lasts until the signal returns to its median (baseline), and at most `TAIL_FACTOR` = 3 times the duration of the event after its end, whichever comes first. The tails are left out of the noise estimate of the second round of the search (before, each large event plus the same duration again), and no small event is accepted with its peak inside a large event (from its start) or in its tail (before, the small events were searched in the whole signal, so the bumps on the tails, where the signal is higher and noisier than the baseline that sets the threshold, became small events).
- **Secondary events, separation at the bin of the secondary search (November 2026):** `Optimal_secondaryPeak` builds the population signal of each gap at the bin of the secondary search (the same type of signal and smoothing as the search). The search starts where the large event that opens the gap falls to 50% of its height above the median of the gap (START_LEVEL_FINE), instead of a point set by the bin of the large events (after their end, at their end, or halfway down their decay, which missed the first reverberations). Each secondary event must be separate from its large event: between the two peaks, the signal must fall below 50% of the height of the secondary peak above the median (SEPARATION_LEVEL; the 50% rule used to merge the events of option 4) and, for option 4, by at least 3 times its binomial sampling noise (MIN_DROP_Z; with few units the participation is noisy). The drop is reported (`Secondary.IndivBurst.dropOwn`, `dropZ`, and the start of the search, `searchStart`). In `NetworkBurstDetection`, these events are no longer compared with the end of their large event, measured at the bin of the large events, where the large event and its reverberations blur into one mountain.
- **Secondary events, search of the reverberations (November 2026):** on recordings with clear reverberations after the large events (2-3 bumps each), option 3 found almost none. (1) The minimum prominence of the search is 3 x the standard deviation of the search window, and the window starts on the decay of the large event, which raised that standard deviation 1.4-4.4 times: the decay, until the signal first falls to the median of the window, is now left out of the noise (`PopulationActivityDetection`, `SECONDARY_NOISE_TAIL` = 1). (2) Consecutive secondary events of a gap that the separation rule does not separate (the signal between them falls below `SEPARATION_LEVEL` of the lower peak, by `MIN_DROP_Z` times its noise) are one event, with the higher peak and edges spanning both (`Optimal_secondaryPeak`, `MERGE_SECONDARY` = 1); before, one reverberation split into several events at the fine bins. (3) Options 1-3: each secondary event must rise above the activity around it (median of the signal 0.5-2 s before and after its peak, outside the large event) by at least `MIN_RISE_Z` = 3 times the Poisson noise of the two; on a long elevated tail, the fluctuations recruit more units than chance but are not distinct events. (4) The finest bin of the search is 15 ms (`BIN_MIN`, 10 ms before). (5) The recruitment check of the secondary events (`NetworkBurstDetection`) measures their core, between their 50% edges (`RECRUIT_SECONDARY_CORE` = 1; before, from start to last, where on wide events the chance level of the decay rose with the window). (6) `GAP_WITH_PEAK` is a minimum: exactly 40% of the gaps is enough. On the 18 test recordings and two more, the secondary events changed only where the user judged them: a recording with clear reverberations went from 0 to 7 (70% of its large events followed), one whose 5 secondary events were on the falling edge of their large events (not real) went to 0, one with a regular reverberation rhythm went from 40 to 44 (all 40 kept), and the others were unchanged; the large events did not change (peaks within 0.1 s). Items 1, 2, 3 and 5 can be switched off with their parameter (0), and `BIN_MIN` set back to 0.010.
- **Network-burst events, no duplicate after the refinement (November 2026):** `refinePeakDetection` searches the peak of each event in its window extended by half its width on each side (`WINDOW_FACTOR` = 0.5). A broad event next to a higher one (e.g., a small event 2-4 s after a large one, or a large event next to a larger one) could take the peak and the edges of its neighbour: the same event was reported twice and the event itself was lost. `PopulationActivityDetection` now refines again only such an event (its refined peak outside its own window and inside another event), with its peak kept inside its own window (`REFINE_OWN_PEAK` = 1; new optional input `peak_inside` of `characterizeEvent` and `refinePeakDetection`); every other event is computed as before. Keeping every peak inside its own window was tested and not adopted: the refinement runs in every round of the detection, and it changed the events of 6 of 20 test recording-conditions that had no duplicate. On the 18 test recordings and one more, 19 of 21 recording-conditions are identical; the two with duplicates lose them (8 and 1 pairs), with the same large events (peaks, edges and durations). Among 55 recordings analyzed before the change, 4 had duplicates, and only one of them reports network bursting.
- **Bin size of the network bursts, plateau selection (November 2026):** `DynamicPopActDetect` runs the detection at a series of bin sizes (as before) and compares the level-1 events of each pair of neighbouring bin sizes: they are matched one to one (peaks within the larger of the two bin sizes, at least `MATCH_MIN` = 0.1 s), and their agreement is matched / all the different events. A plateau is a run of neighbouring bin sizes, within `BIN_MIN` and `BIN_MAX` of the method, whose events agree at least 90% (`AGREEMENT`; if there is none, 80%, then 75%), with at least `MIN_EVENTS` = 2 level-1 events each. The widest plateau (in log scale; then the most stable, then the one with more events) gives the bin size at its geometric centre; with no plateau, the bin size is not identified (default bin size, as before). In `NetworkBurstDetection`, the final bin size is the centre of the plateau of the refined search (before, a weighted average of the two rounds). Before, the bin size was the one with a small change of the number of level-1 events and a low CV of their intervals: the ties between bin sizes with the same events were broken towards the finest one, at the edge of the plateau, so small changes of the data moved the choice (e.g., from 0.10 to 0.025 s). With 10% of the units silenced (3 random draws), the plateau bin size was identical on the test recordings with a plateau, while the previous rule changed it on 4 of 5; the large events of options 3 and 4 agreed as well or better. The durations of the events depend on the bin size, so they differ from the previous results. `BIN_SELECTION = 'legacy'` restores the previous rule (not the weighted average).
- **Small events of options 1-3, robust noise and recruitment (November 2026):** in the second round of the search (`PopulationActivityDetection`), the noise of the signal without the large events and their tails is 1.4826 x its median absolute deviation, with a floor at `NOISE_FLOOR_SD` = 0.25 x its standard deviation (`NOISE_ROBUST` = 1; before, the standard deviation, which the small events themselves inflated, so on recordings with many moderate events most of them were missed). In `NetworkBurstDetection`, a small event is then kept only if, in a window of `RECRUIT_PEAK_BINS` = 2 detection bins centred on its peak, (1) the fraction f of the participating units with at least one spike is above the `RECRUIT_PERCENTILE` = 99th percentile of surrogates (each unit circularly shifted; built outside the large events and their tails) and (2) its relative recruitment r = (f - b)/(F_large - b) is at least `RECRUIT_REL_MIN` = 0.20, with F_large the median f of the large events of the same recording and b the local baseline of the event (r = 1, as a typical large event; r = 0, as the baseline). With many units, a small but coherent drift of the rate is far above chance; r measures its size against the network bursts of the recording. The baseline b is the same fraction in the window sliding over the `RECRUIT_LOCAL_BINS` = 20 bins before and after the event (outside the large events and their tails), on its quieter side (the lower of the two medians, so that a neighbouring event does not raise it); with fewer than 5 such windows, or with `RECRUIT_LOCAL_BINS` = 0, the chance level f0 (median of the surrogates). With f0 (and a minimum of 0.15, before), the bumps on the slow decay after a large event passed, because the units there still fire above their average. The whole event (start to last) is not used for the small events: their edges are far from the peak, and over long windows most units fire by chance. Secondary events are tested within the event against surrogates built on the decay of the large events.
- **Robust units, by class of events (November 2026):** in the first screening of `NetworkBurstPlayers` (option 1, robust vs weak), the robust units are found separately in the solid large (level-1) events and in the solid small (level-0) events, each counted over its own events, and a unit is robust if it is robust in either: `MIN_PARTICIPATION_TONIC`, `MIN_PARTICIPATION_PHASIC` and `MIN_SOLID_BURST` (the minimum number of events of the count) for the large events, and `MIN_PARTICIPATION_TONIC_SMALL`, `MIN_PARTICIPATION_PHASIC_SMALL` and `MIN_SOLID_BURST_SMALL` for the small events, which are used only when there are at least `MIN_SOLID_BURST_SMALL` solid small events. The thresholds are 0.75 (tonic units: fire in more than 75% of the solid events of the class) and 0.50 (bursting units: burst in more than 50% of them) for both classes, and 3 events (user's choices; before, 0.65 and 0.25). Before, both classes were pooled, so many small events (e.g., bumps on the tails of large events) diluted the participation in the large events and moved the robust units (on a test recording, 118 vs 127 units with and without those bumps). The units that take part only in small events are still found. On the test recordings (with the other November changes), the robust units changed by at most 2, and no unit was robust only through the small events. The activity filter (`findChannel_in_Peak`, `ACTIVITY_THRESHOLD`) and the labelling of the solid events still use all the events together.
- **Families of network burst events, by their participating units (November 2026; user's choices):** `cluster_neuronal_events` groups the events (primary, and secondary when there are more than 5) with the same composite distance as before, 0.65 x the Jaccard distance between their participating units + 0.35 x the distance between their 7 scalar metrics (duration, spikes, active units, ISI, time since the previous event, tail, amplitude), with three changes. (1) A unit participates in an event when it fires significantly more spikes within the event (start to last) than expected from its own rate outside all the events of the condition (Poisson, `PARTICIPATION_ALPHA` = 0.01, and at least `PARTICIPATION_MIN_SPIKES` = 2 spikes); before, when it fired at least one spike within the window of the event (plus the extra time of `eventSpikes`): in long events most units fired by chance, so events driven by different units looked alike (two small events with no unit in common shared 43% of their units), and every large event had almost all the units. (2) The families come from an average-linkage tree of the composite distance, cut where the events of a family are, on average, at least `MIN_SIMILARITY` = 0.50 similar, and a family needs at least `MIN_FAMILY_EVENTS` = 3 events; the events of smaller groups belong to no family (label 0). Before, Ward linkage on an MDS projection, with 2 to `MAX_K` = 4 families chosen by the silhouette score, a bootstrap of the events and a permutation test (`bootstrap_optimal_k`, `permutation_test_clustering`, no longer used: moved to the folder `legacy`): the identical large events gave a high silhouette to any split of the other events, the permutation test shuffled the order of the lists of units, which leaves their distances unchanged, and there was no minimum size, so families of 1-2 events were common. (3) The stability of each family (`clustering.stability`) is the mean proportion of `N_RESAMPLE` = 200 resamples of the participating units (with replacement) in which two of its events are again in the same family; `clusteringQuality` is the stability of the least stable family (0 when there is no family). `clustering.participation` and `clustering.units` keep the participating units of each event. The new helper `primaryFamilies` gives the groups of primary events used by `Analyze_spikeTrain`, `Compile_NetBurstResults` and `Compile_MEAproject`: the families, then one group with all the remaining events, which belong to no family (label 0; user's choice: one row for them in the per-family results). The new column `Family_Type` of the network-burst tables tells the rows apart: "family", "remaining events (no family)" (not a family: label 0, with its own features and participating units, measured as a family), "all events" (the pooled table), or "no network bursting" (a condition without primary events). When there is no family, the only row is the remaining events (all the primary events). `Number_clusters_Netburst` counts the real families only, and the table of the largest family (NETbursting1) chooses among the real families. The minimum rate of events of each row is now `MIN_NETBURST` = 1/210 per sec (1/200 before), so that a family of 3 events in a 600-s window is measured (before, it had only its number of events). In the network-burst figure, the events in no family are grey and the secondary events magenta. On the 18 test recordings and two more (option 5), the families of 1-2 events disappeared (e.g., 4 families of 14, 2, 1 and 1 events -> 1 family of the 14 large events and 4 small events in no family; 4 families of 16, 5, 5 and 1 -> 16, 5 and 3, and 3 events in no family), and the clear families were kept (e.g., 11 large + 43 small events, stability 0.97-1). On the 14 recording-conditions with families, the 8 that had families of 1-2 events have none now. With the scalar metrics in the distance, weaker large events can fall out of the family of the other large events (2 of 11 on one recording, with half the amplitude and 73-77% of the units active, vs 88-100%; 5 of 79 on another). Test: `tests/test_family_clustering.m`.
- **Recordings with few or no units (November 2026; user's choices):** the pipeline runs end to end when no unit passes the quality control (before, it stopped at the plots of `QualityControl_Manual`, then in `Analyze_spikeTrain`). `QualityControl_Manual` skips its plots when no cluster is left and saves the empty result; `IntegrateDataSets` keeps the other dataset when the Kilosort or the 3Brain dataset has no cluster left (the integration of the two was not tested with real 3Brain data). `Analyze_spikeTrain` analyzes a recording only if it has at least `MIN_UNITS` = 5 good units: with fewer (0 to 4), each condition is saved with its label, dose, window and file only (`PerCondition.too_few_units` = true; false for the analyzed conditions), and no raster plot is saved without units. `Compile_MEAproject` reports every recording with fewer than `MIN_UNITS` good units (also one analyzed before this rule) with 0 active neurons (`Number_Good_Clusters`, `Number_Kilosort` and `Number_3Brain` = 0) and 0 Hz (`Firing_Rate`); every other metric is blank (NaN), in the bursting and network-bursting tables too, with one row per condition and `Family_Type` = "fewer than 5 units". `Compile_MEAproject.MIN_UNITS` is linked to `Analyze_spikeTrain.MIN_UNITS`. Tested on copies of a test recording with 0, 1, 2, 3 and 9 units (including cropping that leaves no unit); on recordings with 5 units or more, every saved result is identical to before, except the new field.
- **Window `CellCircuit` (November 2026; user's choices):** `CellCircuit` opens a window to run `run_pipeline` without writing code, for users without coding experience: tab Recordings (open, build from a folder and edit the recordings table; the recordings of a folder are the folders with files of the pipeline, or with Kilosort output), tab Settings (what to run, which parameters, single values changed for this run, options), tab Run (summary and equivalent `run_pipeline` command, progress table of recordings x stages filled in during the run, messages, logs and Excel results), and a bar with Run and Stop. The parameters have two levels: the common adjustments (`common_parameters.m`: about 15 values that were changed for particular recordings, with a short name and what they do) are shown first, and every other parameter of the parameter file is under "All parameters (advanced)", with a search. Each user can move parameters from the full list to the common adjustments and back (buttons "Add to the common adjustments" and "Remove from this list"; at most 15, user's choice); the list is kept in the MATLAB preferences of that user, and "Default list" goes back to `common_parameters.m`. The parameters added by a user are shown with their name and the description of the parameter file. A rate can be shown as its period (column `display` of `common_parameters`, user's choice): `Compile_MEAproject.MIN_NETBURST` (1/210 events per s) is shown and typed in seconds, as the longest average time between the events of a family (210 s), and the window passes 1/seconds to the pipeline (in the full list, the change is written as `1/(seconds)`). The values typed are checked against the type of the value of the file, highlighted, passed as `'changes'`, and can be saved as a new parameter file (a copy of the file with the changed values at its end). The npy-matlab folder and the last table are remembered (MATLAB preferences). `run_pipeline` has a new option `'progress'`: a function called at each step of the run (start, each stage of each recording, each group, end); when it returns true, the run stops after the current stage (the history of the recording says "stopped by the user"; the compilation does not run). Without it, `run_pipeline` runs as before. `tools/package_CellCircuit(version)` builds the install file of the toolbox (`.mltbx`, from a clean copy of the last commit, so uncommitted changes never reach it: the functions, the window in the Apps tab, README and this reference); `CellCircuit_version.txt`, written in the package, is the version of the code saved in the history of each recording when the code is not a git repository. Test: `tests/test_CellCircuit.m` (the window was also run on copies of a test recording: full run, stop after the first recording).
- **Coordinated neuronal ensembles in the pipeline (November 2026; user's choices):** `Analyze_spikeTrain` runs the new step `Analyze_cNE` after the other results of the recording are saved (`CNE_ANALYSIS` = 1, on by default; 0 = not run): the cNEs of each condition on its own analysis window (`cNEDetection` + `cNEAnalysis`), each other condition against the baseline (`cNECompare`, the baseline cNEs detected once), each condition against the condition recorded just before it when that one is not the baseline (e.g., drug 3 against drug 2; new helper `previousCondition`, by the start times of the analysis windows), and the 2nd half of the baseline against its 1st half (the change that time alone brings). The TTX conditions are left out of every cNE step and of the cNE sheets (user's choice). Everything is saved in the new variable `cNEResults` of `XXX_Analyzed.mat` (kept out of `PerCondition`, so that the compilation stays fast); in a rerun, the conditions that were not analyzed again keep their cNEs, and only the comparisons that involve the analyzed ones are computed again. A step that fails is reported (warning and status) and the rest goes on. With `print_Plots`, one figure per condition, per comparison and for the baseline halves (`XXX_cNE_<condition>.fig`, `XXX_cNECompare_<condition>.fig`, `XXX_cNECompare_<condition>_vs_<condition before>.fig`, `XXX_cNECompare_BaselineHalves.fig`; with the index of the condition when two conditions have the same name), hidden in long batches like the other plots. `Compile_MEAproject` writes three new sheets with the short list of `docs/cNE_METRICS.md` (`CNE_COLUMNS`, `CNE_COMPARE_COLUMNS`): "cNEs" (one row per condition, 22 measures and a note when they are blank), "cNEs vs baseline" (one row per condition against the baseline, and one for the baseline halves; 12 measures, the condition compared with, and a note) and "cNEs vs previous" (one row per condition against the condition before it; the same columns); recordings analyzed without cNEs have blank cells and the note "cNEs not computed". Run time: about 1-4 min per recording (up to about 10 min). The settings of `cNEDetection`, `cNEAnalysis` and `cNECompare` are in `pipeline_parameters.m` (section 7f; the same values as in the functions), so that a parameter file, `'changes'` or the window `CellCircuit` (All parameters) can change them; a saved parameter set without them uses the values written in the functions, which are the same. Test: `tests/test_Analyze_cNE.m`.
- **Network-burst plots of conditions with the same label (November 2026):** `Analyze_spikeTrain` saved the plot of each condition as `XXX_NetBurst_<condition><dose>.fig`, so when two conditions had the same label and dose (e.g., two TTX files, both "+TTX"), the second plot overwrote the first one. The labels of every condition are now computed before the analysis (the same code, moved to the local function `label_condition`), and when two conditions share a label the index of the condition is added to the names of their plots (`XXX_NetBurst_+TTX_2.fig`, `XXX_NetBurst_+TTX_3.fig`), as for the cNE figures. The other names, and every result, are unchanged (checked on a test recording). Test: `tests/test_baseline_condition.m`.
- **Cropping from the photo of the chip (November 2026; user's choices):** `OrganoidsFromPhoto` finds, in a photo of a 3Brain chip (64 x 64 electrodes), the electrode array and the organoids on it, and gives for each organoid its outline on the electrode map, the electrodes to keep and the equivalent ellipse (`x0 y0 radiusx radiusy`, as `crop_ellipse`), so that the ellipse does not have to be set by hand. Only the photo is used. The photo is turned so that the "3D" text of the chip (outside the big circle) is on the right; then electrode (1,1) is the bottom-left corner of the array, x grows to the right and y upwards (as `positionXY + 1`). The electrode array is found from the pattern of the edges around it, the same on every chip and centered on the array (November 2026, revised): the dark square of the package with its 4 rows of leads, its rounded opening, and the chip with the strip on the left of the array and a band above and below it; the pattern (`OrganoidsFromPhoto_chip.mat`, next to the function: mean edges of 49 photos straightened on the electrode coordinates) is searched at all sizes (array from 2% to 32% of the photo) around the largest pale region and around the centers of the concentric circles of the well, then fitted in the photo turned to the reference orientation (size within 1%, tilt within 0.5 degree). The side of the text is the turn where the pattern fits best (strip on the left of the array, slot on its right), checked against the gold landmarks of the board around the well (the connector, 2 columns of contacts on the left when the text is on the right; the triangle, top-right corner; the text itself); a landmark that disagrees is reported only when the pattern does not settle the side. `TEXT_SIDE` is optional: when given, it is used, and a mismatch with the pattern is reported (warning "ORIENTATION MISMATCH"). A poor fit of the pattern (`MIN_FIT`, 0.30) is reported: then give the corners of the array (`ARRAY_CORNERS`). On 52 photos (side and number of organoids known): the array was right in all 52 (fit 0.48 to 0.70; the earlier search, from the big circle of the well and the 7 straight edges of the chip, failed in 11, mostly distant photos or photos where the circle or the organoid was misread, and took in the strip or a band in a few others) and the side of the text was right in all 52 without giving it (the landmarks alone: 34 right, 2 wrong, none visible in 16). The organoids are the pale (unsaturated) regions on the array (at least 4% of its area), without their parts narrower than `MIN_WIDTH` (8 electrodes: pale strips along the edges of the array, where the package or the bands look pale, which added the electrodes along an edge to the crop); touching organoids are split when the border between them is a narrow neck. `N_ORGANOIDS` is optional: when given, the largest regions are kept, and the regions are split (a region holding all of them: Gaussian mixture of N parts) until there are N. The electrodes kept are those under the organoid or within `MARGIN` (0.10, as `CROPPING_SAFETY` = 1.10) x its radius of its outline, but not within `SHARED_GAP` (2) electrodes of another organoid or of the border two touching organoids share (conservative on that border; the margin applies on all the other edges). `SelectElectrodes` uses it as option 3, its default (`file_photo`; `photo_organoid`: number of the organoid in the figure, from the top down; optional `photo_text_side` and `photo_n_organoids`), and saves two figures in the folder of the recording: `<file_data>_<organoid_label>_crop_map.png` (electrode map: electrodes with active units kept or not, outline of the organoid, electrodes kept) and `<file_data>_<organoid_label>_crop_photo.png` (copy of the photo with the same outline). `run_pipeline` reads the columns `picture`, `picture_organoid`, `picture_3D_side` and `picture_n_organoids` (`crop_ellipse` or `crop_excel`, when given, are used instead), and does not start a run with `SelectElectrodes` while a recording with a picture has no `picture_organoid`; `make_recordings_table` and the window `CellCircuit` fill `picture` with the photo found in the folder of each recording (`recording_picture`; none when there are several), and leave out the datasets written by `SelectElectrodes` (`cropped_datasets`: `<recording>_<organoid_label>` next to `<recording>`, whose `_QCmanual.mat` holds `spikeTimestamp_BeforeCropping`; the cropping is run again from the row of the recording), and the window shows the organoids of the picture of a row (button *Organoids in the picture...*). Checked on 4 recordings whose ellipses had been set by hand: the same units were kept (5/5, 5/5, 2/2 and 18/18), except one unit on the border of two touching organoids (17 of 18, left out on purpose); and on 52 photos: without the number, 1 organoid was always found right (29 of 29), 2 touching organoids in 7 of 9, and 3 or 4 touching organoids were not separated (give the number). HEIC photos taken with the phone held sideways are turned as they are shown (property `irot` of the file).
- **Folder `legacy`:** the functions that no stage of the pipeline uses were moved there: `Fix_Analysis_Spiketrain`, `PopulationSilence`, `StructureNetworkBurst`, `bootstrap_optimal_k_new`, `build_distance`, `cluster_neuronal_events_new`, `exploreHDF5`, `hartigan_dip_test`, `makeVideo_MEA`, `organized_Finalresults`, `permutation_test_clustering_new`, `plotMEArecording` (and, before, `Brain3_processedData`). In November 2026, the functions that chose the number of families of network burst events before the revision of `cluster_neuronal_events` (see above): `bootstrap_optimal_k`, `permutation_test_clustering`, `build_mds`, `build_mds_v2` and `safe_linkage`; their parameters were removed from `pipeline_parameters.m` (run by hand, they use the same values, written in their code). They are still on the MATLAB path when it is set with `addpath(genpath(...))`.
- **Neuron scores at baseline (November 2026; user's choices):** the neurons of the organoids do not form classes: their waveform and firing-pattern features form one continuous cloud, no more grouped than cluster-free data of the same shape. Each neuron is therefore described by three continuous scores (`NeuronScores`): waveform score 1 (broader and slower spikes), waveform score 2 (larger second peak and faster recovery; its negative end points toward the stereotyped shape that the quality control removes, so the similarity of each template to that shape is given with it, to be used as a covariate) and the firing-pattern score (higher relative rate, higher CV, lower LV and burst index, all measured outside the network bursts, see `FiringPatternFeatures`). The axes are fixed (computed once on 7,963 somatic units of 58 organoids of all the groups, every organoid weighted equally; user's choice, 13 November 2026), so that the scores of different recordings, groups and projects are comparable and recordings can be added at any time; the version of the axes is given with the scores; the units whose template resembles that shape (correlation > 0.95) and the non-somatic units are not scored. `Analyze_spikeTrain` runs the step after the results of the recording are saved (`NEURON_SCORES` = 1, on by default; 0 = not run): the waveform from the Kilosort templates of the recording folder (`loadKilosortTemplates`, `WaveformFeatures`), the firing pattern from the baseline; everything is saved in the new variable `NeuronScore` of `XXX_Analyzed.mat`. Without templates (e.g., 3Brain units), only the firing-pattern score; a step that fails is reported and the rest goes on. `WAVEFORM_WEIGHT` (2/3, user's choice) sets the weight of the waveform when the scores are combined into one distance (`weightedScores`); the scores themselves do not depend on it. `Compile_MEAproject` writes two new sheets: "Neuron scores" (one row per unit: the three scores, the flags, the 10 core waveform measures and the 4 firing-pattern features) and "Neuron composition" (one row per recording: the units scored, the mean scores and the version of the axes, with a note when they are blank). On the 59 recordings of the reference set, the step reproduces the scores of the analysis that set the axes (8,875 units, differences below 1e-9). Run time: a fraction of a second per recording. Tests: `tests/test_waveform_features.m`, `tests/test_firing_pattern_features.m`, `tests/test_neuron_scores.m`.
- **Comparison of groups by their neurons (November 2026; user's choices):** `CompareNeurons` is run by hand once the groups are compiled, with a copy of `comparison_settings.m` (only the recordings tables and the comparisons must be given; everything else has a default, including a default list of activity measures: bursting yes/no, network-burst rate, participation and duration, firing rate). It always fits the axes of the neuron scores again (`FitNeuronAxes`, every organoid weighted equally) on all the groups listed in `FIT_GROUPS` (default: every group of the tables; groups that only help the fitting can be listed), and refuses to run with fewer than 20 organoids; the units are placed on the new axes in memory (the files of the recordings are not changed; a recording analyzed before `NEURON_SCORES` is scored on the fly from its baseline and templates). Fitting on only the two groups compared (5-17 organoids in the D3/D4 data) moved the group differences by up to 0.48 SD (90th percentile) with sign flips; fitting on all the groups of a dataset moved them by 0.17 at most in 90% of the cases, without flips. Then, for each comparison chosen: composition (mean scores of the organoids by group), activity (for each measure, activity ~ group + composition + covariates, with and without the composition; linear or logistic, mixed when random intercepts are named, e.g. a column "batch" of the recordings table; likelihood-ratio test of the composition, change of the group effect; a note when there are fewer than 10 organoids per term), and roles in the network bursts of the baseline (normal scores within each organoid, mixed model with random intercept and slopes per organoid, Satterthwaite degrees of freedom; all the organoids, each group and their difference). Also new: `AXES_FILE` of `NeuronScores` (the axes of a file written by `FitNeuronAxes`; used by `CompareNeurons`), `NeuronScore.WaveformFeature` kept when the templates are not like those of the axes, and `FiringPatternFeatures` computes again the relative rate of a firing pattern already measured against other reference units. The scores on the fixed axes do not change. On the 58 organoids of the reference set, `FitNeuronAxes` gives back the fixed axes (differences below 4e-10); on the 63 recordings of D3/D4 (51 with units), `CompareNeurons` runs in 41 s. Tests: `tests/test_fit_neuron_axes.m`, `tests/test_compare_neurons.m`.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	ShapeSpikeTrainData	%%% (revised October 2026)



Function that takes the data file coming from the spike sorting algorithm and re-shape it in a standard format for it to be used by the rest of the analysis pipeline

(in principle, it was designed to work along with Kilosort)



INPUT:	- folder:    	path to localize the original data (spike\_times.npy and spike\_clusters.npy)

 	- file\_output:  name of the output file where the variable "spikeTimestamp" will be saved in

PARAMETERS:

 	- unitType\_good:   vector indicating the unit types (classification made by "bombcell") that will be included for downstream analysis



OUTPUT: SAVE in a matlab file called "XXX\_raw.mat", a structure array 'spikeTimestamp' (its length is equal to the number of sorted units), with the following fields

 	.recLength:	length of the recording, in seconds.

 	.nUnits:	total number of identified events or spikes for the unit

 	.timeStamp:	time stamp of all events (nUnits) for the unit, in msec

 	.window:	\[start end] time window for the unit extraction, in sec

 	.position:	position of the electrode assigned to the unit

 	.cluster:	cluster id (unit #; can be used to trace back the waveform data)

 	.t\_step:    	temporal resolution of the recording, in msec



 	- rec\_conditions:     list of the filename for the different conditions that were concatenated for the spike sorting

 	- Indiv\_Cond:         cell with individual "spikeTimestamp" structure arrays, one per recording condition



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	QualityControl\_General		%%% (revised November 2025)



Function to run Quality Control assessment in a MEA recording, based on general features of the activity pattern

It removes clusters based on: 1) too many spikes that violate the refractory period criteria, 2) excessive (tonic) firing level, 3) too low level of firing, 4) responsiveness to TTX treatment, 5) cropping of the electrode map based on a picture

It calls the functions "BurstDetection" and "BurstFeatures"



INPUT:	- folder: 		name of the folder to load (and save) the data

 	- condition:		recording condition to process

 	- option\_plot:		option to plot the results

 	- option\_picture:  	option to crop the electrode map using a picture of the tissue sample

 	- option\_pause:  	=1, wait for a key press before moving on to the next file; =0, do not wait (batch processing)

PARAMETERS:

 	- REFRACTORY\_PERIOD:	absolute refractory period for neuronal firing, in msec (default: 1.0)

 	- REF\_PERIOD\_THRESHOLD: threshold for the percentage of spikes that violate the refractory period criteria; above this threshold, the cluster will be discarded (default: 5%)

 	- REF\_PERIOD\_MINSPIKES: threshold for the absolute number of spikes that violate the refractory period criteria (default: 3); to account for the low activity units

 	- ACTIVITY\_THRESHOLD:   minimum level of spiking frequency for an individual electrode to be considered active, in Hz (default: 0.03 Hz = 1/33 sec)

 	- CV\_THRESHOLD:         threshold for coefficient of variation to separate phasic vs tonic activity (default: 3.5)

 	- LVR\_THRESHOLD:        threshold for local variation ratio to separate phasic vs tonic activity (default: 1.00)

 	- PERCENTILE:           percentile to identify the outliers with too high firing rate, among ONLY the clusters with tonic activity (default: 98)

 	- PERCENTILE\_TOP:       percentile to identify the outliers with too high firing rate, among ALL clusters (default: 99.5)

 	- ISI\_RATIO:            threshold for the ratio of isi during bursting vs non-bursting; the lower, the more permissive (default: 8)

 	- N\_SPIKES\_BURST:       threshold for the proportion of spikes in bursts; if the cluster is mostly bursting, probably is a neuron (default: 0.90)

 	- TTX\_EFFECT:           minimum effect of TTX to be considered; proportion relative to the maximum frequency (0.10 - 0.20)

 	- BURST\_METHOD:         method for the burst detection

 	- MIN\_LENGTH:           minimum duration of a recording condition, in sec (default: 90)



OUTPUT: SAVE in a matlab file called "XXX\_cleaned.mat" the following variables:

 	- a structure array 'spikeTimestamp' (same structure than the output generated by "ShapeSpikeTrainData"), ONLY with good clusters

 	- rec\_conditions:       list of the filename for the different conditions that were concatenated for the spike sorting

 	- Indiv\_Cond:           cell with individual "spikeTimestamp" structure arrays (original, before AC), one per analyzed recording condition

 	- removedSpikes:   	percentage of spikes removed for violating the refractory period; since it refers to the originally set of clusters, it needs to be used together with the variable "spikeTimestamp" generated by 					ShapeSpikeTrainData"

   	- echo\_cluster:	    	it indicates if cluster k was discarded for having too many echo spikes (=1, discarded; =0, good)

 	- silent\_cluster:  	it indicates if cluster k was discarded for having low level of activity (=1, discarded; =0, good)

 	- noisy\_cluster:   	it indicates if cluster k was discarded for having high level of tonic activity (=1, discarded; =0, good)

   	- ttx\_cluster:		it indicates if cluster k showed effect upon TTX treatment

   	- picture\_cluster:     	it indicates if cluster k where positively cropped based on a picture of the setup

 	- discarded\_cluster:   	it indicates if cluster k was finally discarded (=1, discarded; =0, good)

 	- n\_discarded:          total number of discarded clusters (independent of the reason)

 	- n\_cluster0:           total number of originally identified clusters

   	- ellipse\_cropped:      parameters that define the ellipse to be cropped

   	- discarded\_signal:     time stamp of units from discarded clusters



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	QualityControl\_HyperSynch	%%%	(revised December 2025)



Function to run Quality control assessment in a single MEA recording, based on the cross-talk between channels

It analyzes hypersynchronicity by calling the function "HyperSynTrain" and using the output generated by "ShapeSpikeTrainData" and "QualityControl\_General"

It also calls the functions "BurstDetection" and "BurstFeatures"



INPUT:	- folders:		name of the folder to load (and save) the data

 	- file\_data: 		name of the file to load (and save) the data

PARAMETERS:

 	- N\_ITERATION:                  number of iterations for the shuffling (default: 500-1000)

 	- N\_SPIKES\_MINIMUM:		minimum number of hypersynchronized spikes to trigger assessment

 	- N\_SPIKES\_BURST:               threshold for the percentage of hypersynchonized spikes within burst, in % (default: 8%)

 	- N\_SPIKES\_NONBURST:            threshold for the percentage of hypersynchonized spikes outside burst, in % (default: 8%)

 	- N\_SPIKES\_TOTAL:               threshold for the percentage of hypersynchonized spikes in all categories, in % (default: 8%)

 	- FACTOR\_MINORITY:              factor to scale up the threshold when analyzing the minority category % (default: 2)

 	- PERCENTILE\_BURST:             threshold percentile based on proportion of hypersynchonized spikes within burst, in % (default: 95%)

 	- PERCENTILE\_NONBURST:          threshold percentile based on proportion of hypersynchonized spikes outside burst, in % (default: 95%)

 	- PERCENTILE\_TOTAL:             threshold percentile based on proportion of hypersynchonized spikes in all categories, in % (default: 95%)

 	- N\_SPIKES\_BURST\_CUT:		threshold for the percentage of hypersynchonized spikes within burst to make the first cut and proceed with the shuffling, in % (default: 4%)

 	- N\_SPIKES\_NONBURST\_CUT:	threshold for the percentage of hypersynchonized spikes outside burst to make the first cut and proceed with the shuffling, in % (default: 4%)

 	- N\_SPIKES\_TOTAL\_CUT:		threshold for the percentage of hypersynchonized spikes to make the first cut and proceed with the shuffling, in % (default: 4%)

 	- METHOD\_SHUFFLING:		method for the spike-train shuffling

 	- OPTION\_METHOD\_BURST:		method for the bursting detection

 	- REFRACTORY\_PERIOD:		refractory period for the shuffling, in msec (default: 0)



OUTPUT:	SAVE the following variables in the file sprintf('%s\\%s\_QC.mat',folder,file\_data) (revised October 2026: compact file):

 	- hs\_master(i,j):		=true, if the pair of clusters (i,j) has significant hypersynchronicity (final master score); used by "QualityControl\_Manual"

 	- hs\_final.burst / .nonburst / .total(i,j):	=true, if the final score of the pair (i,j) is significant within bursts / outside bursts / for all spikes

 	- hs\_percentage.burst / .nonburst / .total(i,j):	percentage of hypersynchronized spikes within bursts / outside bursts / for all spikes

 	- red\_flag(i,j):		=2, super low level of hypersynchronicity (no further assessment); =1, very low level; =0, the statistical assessment was applied

 	(only the pairs i<j are filled). With the parameter option\_save\_rawQC = 1, it also saves "rawQC" and "scores" (the raw results of "HyperSynTrain" and all the scores, as before; it can take several GB for ~800 clusters)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	QualityControl\_Manual	%%%	(revised March 2026)



Function to run Manual Quality control assessment in a MEA recording

It is the final round in a sequence of QC steps, so we recommend a conservative approach (in particular if there is TTX control)

It uses the output generated by "QualityControl\_HyperSynch" and "QualityControl_General"


INPUT:	- folders:		name of the folder to load (and save) the data

 	- file\_data: 		name of the file to load (and save) the data

    - to_discard:       list of clusters to be discarded

    - option_main:      option to display the data

    - option_saveComplete:    option to save the results and end the manual curation

    - minimum_activity:    threshold for the activity level, in Hz (default: 0.04 Hz -> 1/25sec)

    - cond_low_activity:   indexes of the conditions to analyze the minimum activity level
    
PARAMETERS:

 	- METHOD_BURST:        method for the burst detection (default: 

OUTPUT: SAVE the following variables in the file sprintf('%s\\%s_QCmanual.mat',folder,file_data):
    - a structure array 'spikeTimestamp' (same structure than the output generated by "ShapeSpikeTrainData"), ONLY with good clusters
    - rec_conditions:           list of the filename for the different conditions that were concatenated for the spike sorting, and the time window for each condition (in sec)
    - n_discarded_manual:       total number of clusters discarded by "QualityControl_Manual"
    - n_discarded_HS:           total number of clusters discarded on their Hypersinchronicity 
    - minimum_activity:         threshold for the activity level, in Hz
    - idx_abnormalActivity:     indexes of clusters (manually) discarded because of their abnormal pattern of activity  
    - good_cluster:             it indicates if cluster k was accepted or discarded 
    - n_clusterOriginal:        total number of originally identified clusters (before any quality control process)
    - n_discarded_QCgeneral:    total number of clusters discarded by "QualityControl_General"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%. SelectElectrodes   %%%     (revised April 2026)

Function to select the data from a specific set of electrodes (using a list or a picture)
It uses the output generated by "QualityControl_Manual"

INPUT:	- folders:		name of the folder to load (and save) the data

 	- file\_data: 		name of the file to load (and save) the data

    - organoid_label:   specific label for the subset of selected electrodes

    - option:           =1, electrode map in an excel file (3Brain software); =2, ellipse set by hand (ellipse_cropped, x0 y0 radiusx radiusy
                        in electrodes, enlarged by CROPPING_SAFETY); =3 (November 2026, the default), photo of the chip ("OrganoidsFromPhoto"):
                        file_photo, photo_organoid (number of the organoid in the photo), photo_text_side (optional: side of the "3D" text;
                        "": found), photo_n_organoids (optional: 0, found)

OUTPUT: SAVE the following variables in the file sprintf('%s\\%s_%s_QCmanual.mat',folder,file_data, organoid_label):
    - the same variables than the ones generated by "QualityControl\_Manual"
    - spikeTimestamp_BeforeCropping:    structure array 'spikeTimestamp' but before the cropping (it is the original output generated by "QualityControl\_Manual")
    - option:                           option used for the cropping
    - file_excel:                       excel file used for the cropping
    - ellipse_cropped:                  ellipse used for the cropping [x0 y0 radiusx radiusy] (option 3: equivalent ellipse of the organoid, for reference)
    - (option 3) file_photo, photo_organoid, photo_text_side, photo_n_organoids, photo_orientation (how the side of the "3D" text was set,
      or the mismatch), crop_outline (outline of the organoid, [x y] in electrodes),
      crop_electrodes (electrodes kept, [x y]); and the figures <file_data>_<organoid_label>_crop_map.png and _crop_photo.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	IntegrateDataSets	%%%	(revised October 2026)

Function to integrate the datasets of a recording spike-sorted with Kilosort and/or detected with the 3Brain software (it also runs with a single dataset)
It uses the output generated by "QualityControl\_Manual" or "SelectElectrodes" ("XXX\_QCmanual.mat")

INPUT:	- folder:		name of the folder to load (and save) the data

 	- file\_data\_Kilosort / file\_data\_3Brain:	name of the file of each dataset

 	- include\_Kilosort / include\_3Brain:	=1, include the dataset

 	- unitSelect\_Kilosort / unitSelect\_3Brain:	unit types to keep (bombcell classification / waveform-based classification)

 	- option\_criteria:	criteria to identify the same unit in both datasets (spatial superposition and/or spike train similarity)

OUTPUT: SAVE the file "XXX\_curated\_QCmanual.mat" (a single dataset) or "XXX\_integrated\_QCmanual.mat" (both datasets), with the same variables as "QualityControl\_Manual" plus the information about the integration (n\_cluster\_Kilosort, n\_cluster\_3Brain, ...)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	 Analyze\_spikeTrain	%%%     (revised October 2026)



Function to analyze a single spike-train dataset recorded in a MEA platform. It analyzes the spiking activity, the bursting pattern, and the network bursts.

It uses the output generated by "IntegrateDataSets": "XXX\_curated\_QCmanual.mat" (or "XXX\_integrated\_QCmanual.mat")

It calls the functions "SpikeFeatures", "BurstDetection", "BurstFeatures", "NetworkBurstDetection", "NetworkBurstFeatures" and "NetworkBurstStructure", and (November 2026) "Analyze\_cNE" for the coordinated neuronal ensembles (CNE\_ANALYSIS = 1, default; =0, not run)



INPUT:	- folders / files\_data:	folders and names of the files to analyze (e.g., "XXX\_curated")

 	- condition\_analysis:	list of recording conditions to analyze (=-1, all of them)

 	- window:		analysis time at the end of each recording condition, in sec (=-1, the entire condition)

 	- opt\_rerun:		=1, the results are appended to an existing "XXX\_Analyzed.mat" file (e.g., to analyze other conditions); =0, new file

 	- NETWORK\_BURST\_DETECT\_METHOD:	=1, population firing rate with a fixed bin size; =2, normalized population firing rate with dynamic bin size; =3, method 2 + clustering of the events in families (default)

 	- NETWORK\_BURST\_ANALYSIS\_METHOD:	=1, single family of network bursts; =2, one set of features per family (default)



OUTPUT: SAVE the following variables in the file "XXX\_Analyzed.mat"

 	. spikeTimestamp, n\_clusterOriginal, n\_conditions, n\_cluster\_Kilosort, n\_cluster\_3Brain

	. (November 2026) cNEResults: the coordinated neuronal ensembles of each condition and their comparison with the baseline (see "Analyze\_cNE"; only with CNE\_ANALYSIS = 1)

 	. PerCondition(cond), one per recording condition, with the fields: spikeTrain, SpikeFeature, BurstFeature, BurstDetected, BurstFeature\_Population, TonicFeature\_Population, NetworkBurstDetected, NetworkBurstFeature, NetworkBurstFeature0, NetworkBurstStructure, condition, dose, window, file

 	(NetworkBurstFeature{k}: features of the family k of primary events, with its own participating units; the last one: secondary events. NetworkBurstFeature0: all the primary events together)



Network bursts ("NetworkBurstDetection", method 3): 1) preliminary detection with all the units (normalized and non-normalized population signal; the one with more events is kept); 2) screening of the robustly participating units; 3) detection with those units, with an optimized bin size, of the primary events (level 1: large; level 0: small) and of the secondary events (after each level-1 event); 4) statistical assessment of the participation of each unit (500 shuffles per unit) and final detection; 5) quality control (a level-0 event that coincides with a secondary event is kept only as secondary); 6) clustering of the events in families.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	 Compile\_MEAproject	%%%	(revised October 2026)



Function to compile the results from a project based on data recorded in a MEA platform.

Assessment of the spiking activity, the bursting pattern, and the network bursts. Each recording was previously analyzed using the function "Analyze\_spikeTrain" and the files corresponding to a single group ("XXX\_Analyzed.mat") were collected in a single folder

It calls the function "Compile\_NetBurstResults". The recording conditions that were not analyzed yet are skipped (with a warning)



INPUT:	- folderPath:		path of the project folder

 	- group:		name of the group folder (inside folderPath) with the "XXX\_Analyzed.mat" files



OUTPUT: SAVE the following tables in the file sprintf('%s/%s\_Compiled.mat', folderPath, group)

 	. Firing, Bursting: one row per recording condition

	. (November 2026) cNEs: one row per recording condition, the short list of measures of its coordinated neuronal ensembles (CNE\_COLUMNS; see docs/cNE\_METRICS.md), and a note when they are blank;   cNEsVsBaseline: one row per condition against the baseline, and one for the 2nd half of the baseline against its 1st half (CNE\_COMPARE\_COLUMNS);   cNEsVsPrevious: one row per condition against the condition recorded just before it, when that one is not the baseline. Sheets "cNEs", "cNEs vs baseline" and "cNEs vs previous"; the TTX conditions are left out; every other measure stays in "cNEResults" of each "XXX\_Analyzed.mat" file

 	. NETbursting0: network bursts, all the primary events together;   NETbursting1: only the largest family of primary events;   NETbursting2: one row per family of primary events (all of them include the secondary events; since November 2026, plus one row with the remaining events, which belong to no family: label 0, `Family_Type` = "remaining events (no family)")

 	It also saves the tables in an Excel file with the same name (the first row of each sheet has the units), and the file "group\_Grouped&Compiled.xlsx" with the rows grouped by treatment

 	UNITS USED BY THE NETWORK-BURST METRICS (added November 2026). Sheets "Network Bursting Pooled" (NETbursting0), "Network Bursting Primary" (NETbursting1) and "Network Bursting ALL clusters" (NETbursting2):

 	. Participating units: in NETbursting0, the units participating in the primary events (NetworkBurstDetected.Channels; the units that take part only in secondary events are not included); in NETbursting1 and NETbursting2, the units of the family ("familyParticipants": the participating units that fire in more than 20% of the family's events and in at least 2 of them). Used by Number\_NETburst\_Neurons, Number\_Spikes\_NETburst, NumberNorm\_Spikes\_NETburst (divided by all the spikes of the same units), Participation\_NETburst and its CV (fraction of these units active in each event), Amplitude\_NETburst and its CV (inter-spike intervals within the events), Distance\_NETburst and its CV, and the metrics of "NetworkBurstStructure": tempStructNeuron\_NETburst, Onset\_NETburst, Peak\_NETburst and their standard deviations, Correlation\_NETburst, PeakFire\_NETburst and PeakFireNorm\_NETburst

 	. All the recorded units, only as the denominator: NumberNorm\_NETburst\_Neurons (participating units / all the units), NumberNorm\_SpikesTotal\_NETburst (spikes of the participating units within the events / all the spikes of all the units), NumberNorm\_SecondaryBurst\_Neurons (units participating in the secondary events / all the units)

 	. Timing of the events, from the network-burst detection, whose population signal uses the units participating in the primary or the secondary events (all the families): Duration\_NETburst, Duration20\_NETburst, Duration50\_NETburst and their CVs, Tail\_NETburst and its CV, and, from the times of the peaks, Rate\_NETburst, RateCV\_NETburst, FFibi\_CV2\_NETburst, FFibi\_VarMean\_NETburst, FFcount\_NETburst and MeanCount\_NETburst. In NETbursting1 and NETbursting2 they are computed over the events of the family (the intervals are between consecutive events of the family), but each event was timed on the signal of all the participating units

 	. Secondary events (the same in every row: they are not split by family): Number\_SecondaryBurst\_Neurons and Amplitude\_SecondaryBurst use the units participating in the secondary events (NetworkBurstDetected.Secondary.Channels)

 	. No units (numbers of events): Number\_NETburst\_Events (in NETbursting1 and NETbursting2, the events of the family), Number\_clusters\_Netburst, Clustering\_Quality, NumberNorm\_SecondaryBurst, ParentNorm\_SecondaryBurst; labels: label\_cluster\_Netburst (0 = the remaining events, in no family) and Family\_Type (text, November 2026)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%				Complementary scripts						 %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	SpikeFeatures		%%%



Function to analyze the spiking activity (bursting is analyzed by a different function) by calculating XX different metrics on MEA-based recordings.
It is called by the function "Analyze_spikeTrain"

INPUT:	variable "spikeTimestamp" generated by "QualityControl\_Manual" (or "ShapeSpikeTrainData")

PARAMETERS:
    - N_BIN_ENTROPY:            temporal bin to calculate the entropy and the mutual information, in msec (default: 100)
    - COUNT_WINDOW:             length of the windows to count the spikes (count Fano factor), in sec (default: 60)

OUTPUT:	variable "SpikeFeatures"



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	BurstDetection		%%%	(revised November 2025)



Function to detect bursts in a spike train data, using the method specified



INPUT:	variables

 	- signalSpike:		column array with the time stamp for the spike train

 	- option\_method:	method to be used to detect the burst events



PARAMETERS:

 	- RATIO\_PEAK:           proportion of spikes in the burst that are grouped to look for the peak of the burst (default: 0.5)

        - MIN\_SPIKE\_BURST:      minimum number of spikes within a burst (default: 3)

        - ISI\_START\_BURST:      threshold (maximum) ISI to identify the beginning of a burst, in ms (default: 200 -> 5Hz)

        - ISI\_PEAK\_BURST:       minimum peak ISI to accept a burst, in ms (default: 150 -> 6.7Hz)

        - ISI\_END\_BURST\_INI:    initial threshold (maximum) ISI to identify the end of a burst, in ms (default: 300 -> 3.3Hz)

        - FACTOR\_MATURATION:    factor to increase the ISI\_END\_BURST as the burst has more units (default: 1.05 -> 5% increase per unit, starting from MIN\_SPIKE\_BURST)

        - ISI\_END\_BURST\_MAX:    maximum threshold (maximum) ISI to identify the end of a burst, in ms (default: 500 -> 2Hz)

        - MIN\_IBI:              minimum inter-burst interval to identify two consecutive bursts, in ms (default: 500); from finish to start of the next one

        - MIN\_BURST\_DURATION:   minimum burst duration, in ms (default: 25 -> with 3 spikes, it means 80 Hz)



OUTPUT:	variable "BurstDetection", a structure with the following fields:

 	- BurstDetection.spikes(i,1):		time stamp of spike i, in ms

               		.spikes(i,2):           isi corresponding to the spike i, in ms     \[isi(j)=T(j+1)-T(j)]

               		.spikes(i,3):           # of burst to which the  spike i belongs (=0 if spike i is NOT part of a burst)

 	- BurstDetection.IndivBurst(j).start:   time stamp of the first spike of burst j, in ms

 				      .last:    time stamp of the last spike of burst j, in ms

 				      .peak:    time stamp of the peak of burst j, in ms

 				      .median:  time stamp of the median/middle spike of the burst j, in ms

 				      .isi(i):  isi corresponding to the spike i of the burst j

 	- BurstDetection.parameter:             structure with information about the burst detection method

 	- BurstDetection.nBurst:                total number of bursts





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	ShuffleSpikeTrain	%%%	(revised November 2025)



Function to shuffle a spike train

It calls the functions "BurstDetection" and "BurstFeature"



INPUT:  - signal: 		time stamp for the spike train, in msec

 	- range: 		range to generated the new spike train (range(1)=start; range(2)=end), in msec

 	- t\_step: 		temporal resolution - sampling rate, in msec

 	- n\_shuffle: 		number of iterations

 	- refractory\_period:	absolute refractory period, in msec (usually =0 because I do not care)

 	- BurstDetected0:       output from the function "BurstDetection" corresponding to the signal to be shuffled

 	- BurstFeature0:        output from the function "BurstFeatures" corresponding to the signal to be shuffled

 	- method: 		method for spike distribution (=1, uniform; =2, biphasic \[burst \& non-burst]; =3, shuffling of burst location and spikes)

 	- seed:			(optional) seed of the random number generator: same seed + same inputs = exactly the same shuffled spike trains (revised October 2026)



OUTPUT: - ShuffleSpikeTrain{j}(i,1): 	Cell array with the time stamp of the spike i for newly shuffled spike train, iteration j, in msec

 			      (i,2): 	indicates if the spike i belongs to a burst

 	- shuffle\_bins:		(optional, method 3) the same shuffled spike trains, as bin indexes (used by "shuffleCache")



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	HyperSyncPair		%%%	(revised November 2025)



Function to identify hypersynchronized pairs of spikes from the activity recorded in two channels of an MEA system

It calls the function "BurstDetection"



INPUT:	- signal1: 		time stamp recorded in channel 1, in msec (reference channel)

 	- signal2: 		time stamp recorded in channel 2, in msec (study channel)

 	- t\_step: 		temporal resolution - sampling rate, in ms



PARAMETERS:

 	- JITTER:               threshold to define the hypersynchronicity, relative to t\_step (default: 3.0)

 	- MAX\_JITTER:           threshold to define the hypersynchronicity, in msec (default: 1.0)



OUTPUT: structure HS\_spikes, with the following fields

 	.data{k}: 		cell array, where (i,1): time stamp of the spike i, in channel k

                                    		  (i,2): # of hypersynchornized pair for spike i, channel k (=0 if it is not a HS spike)

                                    		  (i,3): =1 if spike i belongs to a burst in channel k, =0 if it is outside bursts in channel k

 	.threshold: 		threshold used to identify hypersynchronized spikes, in msec

 	.burst / nonburst / total: correponding to units within bursts / outside of bursts / entire recording. Each field has a set of subfield:

 		.nspikes\_norm: proportion of spikes displaying hypersynchronicity

 		.nspikes\_hs: total number of spikes displaying hypersynchronicity, ch1 + ch2

 		.nspikes\_total: total number of spikes, ch1 + ch2



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	HyperSynTrain		%%%	(revised November 2025)



Function to assess the statistical significance on hypersynchronicity - Analysis by comparing two channels.

It calls the functions "HyperSyncPair" and "ShuffleSpikeTrain"



INPUT:	- signal1: 		time stamp recorded in channel 1, in msec (reference channel)

 	- signal2: 		time stamp recorded in channel 2, in msec (study channel)

 	- t\_step: 		temporal resolution - sampling rate, in ms

 	- n\_iteration: 		number of iteration for the shuffling

 	- method: 		method for the shuffling of the spike train

 	- range: 		range to generated the new spike train (range(1)=start; range(2)=end), in msec

 	- refractory\_period: 	absolute refractory period, in msec

 	- BurstDetected2:       output from the function "BurstDetection" corresponding to signal2

 	- BurstFeature2:        output from the function "BurstFeatures" corresponding to signal2

 

OUTPUT:	structure QC, with the following fields

 	.threshold: 		threshold used to identify hypersynchronized spikes, in msec

 	.HS\_spikes\_Indv		result for hypersynchronicity assessment done by "HyperSyncPair" for the multiple iterations

 	.burst / nonburst / total: correponding to units within bursts / outside of bursts / entire recording. Each field has a set of subfields:

 		.nspikes\_hs: proportion of spikes displaying hypersynchronicity

 		.nspikes\_total: total number of spikes, ch1 + ch2

 		.percentile: percentile/chance for the measured hypersynchronicity (%)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	BurstFeatures		%%%	(revised October 2026)



Function to analyze the bursting activity (simple spiking is analyzed by a different function, "SpikeFeatures") by calculating different metrics on MEA-based recordings



INPUT:	variable "BurstDetection" (generated by the function "BurstDetection"), a structure with the following fields:

 	- BurstDetection.spikes(i,1):		time stamp of spike i, in ms

               		.spikes(i,2):           isi corresponding to the spike i, in ms     \[isi(j)=T(j+1)-T(j)]

               		.spikes(i,3):           # of burst which the spike i belongs to (=0 if spike i is NOT part of a burst)

 	- BurstDetection.IndivBurst(j).start:   time stamp of the first spike of burst j, in ms

 				      .last:    time stamp of the last spike of burst j, in ms

 				      .peak:    time stamp of the peak of burst j, in ms

 				      .median:  time stamp of the median/middle spike of the burst j, in ms

 				      .isi(i):  isi corresponding to the spike i of the burst j

 	- BurstDetection.nBurst:                total number of bursts

 	- tRange (optional):			[start end] of the analysis window, in sec, to count the bursts in windows (without it, the count metrics are NaN)



PARAMETERS: COUNT_WINDOW, length of the windows to count the bursts (count Fano factor), in sec (default: 60)



OUTPUT: structure "BurstFeatures" with the following fields:

 	- BurstFeatures.nBurst:						total number of burst

 	- BurstFeatures.neuronMetric.nonBursting.isi:          		characterization of the inter-burst interval outside of bursts, in ms (fields: mean, median, cv, and the Fano factor of the intervals as FF_cv2 = cv^2 and FF_varmean = variance/mean)

 	- BurstFeatures.neuronMetric.nonBursting.nSpikes\_norm:      	proportion of non-bursting spikes (ratio of spikes non-participating in bursting)

 	- BurstFeatures.neuronMetric.nSpike:	                   	characterization of the number of spikes per burst (fields: mean, median and coefficient of variation)

 	- BurstFeatures.neuronMetric.nSpikes.norm:                   	total number of spikes within burst normalized by the total number of spikes in the recording (ratio of spikes participating in bursting)

 	- BurstFeatures.neuronMetric.duration:                  	characterization of the burst duration, in sec (fields: mean, median and coefficient of variation)

 	- BurstFeatures.neuronMetric.isi:                       	characterization of the isi within bursts, in ms (fields: mean, median, cv, and the Fano factor of the intervals as FF_cv2 = cv^2 and FF_varmean = variance/mean)

 	- BurstFeatures.neuronMetric.nSpikes.mean:       		mean number of spikes per burst

 	- BurstFeatures.neuronMetric.nSpikes.cv:         		coefficient of variation for the number of spikes per burst

 	- BurstFeatures.neuronMetric.nSpikes.norm:       		total number of spikes within burst normalized by the total number of spikes in the recording (ratio of spikes participating in bursting)

 	- BurstFeatures.neuronMetric.duration:      			characterization of the burst duration, in sec (fields: mean and coefficient of variation)

 	- BurstFeatures.neuronMetric.ibi:           			characterization of the inter-burst interval, in sec (fields: mean, median, cv, and the Fano factor of the intervals as FF_cv2 = cv^2 and FF_varmean = variance/mean)

 	- BurstFeatures.neuronMetric.count:          		bursts counted in windows of COUNT_WINDOW sec (fields: FFcount = variance/mean of the counts, meanCount = mean number of bursts per window, nWindows)

 	- BurstFeatures.indivMetric(i).nSpikes:          		total number of spikes in burst i

 	- BurstFeatures.indivMetric(i).duration:         		duration of the burst i, in sec

 	- BurstFeatures.indivMetric(i).isi:         			characterization of the isi for the burst i, in ms (fields: mean, median, cv, and the Fano factor of the intervals as FF_cv2 = cv^2 and FF_varmean = variance/mean)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	Reshape_spikeTimestamp		%%%

Function to reshape the structure "spikeTimestamp" for the specified temporal window (in sec)

INPUT:  - structure 'spikeTimestamp' with the original data
    - window:       temporal window to extract and reshape the structure 'spikeTimestamp'

OUTPUT: - newly reshaped structure 'spikeTimestamp' 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	plotMEArecording		%%%	(in the folder "legacy", October 2026)



Function to plot data from MEA recording



INPUT:	- folder:    	path to localize the data file to plot

 	- file\_output:  name of the data file to plot (after "Manual QC")

 	- option:	select one out of multiple forms of plotting (and set the specific parameters)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	cNEDetection		%%%	(new, October 2026; robust definition, November 2026; run by the pipeline through `Analyze_cNE`, November 2026)

Function to detect coordinated neuronal ensembles (cNEs), groups of units that repeatedly fire together within the same 15-ms bin. Reimplementation of the method of See et al. 2018 (eLife 7:e35587), based on Lopes-dos-Santos et al. 2013 (J Neurosci Methods 220:149-166), whose code was not shared: z-scored spike counts -> number of cNEs = eigenvalues of the correlation matrix above a bound -> fastICA on the significant principal components (IC weights) -> cNE activity (projection of each IC onto the z-scored counts). Membership (null mean +/- 1.5 sd of the IC weights) and activity (99.9th percentile) thresholds come from 100 matrices in which each unit's binned spike train is circularly shifted.

Robust definition (November 2026, default parameters), so that the cNEs are reproducible, also on large recordings (up to 1059 units tested). On 17 recording-conditions, the cNEs of the odd and of the even 30-s blocks of the window were found in both halves in 86% of the cases and were nearly identical in 68% (63% and 23% with the steps of See et al.):

 	- units: only the units that fire at least 0.5 Hz in the window (MIN_RATE). The rarely firing units made spurious components (a single spike gives them a very large z-score)

 	- number of cNEs: eigenvalues above the 99th percentile of the largest eigenvalue of 100 circularly shifted matrices (N_CNE_RULE "shift"), instead of the Marcenko-Pastur bound ("MP"), which let through 3-7 times more components. Known limitation: a group of units that share slow rate changes (e.g., seconds-long high-rate states) raises this bound for every component, which can hide a weaker ensemble (synthetic example in the test)

 	- consensus: the detection is run on 10 resamples of the window (N_RESAMPLE), each made of 80% of its 30-s blocks; the components with at least 3 members are clustered by weight similarity (cosine, cut at 0.5), and a cNE is a cluster found in at least 7 of the 10 resamples, represented by its centrotype (weights and members of that component). cNE.stability = fraction of the resamples where it was found. N_RESAMPLE = 0: a single run on the whole window

 	- activation: a cNE is active in a bin when (1) its activity is above the 99.9th percentile of its null activity (See et al.) and (2) more of its members fire in the bin than expected by chance, with at least 2 members firing (MEMBER_COUNT_RULE). Rule (1) alone was crossed in bins where 0 or 1 member fired: the activity squares the projection, so many non-members with small negative weights firing together also raise it

 	- chance for the rule (2): interval jitter (COUNT_NULL "jitter"). The bins are grouped in 5-s windows (JITTER_WINDOW), and each spike of a member is moved to a random bin of its window, so the rate changes slower than 5 s are kept; a bin passes when the probability of at least the observed number of members firing is below 0.001 (computed exactly, without surrogates). The first version compared the count with circularly shifted trains (COUNT_NULL "shift"; 99.9th percentile of 200 surrogates), which also removes the slow rate changes: cNEs whose members fire at high rates together for seconds (e.g., after network bursts) were "active" in most bins of those periods. The window sets which rate changes count as a shared state (slower than the window) and which as coordination (faster); network bursts count as coordination with any window of 2-5 s. The window favors the activations inside network bursts, which the cNEs are meant to describe: with 5-s windows, 35 of the 37 cNEs of the recordings with large network events stayed active in as many of them as with the first version, and the activations of the cNEs above more than 1 s from any large event fell by 62-94%. 2-s windows removed more of those (87-98%), but the small cNEs (3-7 members) lost up to all of their large events; windows of 0.5 s or less kept only 1-3 bins per event

The steps of See et al. alone: struct('MIN_RATE', 0, 'N_CNE_RULE', "MP", 'N_RESAMPLE', 0, 'MIN_MEMBERS', 0, 'MEMBER_COUNT_RULE', false, 'BIN_SIZE', 10). Checked by `tests/test_cNEDetection.m` (planted ensembles; a dense implementation of the same method, also of the shifted-data bound and of both versions of the member-count rule; each resample = a single run on its blocks; synthetic slow co-modulation).

Tracking (November 2026; TEMPLATE = the output of `cNEDetection` for another window of the same units, e.g., the baseline before a drug): no new cNEs are detected; the cNEs of the template keep their weights and members, and their activation is computed on this window with the activation settings of the template (bin size, both rules, their percentiles and nulls) applied to this window's counts and spike trains. The units analyzed are those of the template that fire in this window, also below MIN_RATE (a drop of the rate is part of a change); a member silent in the window cannot fire. Tracking on the window of the detection gives its activations exactly (checked by `tests/test_cNECompare.m`). Used by `cNECompare`. A window where fewer than 2 units are analyzed (e.g., a silenced network) now gives no cNEs instead of an error.

Chance level of the activations (November 2026, N_CHANCE = 20): both activation rules are applied, with the thresholds of the real data, to 20 surrogates in which every spike of every unit analyzed is moved to a random bin of its 5-s jitter window (the null of rule (2)); cNE.chanceEvents(k) = mean number of active bins of the cNE k on the surrogates. The rules test each bin at 0.1%, so a cNE whose members rarely fire together can have a large part of its activations at chance level; see `cNEAnalysis` (chance share, excess rate, flag). It does not change the cNEs or their activations (own random generator), and adds little time.

INPUT:	- spikeTimestamp:	structure array with the field .timeStamp (in msec), or a cell array of spike times (in msec)

 	- tRange:		[start end] of the analysis window, in sec (with the consensus, at least 3 blocks of 30 sec)

 	- inputs:		(optional) structure that overwrites the PARAMETERS (BIN_SIZE, MIN_RATE, N_CNE_RULE, N_SHIFT_EIG, EIG_PERCENTILE, N_RESAMPLE, RESAMPLE_BLOCK, RESAMPLE_FRACTION, MIN_RESAMPLES, CLUSTER_SIMILARITY, MIN_MEMBERS, N_SHUFFLE, MEMBER_SD, ACTIVITY_PERCENTILE, MEMBER_COUNT_RULE, COUNT_NULL, JITTER_WINDOW, COUNT_PERCENTILE, N_SHIFT_COUNT, MIN_FIRING, ICA_CONTRAST, ICA_MAX_ITER, ICA_TOL, SEED, N_CHANCE); TEMPLATE for the tracking

OUTPUT:	- cNE.nCNE, cNE.weights(i,k):		number of cNEs, and IC weight of the unit i in the cNE k

 	- cNE.members{k}, cNE.negativeMembers{k}:	units with a weight above / below the membership thresholds (cNE.memberThreshold(k,:))

 	- cNE.stability(k), cNE.similarity(k):	fraction of the resamples where the cNE k was found, and mean weight similarity within its cluster

 	- cNE.activity(k,b), cNE.activityThreshold(k):	activity of the cNE k in the bin b, and the threshold of the rule (1)

 	- cNE.membersNeeded(k):			members that must fire in a bin for the rule (2) with COUNT_NULL "shift" (NaN with "jitter", where it changes from bin to bin)

 	- cNE.events{k}, cNE.eventTimes{k}:		bins in which the cNE k is active (both rules), and their midpoint (in msec); cNE.activityEvents{k} and cNE.memberCountEvents{k}: bins that pass each rule

 	- cNE.chanceEvents(k):			active bins of the cNE k expected by chance (interval-jitter surrogates; NaN with N_CHANCE = 0)

 	- cNE.rate, cNE.active:			firing rate of each unit in the window, and the units analyzed

 	- cNE.resamples(r):			blocks, number of components, bound, weights and members of each resample

 	- cNE.eigenvalues, cNE.lambdaMax, cNE.lambdaMP, cNE.binTime, cNE.nullWeight, cNE.icaConverged, cNE.parameters


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	cNEAnalysis		%%%	(new, November 2026; run by the pipeline through `Analyze_cNE`)

Function to describe the network activity of one recording-condition with its cNEs (output of `cNEDetection`): when each cNE is active relative to the network events of option 5 of `NetworkBurstDetection`, in which order the cNEs are active within the large events, and which cNEs the small and the secondary events recruit. The cNEs do not detect network events; they extend their analysis.

 	- periods: each 15-ms bin of the cNE detection belongs to one period, by precedence: secondary (secondary events) > core (large events, start to end) > small (small events) > tail (end of a large event to the end of its closed zone) > onset (the 0.5 s before the start of a large event, ONSET_WINDOW) > quiet baseline. Per cNE: share of its activations and activations per minute in each period

 	- time of a cNE in an event: 25th percentile of the times of its active bins, from the start of the event (with fewer than 4 active bins: midpoint of the first one and the median). Windows: large events from 0.5 s before their start (after the previous tail) to the end of their tail (before the next onset); small events from 0.5 s before their start (outside the large events, onsets and tails) to their end; secondary events from start to end. Template = median time of each cNE over the large events (rank, and position from 0 = first to 1 = last)

 	- order of the large events: consistency S = sum_ij n_ij |2 p_ij - 1| / sum_ij n_ij over the pairs of cNEs, against a null that permutes the times among the active cNEs of each event (2000 permutations), reported scaled to chance, (S - S0)/(1 - S0) (0 = chance, 1 = the same order in every event), with its p; per-event match = Kendall tau between each event and the template of the other events (median); span of the template (sec, and as a fraction of the duration of the large events); timing of the order (November 2026): for each pair of cNEs active together in at least 3 large events, the lag between them in each of these events (later minus earlier cNE of the template), its median (summary.orderLag) and its median absolute deviation (summary.orderJitter: how precise the timing is), in sec, medians over the pairs

 	- small and secondary events: % with an active cNE, mean % of the cNEs an event recruits, over the events that recruit at least one (November 2026), mean template position of the cNEs they recruit (0 = leading, 1 = late), agreement with the template order (share of the pairs of active cNEs in that order; 0.5 = chance), with its p; also the consistency of their own order, with its p (small events: November 2026, with its own generator of permutations, so that the other p values do not change)

 	- member participation: share of the members of a cNE that fire in each of its active bins (median; also inside and outside the large events). It falls with the number of members, so it is meant to compare a cNE with itself (e.g., across conditions)

 	- roles (a cNE can have several): core (largest share in the core), leader (active in the onset of at least 25% of the large events), tail (at least 30% of its activations in the tails and active in at least 50% of them), reverberation (active in at least 50% of the secondary events), small-event (active in at least 50% of the small events), baseline (largest share in the quiet baseline); the reverberation and small-event roles need at least 3 events of that type

 	- minimums: the order measures need at least 2 cNEs (November 2026; was 3) and 5 large events (with 2 cNEs, the order is which cNE leads: the consistency tells whether the same one leads in every event, and the per-event match stays blank, as it needs 3 active cNEs in an event); the comparisons of the small and secondary events, at least 3 events with an active cNE (position) or with 2 active cNEs (agreement, own order)

 	- chance level (November 2026): activations expected by chance (cNE.chanceEvents of `cNEDetection`) as a rate and as a share of the activations; excess rate = activation rate - chance rate; a cNE is flagged at chance (atChance) when the chance level is at least 50% of its activations (MAX_CHANCE_SHARE). Nothing is left out: the flag tells which activation counts are mostly what the rules let through by chance

 	- reference order (REFERENCE_TIME, for `cNECompare`): with the template of the same cNEs in a reference condition, the large events of this condition are compared with that order: agreement A (as for the small events; at least 2 cNEs and 3 events with 2 active cNEs that have a reference time) with its p, and Kendall tau-b between the two templates (at least 3 cNEs with both times; summary.referenceAgreement, referenceAgreementP, referenceMatch)

 	- events the cNEs do not describe: when fewer than 50% of the large events have an active cNE (summary.notDescribed), the roles, the times of the cNEs and the order and event-type measures are left blank; the periods and the share of the events where each cNE is active are kept. No large event in the window (summary.noBursting): only the cNEs, their activation rates and member participation. Only the events entirely inside the window of the cNE detection are counted

Checked by `tests/test_cNEAnalysis.m` (synthetic activations and network events with known answers) and, on 17 recording-conditions, against the separate scripts with which the measures were designed (same periods, times per event, templates, onset leads, member participation and roles).

INPUT:	- cNE:			output of `cNEDetection` for the condition

 	- NetworkBurstDetected:	output of `NetworkBurstDetection` (option 5) for the same condition ([] or nBurst = 0: no network events)

 	- spikeTimestamp:	(optional) spike trains given to `cNEDetection`, for the member participation ([]: not computed)

 	- inputs:		(optional) structure that overwrites the PARAMETERS (ONSET_WINDOW, TIME_PERCENTILE, MIN_BINS_TIME, N_PERM, SEED, MIN_COVERAGE, MIN_CNE_ORDER, MIN_LARGE_ORDER, MIN_EVENTS_TYPE, MIN_EVENTS_ROLE, LEADER_ONSETS, TAIL_SHARE, TAIL_COVER, REVERBERATION_COVER, SMALL_COVER, REFERENCE_TIME, MAX_CHANCE_SHARE)

OUTPUT:	- report.summary:		one row per condition: units and cNEs, members, % of units in a cNE, activation rate, excess rate and number of cNEs at chance, member participation; number of events of each type, % of the large events with an active cNE and the flags; % of time and median share of the activations in each period; number of cNEs per role and the lead of the leaders; order consistency (scaled, p), per-event match, span, lag and jitter between the cNEs; small and secondary events (% with a cNE, mean % of the cNEs recruited by the events with a cNE, position, agreement, own order, and their p); with REFERENCE_TIME, the large events against the reference order

 	- report.cNEs:			table, one row per cNE: members (and their cluster ids, to track the cNE across conditions), stability, activation rate, chance rate, excess rate, chance share and flag, member participation, share and rate in each period, % of the large events, onsets, tails, small and secondary events where it is active, onset lead, template time, rank, position, roles

 	- report.events:		per event (large, small, secondary): limits, time of each cNE, active bins; large events, also the per-event match and the share of the events where each cNE is active at each lag from their start

 	- report.binPeriod, report.periods, report.parameters:	period of each bin of the cNE detection, names of the periods, parameters used

 	- fig:			(only if requested) shares of the periods with the roles, activity of each cNE around the large events, and participation in each type of event


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	cNECompare		%%%	(new, November 2026; run by the pipeline through `Analyze_cNE`)

Function to compare the cNEs of two conditions of the same recording, e.g., the baseline and a drug, or the two halves of the baseline (the change that time alone brings, to compare a drug with). Two complementary views:

 	- tracking: the cNEs of the reference condition are followed in the test condition (same weights and members; activation on the test condition's own data, `cNEDetection` with TEMPLATE), and both conditions are described with `cNEAnalysis` (each with its own network events; the test condition also against the order of the reference). Per cNE: activation rate, mean firing rate of its members (if the activations fall more than the members' rates, the members fire less together, not only less), member participation, share of the large events where it is active, roles and position in the order, in both conditions

 	- separate detection: the cNEs of the test condition are detected on their own and matched one to one with the reference cNEs by weight similarity (cosine; greedy from the most similar pair; at least 0.5, MATCH_SIMILARITY, the cut of the consensus of `cNEDetection`): reference cNEs found again, lost, and new cNEs

Ratios are test/reference (NaN when the reference value is 0). The excess rates over chance are compared too (excess ratio; the excess of the test condition counted as 0 when negative), and the reference cNEs whose activations are mostly at chance are flagged and left out of the median excess ratio. Checked by `tests/test_cNECompare.m` (synthetic ensembles kept, uncoordinated at the same rates, silenced and new; no change with large network events in a fixed order, and the order reversed; silenced test window).

INPUT:	- reference, test:	one structure per condition: spikeTimestamp (the same units, in the same order, in both; checked with .cluster), tRange ([start end], in sec), and optionally NetworkBurstDetected (option 5), name, and cNE (a detection of that window, to reuse it)

 	- inputs:		(optional) structure that overwrites the PARAMETERS (MATCH_SIMILARITY; DETECTION and ANALYSIS: inputs passed to `cNEDetection` and `cNEAnalysis`)

OUTPUT:	- comparison.summary:		reference cNEs and cNEs detected in the test condition; found, lost, new, % found, median best similarity; % of the reference cNEs active in the test condition; median ratios of the activation rates, of the excess rates (cNEs not at chance) and of the member firing rates; reference cNEs at chance; median change of member participation and of the share of the large events; roles kept, lost and gained (NaN when a condition has no roles); the large events of the test condition against the reference order (agreement, p, template match); large events and flags of both conditions

 	- comparison.cNEs:		table, one row per reference cNE: members, cluster ids, silent members in the test window, activation rate, chance share, excess rate (and their ratios), flag at chance, member firing rate, member participation, share of the large events, position and roles in both conditions, ratios, matched test cNE and its similarity, best similarity

 	- comparison.similarity, comparison.newCNEs:	weight similarity of each reference cNE with each test cNE; test cNEs not matched

 	- comparison.reference, comparison.tracked, comparison.detected, comparison.test, comparison.parameters:	the cNEs and `cNEAnalysis` reports of the reference, of the tracked cNEs in the test window and of the separate detection; names and windows; parameters used

 	- fig:			(only if requested) activation rate and member firing rate of each reference cNE in both conditions, and the weight similarity with the cNEs detected in the test condition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	Analyze\_cNE		%%%	(new, November 2026)

Function of the pipeline, called by `Analyze_spikeTrain` (CNE\_ANALYSIS = 1): the coordinated neuronal ensembles (cNEs) of each condition of a recording, and their comparison with the baseline and with the condition recorded before. 1) Each condition: its cNEs on its own analysis window (`cNEDetection`) and their description against its network events (`cNEAnalysis`). 2) Each other condition against the baseline (`cNECompare`; the baseline cNEs of 1, detected once, followed in the condition, and the cNEs of the condition, also from 1, matched to them). 3) Each condition against the condition recorded just before it, when that one is not the baseline (`previousCondition`: the condition whose analysis window starts last before the start of its window; e.g., drug 3 against drug 2), in the same way. 4) The 2nd half of the baseline against its 1st half (`cNECompare`): the change that time alone brings, to compare the drugs with. The TTX conditions (label with "TTX") are left out of every step (status "excluded (TTX)"; user's choice). In a rerun, the results of the conditions that were not analyzed again are kept (and computed if they are missing), and only the comparisons that involve the analyzed conditions are computed again (all of them when the baseline is analyzed again, or is another condition). A step that fails is reported with a warning and its status ("error: <message>"); the other steps go on.

INPUT:	- PerCondition, base\_cond:	results of `Analyze_spikeTrain` for the recording, and the index of the baseline

 	- analyzed:		conditions analyzed in this run

 	- cNEResults:		(optional) results of a previous run ([] = none)

 	- options:		(optional) plots, folder, file\_data, hide (figures, as the other plots of `Analyze_spikeTrain`); detection, analysis, compare (inputs of `cNEDetection`, `cNEAnalysis` and `cNECompare`)

OUTPUT:	- cNEResults.condition(c):	name, window, status ("ok", "too few units", "not analyzed", "error: ..."), cNE (`cNEDetection`), report (`cNEAnalysis`)

 	- cNEResults.comparison(c):	each condition against the baseline (empty for the baseline): name, status, result (`cNECompare`)

 	- cNEResults.previous(c):	each condition against the condition recorded just before it: reference (index of that condition; 0 = no such comparison), name, status, result (`cNECompare`)

 	- cNEResults.halves:		the 2nd half of the baseline against its 1st half: name, status, result (`cNECompare`)

 	- cNEResults.baseline, cNEResults.created:	index of the baseline; date and time of the run

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	NeuronScores		%%%	(new, November 2026)

Function of the pipeline, called by `Analyze_spikeTrain` (NEURON\_SCORES = 1) for the baseline condition: three continuous scores per neuron, on fixed axes. Waveform score 1 = first principal component of the 10 core measures of `WaveformFeatures` (44% of their variance; + broader and slower spikes); waveform score 2 = second component (29%; + larger second peak and faster recovery, - larger first peak; its negative end points toward the QC noise shape, so it is used with the similarity to that shape as a covariate, and its results are repeated on the "clean" units, user's choice); firing-pattern score = first component of the 4 features of `FiringPatternFeatures` (47%; + higher relative rate, higher CV, lower LV and burst index). The axes were computed once (13 November 2026) on 7,963 somatic units of 58 organoids without the look-alikes (waveform) and on the 5,180 of them with the 4 features (firing pattern), every organoid weighted equally (each of its n units counts 1/n): each measure is transformed (log of the skewed ones), standardized with the weighted means and SDs of these units and projected on the principal axes of their weighted correlation matrix, so each score is in SD of the reference units (mean 0). The axes must be fitted on a broad set that includes every group to be compared, never on one of these groups: fitted on a single group of 5-10 organoids, they move the differences between groups by up to 0.3-0.65 SD of the organoid means. Their version (REFERENCE\_VERSION, in NeuronScore.reference.version and in the sheet "Neuron composition") tells which axes the scores are on. The reference values are in the REFERENCE section of the function (they are not parameters). The waveform axes only apply to templates like those of the reference (Kilosort 2.5 export, 61 samples after the padding, about 19.75 kHz). (New, November 2026) AXES\_FILE points the scores to other axes, written by `FitNeuronAxes` (used by `CompareNeurons`, which fits the axes again on the groups it compares; default "": the fixed axes); an axes file whose templates are not like the QC noise shape has no noise shape (no look-alikes, no "clean" units). Units scored: somatic, not look-alikes (similarity > LOOKALIKE\_SIMILARITY) and with the 10 measures; the firing-pattern score for these units when they have the 4 features. Without usable templates, the firing-pattern score of every unit with the 4 features (relative rate against all units).

INPUT:	- spikeTrain, tRange, netBursts:	baseline units, window and network bursts (as `FiringPatternFeatures`), or the output of `FiringPatternFeatures`

 	- waveforms:		templates (as `loadKilosortTemplates`) or the output of `WaveformFeatures`; [] without templates

 	- samplingRate, unitType:	as `WaveformFeatures`

OUTPUT:	- NeuronScore.waveformScore1, waveformScore2, firingPatternScore:	the scores of each unit (NaN when not scored)

 	- NeuronScore.weightedScores:	the three scores scaled for a combined distance with WAVEFORM\_WEIGHT (the waveform gets 2*WAVEFORM\_WEIGHT of variance along its first axis, the firing pattern 2*(1-WAVEFORM\_WEIGHT))

 	- NeuronScore.noiseSimilarity, lookAlike, clean, nonSomatic, waveformScored, firingScored, templatesUsed:	information about each unit

 	- NeuronScore.composition:	nUnits, nSomatic, nLookAlike, nWaveformScored, nFiringScored, nClean and the mean scores of the scored units (NaN with less than MIN\_COMPOSITION\_UNITS of them)

 	- NeuronScore.WaveformFeature, FiringPattern, scores, reference (with its version and file), parameters (with `Analyze_spikeTrain`, also cluster and condition); WaveformFeature is also kept when the templates are not like those of the axes (revised November 2026)

PARAMETERS:	WAVEFORM\_WEIGHT = 2/3 (user's choice), LOOKALIKE\_SIMILARITY = 0.95, CLEAN\_SIMILARITY = 0.7, MIN\_COMPOSITION\_UNITS = 10, AXES\_FILE = "" (new, November 2026). Tests: `tests/test_neuron_scores.m`, `tests/test_fit_neuron_axes.m`.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	WaveformFeatures		%%%	(new, November 2026; run by the pipeline through `NeuronScores`)

Function to measure the shape of the spike waveform of each unit on its Kilosort 2.5 template (peak electrode; whitened and of unit norm, so only the shape is measured, on the template scaled to trough = -1 and upsampled x10 with a cubic spline). The 10 core measures: trough-to-peak, half-width, P2/|T|, P1/|T|, asymmetry (P2-P1)/(P2+P1), steepest descent, steepest recovery, their ratio, the derivative trough-to-peak (CellExplorer) and the area of the phase plot (= (2*pi*RMS frequency)^2); plus 2 spectral measures (spectral centroid, power above 1 kHz), measured but left out of the scores. A waveform whose largest deflection is positive is flipped and flagged as non-somatic, as the units of bombcell type 3.

INPUT:	- waveforms:		template of each unit on its peak electrode (82 samples, the first 21 are zeros), cell array or matrix

 	- samplingRate:		in Hz

 	- unitType:		(optional) bombcell type of each unit

OUTPUT:	- WaveformFeature.<measure>:	the 12 measures (NaN when they cannot be measured, e.g., the second peak after the end of the template)

 	- WaveformFeature.polarity, unitType, nonSomatic, nPeaks, nTroughs, edgePeakBefore, edgePeakAfter, waveform, time, measures, parameters

Test: `tests/test_waveform_features.m`.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	loadKilosortTemplates		%%%	(new, November 2026)

Function to read the Kilosort template of each unit of spikeTimestamp: from "waveforms.mat" of the recording folder (or, if missing, from the column "waveform\_template" of "metadata.csv"), matched by spikeTimestamp.cluster to the column "new\_cluster\_id". Kilosort units only: it stops if a unit is not in "metadata.csv" (e.g., 3Brain units; `Analyze_spikeTrain` then computes only the firing-pattern score).

INPUT:	- folder, spikeTimestamp

OUTPUT:	- templates (one row per unit), unitType (bombcell), samplingRate (Hz, from spikeTimestamp.t\_step)

Test: `tests/test_waveform_features.m`.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	FiringPatternFeatures		%%%	(new, November 2026; run by the pipeline through `NeuronScores`)

Function to measure the firing pattern of each unit outside the primary network bursts (each burst from its first to its last spike, widened by GUARD = 0.5 s on each side): the relative firing rate (log10 of the rate outside the bursts minus the median of the reference units), the local variation (LV), the CV of the inter-spike intervals and the burst index (Royer et al. 2012: autocorrelogram at 3-5 ms against 200-300 ms). Measured on the whole baseline, these features mostly restate the role of each unit in the network bursts (about half of its spikes fall inside them); outside the bursts, LV, CV and the burst index are independent of those roles. No interval or pair of spikes spans a burst. The verdict of the pipeline is followed (user's choice, November 2026): when `NetworkBurstDetected.nBurst` is NaN (the candidate bursts were rejected by the final check of `NetworkBurstDetection`) or 0, the recording has no network bursts and the whole window is used, even if IndivBurst still holds the candidates (before this choice, their windows were left out: 9-70% of the baseline in 8 of 59 recordings).

INPUT:	- spikeTrain, tRange:	units and analysis window (sec) of the condition; (new, November 2026) or a firing pattern already measured: only its relative rate is computed again, against the units of "reference"

 	- netBursts:		NetworkBurstDetected, or [start last] in sec, or [] (the whole window)

 	- reference:		(optional) units that set the median rate of the organoid (all units if omitted)

OUTPUT:	- FiringPattern.relativeRate, lv, cv, burstIndex:	the 4 features (NaN with less than MIN\_SPIKES = 50 spikes outside the bursts)

 	- FiringPattern.rateOutside, nSpikes, nSpikesOutside, enoughSpikes, timeOutside, periods, features, parameters

Test: `tests/test_firing_pattern_features.m`.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	FitNeuronAxes		%%%	(new, November 2026; run by `CompareNeurons`, or by hand)

Function to fit a new set of axes for `NeuronScores` on the `NeuronScore` of many analyzed recordings of ALL the groups to compare (axes fitted on a single group of 5-10 organoids move the differences between groups by up to 0.3-0.65 SD of the organoid means). `CompareNeurons` runs it on all the groups of a comparison before comparing them. The axes are fitted as the fixed ones: the somatic units that are not look-alikes and have the 10 measures (waveform), the same units with the 4 features (firing pattern, relative rate against them; every unit of a recording without templates), every organoid weighted equally, weighted standardization and principal axes of the weighted correlation matrix, scores of weighted mean 0 and SD 1; signs: waveform score 1 with the half-width, score 2 with the asymmetry, the firing-pattern score with the CV. The templates of all the sources must be alike (same number of samples, sampling rates within 1%); if they are not like the QC noise shape, the axes have no noise shape. On the 58 organoids of the reference set it gives back the fixed axes (differences below 4e-10; split halves |cos| 0.98, 0.98 and 0.998).

INPUT:	- sources:		"\_Analyzed.mat" files with NeuronScore (one per organoid), or a cell array of NeuronScore structures

 	- groups:		group of each source

 	- axesFile:		(optional) .mat file to write (variable NeuronAxes; an existing file is not replaced)

 	- version:		(optional) name of the new axes (default: date, organoids and groups); it must differ from the version the sources were scored with

OUTPUT:	- NeuronAxes:		the fields of the REFERENCE section of `NeuronScores`, and NeuronAxes.fit (sources, groups, units and organoids used, all the eigenvalues, stability between halves of the organoids)

PARAMETERS:	MIN\_ORGANOIDS = 20 and MIN\_GROUPS = 2 (otherwise nothing is fitted), MIN\_ORGANOID\_UNITS = 1, N\_SPLITS = 50, MIN\_STABILITY = 0.9 (warning below). Test: `tests/test_fit_neuron_axes.m`.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%	CompareNeurons		%%%	(new, November 2026; run by hand once the groups are compiled)

Function to compare groups of organoids by the neurons they are made of (user's choices, 13 November 2026). The organoids of a group are the "\_Analyzed.mat" files of its folder `<project>/<group>` (as for `Compile_MEAproject`); the activity comes from `<project>/<group>_Compiled.mat`, at the baseline (identified by its file name) or at the condition chosen. Steps: (1) the axes of the scores are fitted again on all the organoids of FIT\_GROUPS (`FitNeuronAxes`; at least 20 organoids of at least 2 groups) and every unit is placed on them in memory (the files are not changed; a recording without saved scores is scored on the fly from its baseline and the Kilosort templates of its folder); (2) for each comparison (2 or more groups, the first one is the reference): composition (mean score ~ group, with the similarity to the QC noise shape for waveform score 2, + random intercepts), activity (activity ~ group + covariates, with and without the composition, + random intercepts; linear, or logistic for a yes/no measure; mixed when random intercepts are named; likelihood-ratio tests of the group and of the composition, change of the group effect in %, composition effect per SD of the organoid means) and roles (onset and participation in the primary network bursts of the baseline, organoids with at least MIN\_BURSTS bursts and MIN\_ROLE\_UNITS units; normal scores within each organoid; mixed model with random intercept and random slopes of the scores per organoid, Satterthwaite degrees of freedom, or organoids - 1 when they cannot be computed; across the organoids, in each group, and the difference between the groups).

INPUT:	- settings:		path of a copy of `comparison_settings.m` (TABLES, PROJECTS, CONDITION, FIT\_GROUPS, COMPARISONS, ACTIVITY, COMPOSITION, COMPOSITION\_MODE, COVARIATES, RANDOM, ROLES, ROLE\_SCORES, ROLE\_COVARIATES, MIN\_BURSTS, MIN\_ROLE\_UNITS, OUTPUT, NAME; see that file), or a structure with these fields; "defaults" returns the defaults

OUTPUT:	- Results:		settings, axes (with their agreement with the fixed axes), Organoids, Composition, Activity, Roles and Units tables; saved in `<OUTPUT>/<NAME>_<date and time>.mat` and `.xlsx` (one sheet per table, with the settings and the axes), the axes in `..._axes.mat`

Default activity measures (ACTIVITY; each row: name, `<table>.<column>` of the compiled tables, transformation none / log / yes/no, organoids all / bursting): bursting (NETbursting1.Number\_NETburst\_Events, yes/no, all), network-burst rate (NETbursting1.Rate\_NETburst, log, bursting), network-burst participation (NETbursting1.Participation\_NETburst, none, bursting), network-burst duration (NETbursting1.Duration\_NETburst, log, bursting), firing rate (Firing.Firing\_Rate, log, all). Default covariate: MeanNoiseSimilarity (also LogUnits, PrimaryBursts, a compiled column or a column of the recordings tables). Test: `tests/test_compare_neurons.m`.
