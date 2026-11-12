# Regression tests (October 2026)

These scripts check that a revised version of the pipeline gives the **same results** as a reference version, and measure how long each stage takes. They never modify the reference recording: every stage runs on a copy of the data, in a fresh MATLAB process.

## Files

| File | What it does |
|---|---|
| `regression_config.m` | Paths and settings: reference recording, npy-matlab folder, temporary folder, number of shuffling iterations and clusters used for `QualityControl_HyperSynch` |
| `make_baseline.m` | Saves a snapshot of the current code as the **old** (reference) version |
| `run_regression.m` | Runs one stage (or `"all"`) with the old and the new code, compares every saved variable and writes a report with the profiler summary |
| `run_stage.m` | Runs one stage with one version of the code (called by `run_regression`) |
| `compare_outputs.m` | Compares all the variables in two `.mat` files |
| `test_core_functions.m` | Function-level test of `BurstDetection`, `BurstFeatures`, `lvr` and `Reshape_spikeTimestamp` on every cluster of a recording plus synthetic edge cases (run it with each version, then compare with `compare_outputs`) |
| `test_spike_features.m` | Function-level test of `SpikeFeatures` (entropy, mutual information, STTC, LV) on the clusters of an analysis window, same procedure |
| `test_hypersynch_functions.m` | Function-level test of `HyperSyncPair` (random pairs + edge cases) and `ShuffleSpikeTrain` (methods 1, 2 and 3), same procedure |
| `test_cNEDetection.m` | Self-check of `cNEDetection` on synthetic spike trains: planted ensembles and events recovered (consensus of resamples, stability), same results as a dense implementation of the same method (also the shifted-data bound, the rate cutoff and the member-count rule with both nulls, interval jitter and circular shifts), each resample = a single run on its blocks, synthetic slow co-modulation (chance-level activations with the interval-jitter null), chance level of the activations (small for planted ensembles, no effect on the cNEs), reproducibility and edge cases (no reference recording needed) |
| `test_cNEAnalysis.m` | Self-check of `cNEAnalysis` (November 2026) on synthetic cNE activations and network events with known answers: period of every bin (against a direct labelling), shares and rates of the periods, onsets and tails, a planted order (consistency 1) and a random one (chance), lag and jitter between the cNEs (planted lags, direct count), small and secondary events vs the template, their own order (mixed order below chance) and the share of the cNEs they recruit (direct count), the large events against a reference order, roles, member participation (planted spikes), events the cNEs do not describe, minimums (1 cNE: no order; 2 cNEs: which one leads), no network bursting, only small events, chance level (shares, flag, excess rate), reproducibility (no reference recording needed) |
| `test_cNECompare.m` | Self-check of `cNECompare` and of the tracking of `cNEDetection` (November 2026) on synthetic spike trains: tracking on the window of the detection = its activations; a "drug" with an ensemble kept, one whose members keep their rates but stop firing together, one silenced and a new one (found / lost / new, activation, excess and member-rate ratios, silent members, chance level); no change with large network events in a fixed order (all found, same rates and roles, agreement with the reference order) and the order reversed; silenced test window, no cNEs, units that differ (no reference recording needed) |
| `test_Analyze_cNE.m` | Self-check of `Analyze_cNE`, the cNE step of `Analyze_spikeTrain` (November 2026), and of `previousCondition`, on a synthetic recording (baseline, a drug that breaks one of two planted ensembles, a second drug, a condition with too few units, TTX): the same results as calling `cNEDetection`, `cNEAnalysis` and `cNECompare` directly (each drug against the baseline, the second drug against the first), the baseline halves, TTX left out, what a rerun keeps and computes again, a step that fails (reported, the rest goes on), and the figures (saved, hidden ones closed and opening visible) (no reference recording needed) |
| `test_family_clustering.m` | Self-check of the families of network burst events (`cluster_neuronal_events`, `primaryFamilies`, November 2026) on synthetic spike trains: participating units of each event (planted units found, none in a long event without ensemble), one family per repeated ensemble, no family for an ensemble seen twice, stability, reproducibility (no reference recording needed) |
| `test_CellCircuit.m` | Self-check of the window `CellCircuit` (November 2026): every parameter of `common_parameters` exists in `pipeline_parameters.m`, the window opens with the common adjustments first, the values typed are checked, the equivalent `run_pipeline` command follows the settings, and missing settings are reported (no reference recording needed) |
| `test_recordings_table.m` | Self-check of the reading of the recordings table by `run_pipeline` (November 2026): a .csv table with file names full of "_" and empty optional columns is read column by column (it used to be read with "_" as the delimiter), as is the same table saved as .xlsx; a table saved by the window for a few recordings of a long list (unticked rows without a group or file name are skipped, a ticked row without a group still stops the run) (no reference recording needed) |
| `test_baseline_condition.m` | Self-check of the baseline found by its file name (`baselineCondition`, used by `Analyze_spikeTrain`, November 2026) on lists of file names, of `Analyze_spikeTrain` on a synthetic recording whose TTX condition is listed first (labels "Baseline" and "+TTX"), and of two conditions with the same label (two TTX files): each keeps its own network-burst plot, with the index of the condition in its name (no reference recording needed) |

## How to use

First, set the reference recording and the npy-matlab folder: edit `cfg.ref_folder`, `cfg.base_name` and `cfg.npy_folder` in `regression_config.m`, or set the environment variables `CELLCIRCUIT_REF_FOLDER`, `CELLCIRCUIT_BASE_NAME` and `CELLCIRCUIT_NPY_FOLDER` (this keeps the paths of your computer out of the repository). Then:

```matlab
addpath('path/to/CellCircuit/tests')
setenv('CELLCIRCUIT_REF_FOLDER', "path/to/reference_recording");          % (or edit regression_config.m)
setenv('CELLCIRCUIT_BASE_NAME', "reference_recording_name");
setenv('CELLCIRCUIT_NPY_FOLDER', "path/to/npy-matlab");
make_baseline();                                  % once, before modifying the code
run_regression("QualityControl_General");         % after modifying the code
run_regression("all");                            % every stage, one after the other
```

The report of each stage is saved in `<work_root>/<stage>/report.txt` (by default `work_root` is a folder in `tempdir`; it can be changed in `regression_config.m` or with the environment variable `CELLCIRCUIT_REGRESSION_DIR`).

## Notes

- Each stage gets the same input for both versions: the output of the previous stage, generated once with the **old** version (folder `<work_root>/chain`).
- `QualityControl_HyperSynch` and `QualityControl_Manual` run on the first `hs_n_clusters` clusters (default: 100) with `n_iteration_hs` shuffling iterations (default: 5); with all the clusters of the reference recording (~800) the original code would take weeks.
- `IntegrateDataSets` is tested with the Kilosort dataset only (the reference recording has no 3Brain dataset).
- `Analyze_spikeTrain` uses the `_curated_QCmanual.mat` file of the reference recording (generated by `IntegrateDataSets`), and only analyzes the conditions in `analyze_conditions` (default: the baseline), because the network burst detection takes most of the time (~1 min per condition with the revised code; ~30 min before).
- Timing differences below ~15% are within the run-to-run noise (especially with the profiler on, or when other programs are running).
- The pipeline functions accept an optional structure `inputs` that overwrites the values of their INPUTS / PARAMETERS sections (see `override_inputs.m`). Without it, they behave exactly as before.
