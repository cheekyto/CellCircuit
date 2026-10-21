# CellCircuit

**A MATLAB pipeline that turns multielectrode-array (MEA) recordings of brain organoids into quantitative measures of neuronal firing, bursting and network bursting, from spike-sorted data to group-level Excel tables, in one command.**

![MATLAB](https://img.shields.io/badge/MATLAB-R2025b-orange) ![License](https://img.shields.io/badge/license-GPL--3.0-blue)

Developed by **Sung Min Yang, Ph.D.**, Harvard University

---

## What it does

High-density MEA chips (e.g., 3Brain, 64 × 64 = 4,096 electrodes) record the electrical activity of hundreds of neurons at once. CellCircuit takes the spikes of every unit, either spike-sorted with Kilosort + BombCell or exported directly from 3Brain BrainWave, and:

1. **Cleans the units**: refractory-period violations, silent and noisy units, response to TTX, and electrical artifacts that appear as hypersynchronous spikes between pairs of units (a shuffle test, 500 iterations).
2. **Detects bursts** in every unit (maximal-interval method) and **network bursts** across the population: population events detected on an optimized bin size, the units that participate in them (shuffle-based statistics), primary and secondary events, and families of events found by clustering.
3. **Measures** firing (rate, variability, entropy, mutual information, spike-time tiling), bursting (rate, duration, spikes per burst, intervals) and network bursting (rate, participation, duration, amplitude, temporal structure of the onset and the peak), for every recording condition (e.g., baseline and each drug).
4. **Compiles** every recording of an experimental group into Excel tables, one row per recording condition, ready for statistics.

```mermaid
flowchart LR
    A["Kilosort + BombCell<br/>or 3Brain BrainWave"] --> B[ShapeSpikeTrainData]
    B --> C[QualityControl_General]
    C --> D[QualityControl_HyperSynch]
    D --> E[QualityControl_Manual]
    E --> F["SelectElectrodes<br/>(optional)"]
    F --> G[IntegrateDataSets]
    G --> H[Analyze_spikeTrain]
    H --> I[Compile_MEAproject]
```

| Stage | What it does | Output |
|---|---|---|
| `ShapeSpikeTrainData` | Groups the spikes of each unit and splits the recording into its conditions | `<name>_raw.mat` |
| `QualityControl_General` | Discards units with refractory violations, no activity, noise or no response to TTX | `<name>_cleaned.mat` |
| `QualityControl_HyperSynch` | Finds pairs of units with artifactual synchrony (shuffle test) | `<name>_QC.mat` |
| `QualityControl_Manual` | Final curation (a list of units to discard, or interactive plots) | `<name>_QCmanual.mat` |
| `SelectElectrodes` | Keeps the electrodes under the organoid (optional) | `<name>_<organoid>_QCmanual.mat` |
| `IntegrateDataSets` | Combines Kilosort and 3Brain units, removing duplicates | `<name>_curated_QCmanual.mat` |
| `Analyze_spikeTrain` | Spiking, bursting and network-burst analysis of every condition | `<name>_curated_Analyzed.mat` |
| `Compile_MEAproject` | Group tables: firing, bursting, network bursting | `<group>_Compiled.xlsx` |

---

## Highlights of the 2026 revision

In October 2026 the whole pipeline was revised, function by function, in nine reviewed phases plus a set of upgrades. Every change was tested on a reference recording.

**Faster**, with identical results (reference recording, ~800 units before quality control):

| Stage | Before | After |
|---|---|---|
| `QualityControl_HyperSynch` | ~850 h (estimated) | ~30 min |
| `Analyze_spikeTrain`, per recording condition | ~50 min | ~1 min |

The speed-ups came from profiling: vectorized spike counting, a cache of the shuffled spike trains, and removing repeated work.

**Correct.** The review found and fixed errors that changed results. Some examples:
- network bursts: the same population event could be counted both as a primary and as a secondary event;
- each secondary event was compared with the wrong parent event, and the spike times of deleted events were kept, so the features and the clustering read the wrong spikes;
- a quality-control metric of the network bursts was always empty, and one of the three detection methods always crashed;
- spikes that fell exactly on the boundary between two recording conditions were lost;
- the TTX conditions were never excluded from the activity threshold of the manual curation;
- the clustering of network bursts depended on the random state of the MATLAB session; it is now reproducible (fixed seeds).

Every change that alters a result was reviewed on its own, with its measured effect on the reference recording, before it was applied.

**Tested.** A regression harness (`tests/`) runs each stage with the old and the new code in fresh MATLAB processes, on copies of a reference recording, and compares every saved variable (relative tolerance 1e-10). It also reports the running time and a profiler summary. Unit tests cover the core functions.

**Reproducible.**
- One command runs the whole pipeline for a list of recordings: `run_pipeline("recordings.csv")`.
- All 232 parameters of the pipeline live in a single file, `pipeline_parameters.m`.
- The parameter set used for each recording is saved next to its data, with the history of its runs and the git version of the code.
- A recording can be rerun with the same set, with a different set, or with a few values changed.

**New analyses.**
- Variability of firing and bursting: CV of the intervals, the interval Fano factor (CV² and variance/mean), and the count Fano factor in 60-s windows with the mean count per window (Rajdl, Lansky & Kostal, 2020).
- A quality score for the families of network bursts.
- A stricter definition of secondary network bursts: they are kept only when more than 75% of the large (level-1) events are followed by one; otherwise those events are counted as primary events.

---

## Quick start

**Requirements.**
- MATLAB (developed and tested with R2025b), with the Signal Processing Toolbox and the Statistics and Machine Learning Toolbox.
- [npy-matlab](https://github.com/kwikteam/npy-matlab), to read the Kilosort output.

**1. List the recordings** in a table (see [`recordings_template.csv`](recordings_template.csv)): one row per recording, with its folder, the common name of its files and its group. Optional columns hold the units to discard, the cropping of the electrodes, a 3Brain dataset to integrate, and whether to include the recording.

**2. Run the pipeline:**

```matlab
addpath("path/to/CellCircuit")
run_pipeline("recordings.csv", 'extra_path', "path/to/npy-matlab")
```

Only the progress is printed; the messages of each stage go to a log file next to the recording. This is the output of a test run on the reference recording (shortened: 5 shuffling iterations instead of 500, and one recording condition; recording and group renamed):

```
=== run_pipeline (2026-10-20 20:16): 1 recording(s); stages: ShapeSpikeTrainData > QualityControl_General > ... > Compile_MEAproject
    parameters: pipeline_parameters.m
    the messages of the stages are saved in the log file of each recording (<file_name>_pipeline_log.txt)

[1/1] recording_01 (group: groupA) - 20:16
      ShapeSpikeTrainData          done (1 s)
      QualityControl_General       done (2 s)
      QualityControl_HyperSynch    done (1.7 min)
      QualityControl_Manual        done (3 s)
      SelectElectrodes             skipped (no cropping for this recording)
      IntegrateDataSets            done (0 s)
      Analyze_spikeTrain           done (3.1 min)

[group groupA] 1 recording(s) in .../groupA - 20:21
      Compile_MEAproject           done (2 s)

=== run_pipeline finished in 5.0 min: 1 recording(s), 0 error(s)
```

**3. Rerun or change things:**

```matlab
run_pipeline("recordings.csv", 'stages', "Analyze_spikeTrain")             % one stage, with its messages on screen
run_pipeline("recordings.csv", 'from', "IntegrateDataSets")                % from a stage to the end
run_pipeline("recordings.csv", 'parameters', "saved")                      % the same parameters as the last run
run_pipeline("recordings.csv", 'parameters', "my_parameters.m")            % a different parameter file
run_pipeline("recordings.csv", 'changes', {'NetworkBurstDetection.method3.THRESHOLD', 5.0})   % one value changed
```

Each stage can also be run by hand, as before (open it, edit its INPUTS section, press Run). In that case it uses the values written in its code.

---

## Outputs

- **For each recording:** the files of each stage (see the table above), `<name>_parameters.mat` (the parameter set and the history of the runs) and `<name>_pipeline_log.txt`.
- **For each group:** `<group>_Compiled.xlsx`, with one sheet per type of activity:
  - *Firing*: rate, CV and Fano factors of the inter-spike intervals, spike-count Fano factor, entropy, local variation, mutual information, spike-time tiling coefficient;
  - *Bursting*: rate, duration, spikes per burst, intervals within and between bursts with their variability, and tonic (non-bursting) units;
  - *Network bursting*: events, participating units, rate and its variability, duration, tail, amplitude, onset and peak structure, families of events and their clustering quality, and secondary events.

  `<group>_Grouped&Compiled.xlsx` has the same tables, with the rows grouped by treatment. The first row of every sheet gives the units.

---

## Repository layout

```
CellCircuit/
├── run_pipeline.m            master function (recordings table → group tables)
├── pipeline_parameters.m     every parameter of the pipeline, in one file
├── recordings_template.csv   example of a recordings table
├── <stage>.m                 the eight stages (ShapeSpikeTrainData ... Compile_MEAproject)
├── <helper>.m                detection, features, statistics and clustering functions
├── tests/                    regression harness and unit tests (see tests/README.md)
├── docs/FUNCTION_REFERENCE.md   detailed documentation of every stage and function, and the list of changes
└── legacy/                   functions that no stage uses any more
```

---

## Development

The revision was carried out in phases (0–8, plus a set of upgrades), each reviewed and committed separately in the lab's development repository; [`docs/FUNCTION_REFERENCE.md`](docs/FUNCTION_REFERENCE.md) lists every change. This public repository starts from the revised version (v2.0).
- Changes that keep the results identical were verified with the regression harness.
- Changes that alter results (bug fixes, new methods) were applied one at a time, each with its measured effect.

The work was done with AI-assisted pair programming ([Claude Code](https://claude.com/claude-code)). The author made the scientific decisions and approved every change that alters a result.

## License

[GNU General Public License v3.0](LICENSE)
