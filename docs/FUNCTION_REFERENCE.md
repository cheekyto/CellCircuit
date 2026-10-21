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
- **Clustering of the network bursts:** `Analyze_spikeTrain` prints `quality_clustering`, and the network-burst tables have the column `Clustering_Quality` (2 = reliable; 1 = it seems to work; 0 = not significant, single family; NaN = not run).
- **Secondary events:** the minimum number of participating units is 6 during the detection (`Optimal_secondaryPeak`, it was 5), as in the final quality control. In the last round of the detection (`NetworkBurstDetection`), after the quality control of each secondary event, the secondary events are accepted only if more than `GAP_WITH_PEAK` = 75% of the gaps after the level-1 primary events have at least `MIN_PEAK2_GAP` = 1 secondary event; otherwise they are all rejected, and the level-0 primary events that coincided with them are kept as primary events (the rule was in `NetworkBurstPlayers`, where it never ran). On the reference recording 98% of the gaps have secondary events: the results do not change.
- **Master function and single parameter file** (`run_pipeline`, `pipeline_parameters`, `load_parameters`, `apply_parameters`, `parameter_set`): the recordings are listed in a table (`recordings_template.csv`), and `run_pipeline` runs the pipeline end to end (Kilosort) or stage by stage (also the 3Brain stages). Every parameter of every function (232) is set in `pipeline_parameters.m`; while `run_pipeline` runs, each function replaces the values written in its code by those of the file (a stage run by hand still uses the values in its code). The set used is saved for each recording (`<file_name>_parameters.mat`, with the history of the runs and the version of the code), and a run can reuse it (`'parameters', "saved"`), use a different file, or change a few values (`'changes'`). `QualityControl_Manual` runs without windows (the clusters in the column `to_discard` are discarded). Only the progress is printed; the messages of the stages go to `<file_name>_pipeline_log.txt`. With the parameter set active, every stage gives exactly the same results as before on the reference recording.
- **Folder `legacy`:** the functions that no stage of the pipeline uses were moved there: `Fix_Analysis_Spiketrain`, `PopulationSilence`, `StructureNetworkBurst`, `bootstrap_optimal_k_new`, `build_distance`, `cluster_neuronal_events_new`, `exploreHDF5`, `hartigan_dip_test`, `makeVideo_MEA`, `organized_Finalresults`, `permutation_test_clustering_new`, `plotMEArecording` (and, before, `Brain3_processedData`). They are still on the MATLAB path when it is set with `addpath(genpath(...))`.



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

OUTPUT: SAVE the following variables in the file sprintf('%s\\%s_%s_QCmanual.mat',folder,file_data, organoid_label):
    - the same variables than the ones generated by "QualityControl\_Manual"
    - spikeTimestamp_BeforeCropping:    structure array 'spikeTimestamp' but before the cropping (it is the original output generated by "QualityControl\_Manual")
    - option:                           option used for the cropping
    - file_excel:                       excel file used for the cropping
    - ellipse_cropped:                  ellipse used for the cropping [x0 y0 radiusx radiusy] 

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

It calls the functions "SpikeFeatures", "BurstDetection", "BurstFeatures", "NetworkBurstDetection", "NetworkBurstFeatures" and "NetworkBurstStructure"



INPUT:	- folders / files\_data:	folders and names of the files to analyze (e.g., "XXX\_curated")

 	- condition\_analysis:	list of recording conditions to analyze (=-1, all of them)

 	- window:		analysis time at the end of each recording condition, in sec (=-1, the entire condition)

 	- opt\_rerun:		=1, the results are appended to an existing "XXX\_Analyzed.mat" file (e.g., to analyze other conditions); =0, new file

 	- NETWORK\_BURST\_DETECT\_METHOD:	=1, population firing rate with a fixed bin size; =2, normalized population firing rate with dynamic bin size; =3, method 2 + clustering of the events in families (default)

 	- NETWORK\_BURST\_ANALYSIS\_METHOD:	=1, single family of network bursts; =2, one set of features per family (default)



OUTPUT: SAVE the following variables in the file "XXX\_Analyzed.mat"

 	. spikeTimestamp, n\_clusterOriginal, n\_conditions, n\_cluster\_Kilosort, n\_cluster\_3Brain

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

 	. NETbursting0: network bursts, all the primary events together;   NETbursting1: only the largest family of primary events;   NETbursting2: one row per family of primary events (all of them include the secondary events)

 	It also saves the tables in an Excel file with the same name (the first row of each sheet has the units), and the file "group\_Grouped&Compiled.xlsx" with the rows grouped by treatment





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

