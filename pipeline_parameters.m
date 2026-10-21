% PARAMETERS OF THE WHOLE PIPELINE, used by "run_pipeline" (revised October 2026)
%
% One field per function: P.<function>.<parameter>. The values are the ones written in the code of each function (INPUTS / PARAMETERS
% sections, October 2026), except the few marked "(run_pipeline: ...)". While "run_pipeline" runs, every function uses the values of
% this file (see "apply_parameters"); a stage run by hand (e.g., pressing Run in "Analyze_spikeTrain") still uses the values in its code.
%
%   - A different set: copy this file (e.g., "parameters_threshold5.m"), edit the values, and run
%         run_pipeline("recordings.csv", 'parameters', "parameters_threshold5.m")
%   - A few values changed, only for one run (on top of this file, or of the saved set):
%         run_pipeline("recordings.csv", 'changes', {'NetworkBurstDetection.method3.THRESHOLD', 5.0})
%   - A parameter removed from this file keeps the value written in its function.
%   - "settings": the options of the INPUTS sections that are not paths (folders, file names and groups come from the recordings table).
%   - P.links (at the end): parameters that must have the same value in two functions; only the first one is written.
%   - The set used for each recording is saved in "<folder>/<file_name>_parameters.mat", with the history of the runs.
%   - Not included: the unit constants (K_MSEC_SEC, K_SEC_MIN), and "Compare_Kilosort_vs_3Brain" (not run by "run_pipeline").

P = struct();

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  1) ShapeSpikeTrainData (Kilosort)

% ShapeSpikeTrainData
P.ShapeSpikeTrainData.unitType_good = [1 2 3];                  % classification of the unit type, performed by "bombcell", that will be considered for downstream analysis
                                                                % =0, noise; =1, good; =2, MUA; =3, non-somatic


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  2) QualityControl_General

% QualityControl_General: settings (INPUTS section of the function)
P.QualityControl_General.condition = -1;                        % =-1, the entire recording (including all conditions)
                                                                % otherwise, the list of conditions (it must be in order)
P.QualityControl_General.option_plot = 0;                       % = 1, plot the raster plot and some metrics
                                                                % (run_pipeline: no plots; in the code, 1)
                                                                % = 0, do not plot (just save the output)
P.QualityControl_General.option_pause = 0;                      % = 1, wait for a key press before moving on to the next file
                                                                % (run_pipeline: no pause; in the code, 1)
                                                                % = 0, do not wait (e.g., batch processing)
P.QualityControl_General.condition_3Brain = 0;                  % = 0, the data does not come directly from 3Brain (e.g., spike sorting done by Kilosort)
                                                                % (run_pipeline: Kilosort data; in the code, 1)
                                                                % = 1, the data comes directly from the 3Brain platform (no intermediate spike sorting)

% QualityControl_General
P.QualityControl_General.REFRACTORY_PERIOD = 1.0;               % absolute refractory period for neuronal firing, in msec (default: 1.0)
P.QualityControl_General.REF_PERIOD_THRESHOLD = 5;              % threshold for the percentage of spikes that violate the refractory period criteria; above this threshold, the cluster will be discarded (default: 5%)
P.QualityControl_General.REF_PERIOD_MINSPIKES = 3;              % threshold for the absolute number of spikes that violate the refractory period criteria (default: 3); to account for the low activity units
P.QualityControl_General.ACTIVITY_THRESHOLD = 0.03;             % minimum level of spiking frequency for an individual electrode to be considered active, in Hz (0.02 Hz = 1/50 sec; 0.03 Hz = 1/33 sec)
P.QualityControl_General.CV_THRESHOLD = 3.5;                    % threshold for coefficient of variation to separate phasic vs tonic activity (default: 3.5)
P.QualityControl_General.LVR_THRESHOLD = 1.00;                  % threshold for local variation ratio to separate phasic vs tonic activity (default: 1.00)
P.QualityControl_General.PERCENTILE = 98;                       % percentile to identify the outliers with too high firing rate, among ONLY the clusters with tonic activity (default: 98)
P.QualityControl_General.PERCENTILE_TOP = 99.5;                 % percentile to identify the outliers with too high firing rate, among ALL clusters (default: 99.5)
P.QualityControl_General.ISI_RATIO = 8;                         % threshold for the ratio of isi during bursting vs non-bursting; the lower, the more permissive (default: 8)
P.QualityControl_General.N_SPIKES_BURST = 0.90;                 % threshold for the proportion of spikes in bursts; if the cluster is mostly bursting, probably is a neuron (default: 0.90)
P.QualityControl_General.TTX_EFFECT = 0.15;                     % minimum effect of TTX to be considered; proportion relative to the maximum frequency (0.10 - 0.20)
P.QualityControl_General.BURST_METHOD = 1;                      % method for the burst detection (1)
P.QualityControl_General.MIN_LENGTH = 90;                       % minimum duration of a recording condition, in sec (default: 90)
P.QualityControl_General.N_PLOT = 20;                           % number of clusters for raster plots (20)
P.QualityControl_General.SHIFT = 0.1;                           % distance between rows in the raster plots (0.1)
P.QualityControl_General.LABEL_POS = 0.6;                       % position for a message in the raster plot (0.6)
P.QualityControl_General.MATRIX_SIZE_X = 64;                    % horizontal length of the electrode array (64*64=4096)
P.QualityControl_General.MATRIX_SIZE_Y = 64;                    % vertical length of the electrode array (64*64=4096)

% lvr
P.lvr.MIN_N = 3;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  3) QualityControl_HyperSynch

% QualityControl_HyperSynch
P.QualityControl_HyperSynch.N_ITERATION = 500;                  % number of iterations for the shuffling (default: 500-1000)
P.QualityControl_HyperSynch.N_SPIKES_MINIMUM = 5;               % minimum number of hypersynchronized spikes to trigger assessment; it will apply to any category (bursting, nonbursting and total) (default: 5)
P.QualityControl_HyperSynch.N_SPIKES_BURST = 8;                 % threshold for the percentage of hypersynchonized spikes within burst, in % (default: 8%)
P.QualityControl_HyperSynch.N_SPIKES_NONBURST = 8;              % threshold for the percentage of hypersynchonized spikes outside burst, in % (default: 8%)
P.QualityControl_HyperSynch.N_SPIKES_TOTAL = 8;                 % threshold for the percentage of hypersynchonized spikes in all categories, in % (default: 8%)
P.QualityControl_HyperSynch.FACTOR_MINORITY = 2;                % factor to scale up the threshold when analyzing the minority category (default: 2)
P.QualityControl_HyperSynch.PERCENTILE_BURST = 95;              % threshold percentile based on proportion of hypersynchonized spikes within burst, in % (default: 95%)
P.QualityControl_HyperSynch.PERCENTILE_NONBURST = 95;           % threshold percentile based on proportion of hypersynchonized spikes outside burst, in % (default: 95%)
P.QualityControl_HyperSynch.PERCENTILE_TOTAL = 95;              % threshold percentile based on proportion of hypersynchonized spikes in all categories, in % (default: 95%)
P.QualityControl_HyperSynch.N_SPIKES_BURST_CUT = 4;             % threshold for the percentage of hypersynchonized spikes within burst to make the first cut and proceed with the shuffling, in % (default: 4%)
P.QualityControl_HyperSynch.N_SPIKES_NONBURST_CUT = 4;          % threshold for the percentage of hypersynchonized spikes outside burst to make the first cut and proceed with the shuffling, in % (default: 4%)
P.QualityControl_HyperSynch.N_SPIKES_TOTAL_CUT = 4;             % threshold for the percentage of hypersynchonized spikes to make the first cut and proceed with the shuffling, in % (default: 4%)
P.QualityControl_HyperSynch.METHOD_SHUFFLING = 2;               % method for the spike-train shuffling (=1, uniform; =2, biphasic <burst and non-bust>)
P.QualityControl_HyperSynch.OPTION_METHOD_BURST = 1;            % method for the bursting detection (=1, Maximal Interval)
P.QualityControl_HyperSynch.REFRACTORY_PERIOD = 0;              % refractory period for the shuffling, in msec (default: 0)
P.QualityControl_HyperSynch.option_save_rawQC = 0;              % =1, also save "rawQC" and "scores" (as before; it can take several GB for ~800 clusters)
                                                                % =0, only save the compact results: "hs_master", "hs_final", "hs_percentage" and "red_flag" (revised October 2026)

% HyperSyncPair
P.HyperSyncPair.JITTER = 3;                                     % threshold to define the hypersynchronicity, relative to t_step (default: 3.0)
P.HyperSyncPair.MAX_JITTER = 1.0;                               % maximum threshold to define the hypersynchronicity, in msec (default: 1.0)

% ShuffleSpikeTrain
P.ShuffleSpikeTrain.BIN = 50;                                   % maximum bin size, for option #3, in ms (default: 50)
P.ShuffleSpikeTrain.BURST_METHOD = 1;                           % burst detection method (=1      Maximal Interval)
P.ShuffleSpikeTrain.DEFAULT_SEED = 0;                           % seed for the random number generator when the input "seed" is not provided


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  4) QualityControl_Manual (run_pipeline: no windows, the clusters listed in the column "to_discard" of the recordings table are discarded)

% QualityControl_Manual: settings (INPUTS section of the function)
P.QualityControl_Manual.option_main = 3;                        % =0, show the clusters that do not meet the threshold for activity level ("minimum_activity")
                                                                % =1, plot raster plots for all units
                                                                % =2, plot firing frequency vs CV of ISI, firing frequency vs LVR, and raster plot for outstanding clusters
                                                                % =3, selet clusters to discard (before waveform_based classification) and plot them againts good clusters
                                                                % =4, run waveform-based classification (only for spike-sorted by 3Brain) and plot the good clusters (WARNING: make sure you run it after completing the other steps)
                                                                % =5, plot the clusters discarded based on hypersynchronicity
P.QualityControl_Manual.option_saveComplete = 1;                % option to save the results and end the manual curation
P.QualityControl_Manual.minimum_activity = 0.04;                % threshold for the activity level, in Hz (default: 0.04 Hz -> 1/25sec)
P.QualityControl_Manual.cond_low_activity = -1;                 % indexes of the conditions to analyze the minimum activity level (=-1, all experimental conditions but TTX)
                                                                % the unit is silent if the activity level in every condition is below the threshold

% QualityControl_Manual
P.QualityControl_Manual.METHOD_BURST = 1;                       % method for the burst detection (default: 1)
P.QualityControl_Manual.N_PLOT = 30;                            % number of clusters to plot in each raster plot (default: 30)
P.QualityControl_Manual.N_OUTSTANDING = 20;                     % number of outstanding clusters (when looking at the CV, LVR and frequency) (default: 20)
P.QualityControl_Manual.SHIFT = 0.1;                            % distance between rows in a raster plot (default: 0.1)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  5) SelectElectrodes (only for the recordings with "crop_ellipse" or "crop_excel" in the recordings table)

% SelectElectrodes: settings (INPUTS section of the function)
P.SelectElectrodes.CROPPING_SAFETY = 1.10;                      % proportion of the actual ellipse to be cropped (default: 1.10)

% SelectElectrodes
P.SelectElectrodes.SHEET = 4;                                   % sheet within the excel file to read the coordinates of the selected electrodes
P.SelectElectrodes.COLUMN = 4;                                  % column within the excel file table that contains the channel id
P.SelectElectrodes.ROW = 5;                                     % row witin the excel file table that contains the first data point
P.SelectElectrodes.MATRIX_SIZE = 64;                            % dimension of the electrode array (64*64=4096)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  6) IntegrateDataSets

% IntegrateDataSets: settings (INPUTS section of the function)
P.IntegrateDataSets.unitSelect_Kilosort = [1; 2; 3];            % classification of the unit type, performed by "bombcell", that will be considered for downstream analysis
                                                                % =0, noise; =1, good; =2, MUA; =3, non-somatic
P.IntegrateDataSets.unitSelect_3Brain = [1];                    % waveform-based classification of the unit type, that will be considered for downstream analysis
                                                                % =0, unknown; =1, good; =2, weak; =-1, noise
P.IntegrateDataSets.option_criteria = 3;                        % =1, only spatial superposition
                                                                % =2, only spike train similarity
                                                                % =3, spatial superposition or spike train similarity
                                                                % =4, spatial superposition and spike train similarity
P.IntegrateDataSets.option_plot = 0;
P.IntegrateDataSets.option_print_list = 0;

% IntegrateDataSets
P.IntegrateDataSets.N_SPIKES_MINIMUM = 20;                      % minimum number of hypersynchronized spikes to trigger assessment; it will apply to any category (bursting, nonbursting and total) (default: 20)
P.IntegrateDataSets.N_SPIKES_BURST_CUT = 50;                    % threshold for the percentage of hypersynchonized spikes within burst, in % (default: 50%)
P.IntegrateDataSets.N_SPIKES_NONBURST_CUT = 20;                 % threshold for the percentage of hypersynchonized spikes outside burst to make the first cut and proceed with the shuffling, in % (default: 20%)
P.IntegrateDataSets.N_SPIKES_TOTAL_CUT = 40;                    % threshold for the percentage of hypersynchonized spikes to make the first cut and proceed with the shuffling, in % (default: 40%)
P.IntegrateDataSets.N_SPIKES_3BRAIN_CUT = 50;                   % threshold for the percentage of hypersynchonized spikes in the 3Brain signal, in % (default: 60%)
P.IntegrateDataSets.OPTION_METHOD_BURST = 1;                    % method for the bursting detection (=1, Maximal Interval)
P.IntegrateDataSets.REFRACTORY_PERIOD = 0;                      % refractory period for the shuffling, in msec (default: 0)
P.IntegrateDataSets.METHOD_SHUFFLING = 2;                       % method for the spike-train shuffling (=1, uniform; =2, biphasic <burst and non-bust>)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  7) Analyze_spikeTrain

% Analyze_spikeTrain: settings (INPUTS section of the function)
P.Analyze_spikeTrain.condition_analysis = -1;                   % list of recording condition to analyze
                                                                % (run_pipeline: all the recording conditions; in the code, 1)
                                                                % =-1, all the recording conditions
P.Analyze_spikeTrain.window = 600;                              % the analysis time, at the end of each recording condition, in sec (default: 300 sec)
                                                                % =-1, the entire recording time, for each recording condition
P.Analyze_spikeTrain.opt_rerun = 0;                             % =1, it is a rerun, an output file already exists and the new data will be appended
                                                                % (run_pipeline: a new "_Analyzed.mat" file; in the code, 1)
                                                                % =0, it is the first time the analysis of this experiment is run
P.Analyze_spikeTrain.opt_distance = 0;                          % option to remove primary network burst events that are too close to the previous one
P.Analyze_spikeTrain.opt_clusterAgain = 0;                      % option to re-run the clustering of network bursting events
P.Analyze_spikeTrain.BURST_DETECT_METHOD = 1;                   % algorithm to be used for individual burst detection
                                                                %  =1       Maximal Interval
P.Analyze_spikeTrain.NETWORK_BURST_DETECT_METHOD = 3;           % algorithm to be used for newtork burst detection
                                                                %  =1       Population firing rate
                                                                %  =2       looking for peaks, combining the regular population firing rate and the normalized population firing; and using a progressive narrowing of the binSize
                                                                %  =3       option 2 + clustering
P.Analyze_spikeTrain.NETWORK_BURST_ANALYSIS_METHOD = 2;         % algorithm to be used for network burst analysis
                                                                %  =1       single family of network bursting, with a single set of participating channels
                                                                %  =2       cluster analysis to identify different families of network bursts
P.Analyze_spikeTrain.opt_checking = 0;
P.Analyze_spikeTrain.print_Plots = 1;                           % option to save the rater plot of all neurons
P.Analyze_spikeTrain.plot_option = 0;
P.Analyze_spikeTrain.plot_network_burst = 1;

% Analyze_spikeTrain
P.Analyze_spikeTrain.MIN_N_BURST = 5;                           % minimum number of burst events to be considered a bursting neuron (default: 5)
P.Analyze_spikeTrain.MIN_WINDOW = 0.4;
P.Analyze_spikeTrain.TREATMENT = {'CGP', 'Gabazine', 'DNQX', 'D-AP5', 'TTX'};          % labels of the treatments, matched to the filename of each condition
P.Analyze_spikeTrain.RATIO_MATCH = 0.60;                        % minimum matching between the treatment labels and the filename


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  7a) Spiking and bursting features

% SpikeFeatures
P.SpikeFeatures.N_BIN_ENTROPY = 30;                             % temporal bin to calculate the entropy and the mutual information, in msec (default: 30)
P.SpikeFeatures.DT_STTC = 15;                                   % dt for the spike train tiling coefficient, in msec (default: 15)
P.SpikeFeatures.COUNT_WINDOW = 60;                              % length of the windows to count the spikes for the Fano factor of the counts, in sec (default: 60)

% BurstDetection
P.BurstDetection.RATIO_PEAK = 0.5;                              % proportion of spikes in the burst that are grouped to look for the peak of the burst (default: 0.5)

% BurstDetection: parameters of the option method1 (they depend on the method)
P.BurstDetection.method1.MIN_SPIKE_BURST = 3;                   % minimum number of spikes within a burst (default: 3)
P.BurstDetection.method1.ISI_START_BURST = 200;                 % threshold (maximum) ISI to identify the beginning of a burst, in ms (default: 200 -> 5Hz)
P.BurstDetection.method1.ISI_PEAK_BURST = 150;                  % minimum peak ISI to accept a burst, in ms (default: 150 -> 6.7Hz [100ms -> 10Hz])
P.BurstDetection.method1.ISI_END_BURST_INI = 300;               % initial threshold (maximum) ISI to identify the end of a burst, in ms (default: 300ms -> 3.3Hz)
P.BurstDetection.method1.FACTOR_MATURATION = 1.05;              % factor to increase the ISI_END_BURST as the burst has more units (default: 1.05 -> 5% increase per unit, starting from MIN_SPIKE_BURST)
P.BurstDetection.method1.ISI_END_BURST_MAX = 500;               % maximum threshold (maximum) ISI to identify the end of a burst, in ms (500 -> 2Hz)
P.BurstDetection.method1.MIN_IBI = 500;                         % minimum inter-burst interval to identify two consecutive bursts, in ms (default: 500); from finish to start of the next one
P.BurstDetection.method1.MIN_BURST_DURATION = 25;               % minimum burst duration, in ms (default: 25 -> with 3 spikes, it means 80 Hz)

% BurstFeatures
P.BurstFeatures.COUNT_WINDOW = 60;                              % length of the windows to count the bursts for the Fano factor of the counts, in sec (default: 60)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  7b) Network burst detection

% NetworkBurstDetection
P.NetworkBurstDetection.MIN_CHAN = 6;                           % minimum number of participating channels/units to classify an organoid as network bursting (default: 6)
P.NetworkBurstDetection.MIN_CHAN_BURST = 5;                     % minimum number of channels/units in a single network burst (default: 5)
P.NetworkBurstDetection.THRESHOLD_BURST_WITH_UNITS = 0.50;      % minimum proportion of potential network burst with MIN_CHAN_BURST units (default: 0.50)
P.NetworkBurstDetection.MIN_CHAN_RATIO = 0.075;                 % minimum number of participating channels/units to classify an organoid as network bursting, as a proportion of the total number of units (default: 0.075)
P.NetworkBurstDetection.MIN_CHAN_BURST_RATIO = 0.01;            % minimum number of channels/units in a single network burst, as a proportion of the total number of units (default: 0.05)
P.NetworkBurstDetection.FACTOR_FIRST_ROUND = 0.80;              % (default: 0.8)
P.NetworkBurstDetection.MIN_CHAN_SECONDARY = 0.25;              % minimum number of participating channels/units to accept the secondary events, in terms of the number of participating units in the primary events (default: 0.25)
P.NetworkBurstDetection.MIN_PEAK2_GAP = 1;                      % minimum number of secondary events per gap/window, after a level-1 primary event (default: 1)
P.NetworkBurstDetection.GAP_WITH_PEAK = 0.75;                   % minimum proportion of gaps/windows (level-1 primary events) followed by MIN_PEAK2_GAP secondary events, to accept
                                                                % the secondary events (default: 0.75)

% NetworkBurstDetection: parameters of the option method1 (they depend on the method)
P.NetworkBurstDetection.method1.BIN = 0.3;                      % size of the temporal bin to calculate the population firing rate, in sec (0.2-0.4)
P.NetworkBurstDetection.method1.THRESHOLD = 5.0;                % in terms of the standard deviation (ventral: 3.0 - 5.0 ; dorsal: 5.0 - 8.0)
P.NetworkBurstDetection.method1.MIN_ROBUST = 4;                 % minimum number of robustly participating clusters (default: 4, as in method 2)

% NetworkBurstDetection: parameters of the option method2 (they depend on the method)
P.NetworkBurstDetection.method2.BIN_MIN = 0.01;                 % minimum size of the temporal bin to calculate the population firing rate, in sec (default: 0.01 sec)
P.NetworkBurstDetection.method2.BIN_MAX = 3.0;                  % maximum size of the temporal bin to calculate the population firing rate, in sec (default: 3.0 sec)
P.NetworkBurstDetection.method2.THRESHOLD = 5.0;                % in terms of the standard deviation (default: 5.0)
P.NetworkBurstDetection.method2.BIN_DEFAULT = 0.20;             % default size of the temporal bin to calculate the population firing rate, in sec (default: 0.25 sec)
P.NetworkBurstDetection.method2.MIN_ROBUST = 4;                 % minimum number of robustly participating clusters (default: 4)

% NetworkBurstDetection: parameters of the option method3 (they depend on the method)
P.NetworkBurstDetection.method3.BIN_MIN = 0.01;                 % minimum size of the temporal bin to calculate the population firing rate, in sec (default: 0.01 sec)
P.NetworkBurstDetection.method3.BIN_MAX = 3.0;                  % maximum size of the temporal bin to calculate the population firing rate, in sec (default: 3.0 sec)
P.NetworkBurstDetection.method3.THRESHOLD = 6.0;                % in terms of the standard deviation (default: 5.0; it must be >3.5)
P.NetworkBurstDetection.method3.BIN_DEFAULT = 0.10;             % default size of the temporal bin to calculate the population firing rate, in sec (default: 0.15 sec)
P.NetworkBurstDetection.method3.MIN_ROBUST = 5;                 % minimum number of robustly participating clusters (default: 5)

% DynamicPopActDetect
P.DynamicPopActDetect.SCALING = 2.0;                            % factor to scale up the bin size to run the dynamic search (default: 2.0)
P.DynamicPopActDetect.SCALING_ADJ = 0.75;                       % factor to progressively change the bin size (default: 0.75)
P.DynamicPopActDetect.N_SHORTLIST = 3;                          % number of pre-candidates to build the shortlist of optimal bin sizes
P.DynamicPopActDetect.CHANGE_THRESHOLD = 0.40;                  % threshold for the metric: change of number of burst events (default: 0.40)
P.DynamicPopActDetect.OPTION_LEVEL = 1;                         % =1, just analyze the level1 peaks to determine the optimal bin size and the secondary peaks
                                                                % =0, analyze all peaks (level0 and level1) to determine the optimal bin size
P.DynamicPopActDetect.MAX_IT = 100;

% PopActivityPeak_Integration
P.PopActivityPeak_Integration.BIN_DEFAULT1 = 0.100;             % default size of the temporal bin to calculate the population firing rate, in sec (default: 0.10 sec)
P.PopActivityPeak_Integration.BIN_DEFAULT2 = 0.150;             % default size of the temporal bin to calculate the population firing rate, in sec (default: 0.15 sec)

% PopulationActivityDetection
P.PopulationActivityDetection.START_LEVEL = 0.20;               % proportion of the peak amplitude of the network burst to identifity the starting point (default: 0.10)
P.PopulationActivityDetection.END_LEVEL = 0.20;                 % proportion of the peak amplitude of the network burst to identifity the ending point (default: 0.10)
P.PopulationActivityDetection.SMOOTH_FACTOR = [3; 5; 7];        % standard smoothing window length, in number of bins (default: 3; odd for symmetric smoothing)
P.PopulationActivityDetection.MIN_POINTS = 3;
P.PopulationActivityDetection.OUTLIER_TOP = 99.5;               % upper percentile to identify outliers, before calculating the noise (dafault: 99)
P.PopulationActivityDetection.PERCENTILE_TOP = 95.0;            % upper percentile to identify the "big mountains" (default: 95)
P.PopulationActivityDetection.CHAIN = 5;                        % total number of consecutive high-value to identify "big mountains" (default: 3)
P.PopulationActivityDetection.FACTOR_MOUNTAIN1 = 2.0;           % (default: 1.8) (in 100 bin and 5 high-values, the chances of having 3 or more consecutives high-values is 1/170)
P.PopulationActivityDetection.FACTOR_MOUNTAIN2 = 2.5;           % (default: 2.5) (in 100 bin and 5 high-values, the chances of having 3 or more consecutives high-values is 1/170)
P.PopulationActivityDetection.BIN_REDUCT = 0.30;                % (default: 0.30)
P.PopulationActivityDetection.MIN_THRESHOLD = 3.0;              % (default: 3.0)
P.PopulationActivityDetection.ACTIVITY_THRESHOLD = 0.03;        % minimum level of spiking frequency for an individual electrode to be considered active, in Hz (default: 0.05)
                                                                % (revised October 2026: 0.05 -> 0.03 Hz, user's choice)
P.PopulationActivityDetection.PEAK_THRESHOLD = 20;              % threshold to run the second round of searching, relative to median firing rate (default: 20)

% refinePeakDetection
P.refinePeakDetection.smoothWindow = 3;                         % smoothing window length, in number of bins (default: 3; odd for symmetric smoothing)
P.refinePeakDetection.WINDOW_FACTOR = 0.5;                      % in terms of the width of the peak (default: 0.5)
P.refinePeakDetection.BIN_FACTOR = [0.25; 0.40; 0.70; 1.0; 1.20]; % in terms of the original bin (default: 0.4, 0.7, 1.0)
P.refinePeakDetection.MIN_SAMPLES = 3;

% tailPopulationPeaks
P.tailPopulationPeaks.WINDOW_FACTOR = 2.0;                      % in terms of the width of the peak (default: 2.0)
P.tailPopulationPeaks.FACTOR_SECONDARY = 0.5;                   % factor to rescale WINDOW_FACTOR when analyzing secondary peaks (default: 0.5)
P.tailPopulationPeaks.MIN_SAMPLES = 3;

% searchActDrop
P.searchActDrop.MOViNG_STEP = 0.25;                             % steps to move the averaging window, in term of the size of the window
P.searchActDrop.THRESHOLD_BIS = 10.0;                           % (default: 10.0)

% evaluatePrimaryBursts
P.evaluatePrimaryBursts.N_BURST_MIN = 8;                        % (default: 8)
P.evaluatePrimaryBursts.MIN_DISTANCE = 4.0;                     % (default: 4.0)

% test_bimodality
P.test_bimodality.RNG_SEED = 0;                                 % seed for the random number generator (change it to get a different, but still reproducible, result)
P.test_bimodality.THRESHOLD_BC = 0.555;                         % (default: 0.555)

% PopulationActivitySize
P.PopulationActivitySize.MIN_CHAN = 6;                          % minimum number of channels involved in a single network burst (default: 5)
P.PopulationActivitySize.MIN_MULTIPLE = 5;                      % minimum number of channels with multiple spikes involved in a single network burst (default: 3)
P.PopulationActivitySize.MULTIPLE_SPIKE = 2;                    % minimum number of spikes to label as multiple (default: 2)
P.PopulationActivitySize.SAFE_K = 0.40;                         % when saving the spikes for each network burst event, the extra time is added on each side (default: 0.4)

% Optimal_secondaryPeak
P.Optimal_secondaryPeak.SMOOTH_ADD = 0;                         % increament of the standard smoothing window length, to look for secondary peaks, in number of bins (default: 1; even for symmetric smoothing)
P.Optimal_secondaryPeak.BIN_MIN = 0.010;                        % minimum size of the bin (default: 0.008 sec)
P.Optimal_secondaryPeak.WINDOW_FRACTION = 1/4;                  % to partition the gap between peaks and assess the general level of activity  right after primary peaks (default: 1/4)
P.Optimal_secondaryPeak.THRESHOLD_BROAD = 1.80;                 % threshold to assess the activity right after primary peaks; ratio to the last partition (default: 1.80)
P.Optimal_secondaryPeak.THRESHOLD_WINDOWS = 0.40;               % minimum number of primary peaks with a signifnicant tail activity (default: 0.50-0.40)
P.Optimal_secondaryPeak.PERCENTILE_2GAP = 0.20;                 % percentile of the distance between potential secondary peaks to determine the "time scale" (default: 0.20)
P.Optimal_secondaryPeak.SEARCH_VALLEY = 1.0;                    % proportion of the "typical_gap" to start the search of secondary peaks (default: 1.0)
P.Optimal_secondaryPeak.THRESHOLD_SECONDARY = 0.60;             % factor to scale the threshold for the primary peaks and use it for the secondary peaks (default: 0.6-0.5)
% MIN_CHAN: the same value as NetworkBurstDetection.MIN_CHAN (see P.links at the end)
% MIN_CHAN_SECONDARY: the same value as NetworkBurstDetection.MIN_CHAN_SECONDARY (see P.links at the end)
P.Optimal_secondaryPeak.BIN_DEFAULT = 0.100;                    % default size of the temporal bin to calculate the population firing rate, in sec (default: 0.100 sec)
P.Optimal_secondaryPeak.BIN2_FACTOR_LOW = [0.40; 0.20; 1.00];   % reduction factors tested for the bin size, when the bin of the primary events is < BIN_DEFAULT (opt_binopti = 1)
P.Optimal_secondaryPeak.BIN2_FACTOR_HIGH = [0.30; 0.12; 0.75; 0.05]; % reduction factors tested for the bin size, when the bin of the primary events is >= BIN_DEFAULT (opt_binopti = 1)
P.Optimal_secondaryPeak.BIN2_FACTOR_FIXED = 0.30;               % reduction factor for the bin size, without optimization (opt_binopti = 0) (default: 0.30)
                                                                % (revised October 2026: the three factors above were written in the "switch" below; same values)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  7c) Participating units (statistical assessment)

% NetworkBurstPlayers
P.NetworkBurstPlayers.SAFE_K = 0.40;                            % when saving the spikes for each network burst event, the extra time is added on each side (default: 0.4)
                                                                % (revised October 2026: in the statistical assessment, SAFE_K*duration before the event and 2*SAFE_K*duration after it, intended)
P.NetworkBurstPlayers.MIN_CHAN_BURST = 5;                       % minimum number of channels involved in a single network burst (default: 5)
P.NetworkBurstPlayers.MIN_CHAN_BURST_RATIO = 0.05;              % minimum proportion of bursting channels involved in a single network burst, relative to the median of bursting participating units (default: 0.4)
P.NetworkBurstPlayers.MIN_MULTIPLE = 3;                         % minimum number of channels with multiple spikes involved in a single network burst (default: 3)
P.NetworkBurstPlayers.MIN_MULTIPLE_RATIO = 0.03;                % minimum proportion of channels with multiple spikes involved in a single network burst, relative to the total number of units (default: 0.03)
P.NetworkBurstPlayers.MIN_PARTICIPATION_TONIC = 0.65;           % proportion of solid network burst in which a non-bursting channel has to participated to be counted (default: 0.95)
P.NetworkBurstPlayers.MIN_PARTICIPATION_PHASIC = 0.25;          % proportion of solid network burst in which a bursting channel has to participated to be counted (default: 0.6)
P.NetworkBurstPlayers.MIN_SOLID_BURST = 3;                      % minimum number of solid burst to classify the channels as solid vs weak (default: 3)
% MIN_PEAK2_GAP: the same value as NetworkBurstDetection.MIN_PEAK2_GAP (see P.links at the end; used only if FAST_MODE = 0)
% GAP_WITH_PEAK: the same value as NetworkBurstDetection.GAP_WITH_PEAK (see P.links at the end; used only if FAST_MODE = 0)
P.NetworkBurstPlayers.MIN_PARTICIPATION_SECONDARY = 0.40;       % minimum proportion of secondary network bursts in which a channel has to participate (default: 0.40)
P.NetworkBurstPlayers.N_PLOT = 20;
P.NetworkBurstPlayers.SHIFT = 1;
P.NetworkBurstPlayers.FAST_MODE = 1;                            % =1, don't search for participating units for secondary peaks until the last round
                                                                % =0, seach for participating units for secondary peaks

% findChannel_in_Peak
P.findChannel_in_Peak.SAFE_K = 0.40;                            % when saving the spikes for each network burst event, the extra time is added on each side (default: 0.4)
P.findChannel_in_Peak.ACTIVITY_THRESHOLD = 3.0;                 % threshold for the ratio between the firing rate within potential network bursts and overall firing rate, for individual neurons (default: 3.0)
P.findChannel_in_Peak.K_ONSET = 0.20;                           % when assessing the statistical significance of network burst participation, the extra time is added at the begining (default: 0.2)
P.findChannel_in_Peak.MIN_MULTIPLE = 2;                         % minimum number of spikes a unit has to display during a network burst to be considered a multiple unit (default: 2)
P.findChannel_in_Peak.N_SHUFFLING = 500;                        % number of shuffling repetitions for statistical assessment (default: 500)
P.findChannel_in_Peak.RNG_SEED = 0;                             % base seed for the shuffling; channel ch uses RNG_SEED + ch (change it to get a different, but still reproducible, set of shuffles)
P.findChannel_in_Peak.STAT_THRESHOLD = 0.95;                    % threshold of statistical significance (default: 0.95)
P.findChannel_in_Peak.MIN_PARTICIPATION = 0.20;                 % minimum proportion of network bursts in which a channel has to participate (default: 0.20)
P.findChannel_in_Peak.MIN_PARTICIPATION_ABS = 2.0;              % minimum number of network bursts in which a channel has to participate (default: 2)
P.findChannel_in_Peak.ACTIVITY_THRESHOLD_FINAL = 1.80;          % (used during statistical assessment) threshold for the ratio between the firing rate within potential network bursts and overall firing rate (default: 2.0, NOT higher)
                                                                % lower it (1.60) when there is high level of population activity (it is not a few individual units with high level of activity)

% familyParticipants
% K_ONSET: the same value as findChannel_in_Peak.K_ONSET (see P.links at the end)
% MIN_PARTICIPATION: the same value as findChannel_in_Peak.MIN_PARTICIPATION (see P.links at the end)
% MIN_PARTICIPATION_ABS: the same value as findChannel_in_Peak.MIN_PARTICIPATION_ABS (see P.links at the end)

% shuffleCache
P.shuffleCache.MAX_CACHE_MB = 1000;                             % maximum memory for the cache, in MB (~175 MB for the reference recording, 87,000 spikes)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  7d) Clustering of the network burst events

% clusterNetburst
P.clusterNetburst.MIN_PRIMARY_EVENTS = 3;                       % the clustering runs when there are MORE than MIN_PRIMARY_EVENTS primary events (i.e., at least 4)
                                                                % (revised October 2026: comment clarified; before, "minimum number of primary events to run the clustering")

% cluster_neuronal_events
P.cluster_neuronal_events.N_SCALAR = 7;
P.cluster_neuronal_events.MIN_SECONDARY_EVENTS = 5;             % minimum number of secondary events to include them
P.cluster_neuronal_events.MIN_SIZE_CLUSTER = 0;                 % minimum size of the clusters to run the boostrapping (default: 3)
P.cluster_neuronal_events.colors_event = {'r', 'b', 'c', 'g', 'm'};
P.cluster_neuronal_events.w_scalar = 0.35;                      % (default: 0.35)
P.cluster_neuronal_events.w_population = 0.65;                  % (default: 0.65)
P.cluster_neuronal_events.MAX_K = 4;                            % maximum number of clusters to evaluate (default: 4)
P.cluster_neuronal_events.LINKAGE_METHOD = 'ward';              % 'ward' | 'average' | 'complete' | 'single'
P.cluster_neuronal_events.ALPHA = 0.05;                         % significance threshold (default: 0.05)
P.cluster_neuronal_events.MIN_DIFF = 0.10;                      % maximum difference in "vote_frac" and "silhouette score" to consider to k-sizes non-different (default: 0.10)
P.cluster_neuronal_events.SIL_THRESHOLD = 0.40;                 % AND have a mean silhouette >= 0.25 (0.25 is a modest but non-trivial structure threshold) (the higher the value, the more conservative)

% bootstrap_optimal_k
P.bootstrap_optimal_k.RNG_SEED = 0;                             % seed for the random number generator (change it to get a different, but still reproducible, result)
P.bootstrap_optimal_k.N_BOOT = 1000;                            % number of bootstrap resamples (>= 1000 recommended)
P.bootstrap_optimal_k.VOTE_THRESHOLD = 0.60;                    % candidate K must win in >= 60% of bootstraps (the higher the value, the more conservative)
P.bootstrap_optimal_k.SIL_THRESHOLD = 0.30;                     % AND have a mean silhouette >= 0.25 (0.25 is a modest but non-trivial structure threshold) (the higher the value, the more conservative)

% permutation_test_clustering
P.permutation_test_clustering.RNG_SEED = 0;                     % seed for the random number generator (change it to get a different, but still reproducible, result)
P.permutation_test_clustering.N_PERM = 1000;                    % number of permutations (>= 1000 recommended)
P.permutation_test_clustering.ALPHA = 0.05;                     % significance threshold

% safe_linkage
P.safe_linkage.FALLBACK = 'average';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  7e) Network burst features

% NetworkBurstFeatures
P.NetworkBurstFeatures.BIN_REF = 0.1;                           % reference bin size to calculate the height of a network burst, in sec (default: 0.1s = 100 msec)
P.NetworkBurstFeatures.COUNT_WINDOW = 60;                       % length of the windows to count the network bursts for the Fano factor of the counts, in sec (default: 60)

% NetworkBurstStructure
P.NetworkBurstStructure.BIN_SIZE_POP = 30;                      % bin size to analyze the population activity, in msec (default: 30)
P.NetworkBurstStructure.BIN_NEURON = 10;                        % bin size to analyze individual neuron activity, in msec (10)
P.NetworkBurstStructure.THRESHOLD_ACT = 0.4;                    % threshold to detect the onset; proportion of burst events with activity (0.4)
P.NetworkBurstStructure.THRESHOLD_ACT_TWIN = 0.25;              % threshold to detect the onset (backup option); proportion of burst events with activity across two consecutive bins (0.25)
P.NetworkBurstStructure.THRESHOLD_TAIL = 0.25;                  % threshold to detect the tail; proportion of the peak firing rate (0.25)
P.NetworkBurstStructure.MIN_BURST_PLOT = 3;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  8) Compile_MEAproject (per group)

% Compile_MEAproject
P.Compile_MEAproject.MIN_NETBURST = 1/200;                      % minimum rate of network burst events, in # per sec (default: 1/150); it is applied to each family of events
                                                                % (revised October 2026: 1/300 -> 1/200; before, this value was not used -only in the disabled block "if 0" below- and
                                                                %  "Compile_NetBurstResults" used its own value, 1/150). Below it, the metrics stay NaN, but the number of events is reported
P.Compile_MEAproject.MIN_SECONDARY_RATE = 0.5;                  % minimum average number of secondary peaks per primary event (default: 0.5)
P.Compile_MEAproject.MIN_PARENT = 0.25;                         % minimum proportion of primary peaks with secondary peaks (default: 0.25)
P.Compile_MEAproject.TREATMENT = {'BASELINE', 'Gabazine', 'CGP', 'DNQX', 'D-AP5', 'TTX', 'unknown'};    % labels of the treatments
P.Compile_MEAproject.RATIO_MATCH = 0.70;                        % minimum matching between the treatment labels and the filename
P.Compile_MEAproject.COL_COND = 4;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  3Brain stages (step by step only)

% ShapeSpikeTrainData_3Brain: settings (INPUTS section of the function)
P.ShapeSpikeTrainData_3Brain.suffix = '_Baseline.bxr';
P.ShapeSpikeTrainData_3Brain.option_explore_conditions = 0;     % =1, to check the original order of the conditions (these files might have been saved in no chronological order)
P.ShapeSpikeTrainData_3Brain.condition_order = [2;1;3;4];       % set the temporal order of the conditions
P.ShapeSpikeTrainData_3Brain.option_version = 2;                % =1, BrainWave v5.x or older
                                                                % =2, BrainWave v6.x

% ShapeSpikeTrainData_3Brain
P.ShapeSpikeTrainData_3Brain.MATRIX_SIZE = 64;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Linked parameters: {parameter, the same value as}. "run_pipeline" copies the second one into the first one (after the "changes")

P.links = {
    'familyParticipants.K_ONSET',                   'findChannel_in_Peak.K_ONSET';
    'familyParticipants.MIN_PARTICIPATION',         'findChannel_in_Peak.MIN_PARTICIPATION';
    'familyParticipants.MIN_PARTICIPATION_ABS',     'findChannel_in_Peak.MIN_PARTICIPATION_ABS';
    'Optimal_secondaryPeak.MIN_CHAN',               'NetworkBurstDetection.MIN_CHAN';
    'Optimal_secondaryPeak.MIN_CHAN_SECONDARY',     'NetworkBurstDetection.MIN_CHAN_SECONDARY';
    'NetworkBurstPlayers.MIN_PEAK2_GAP',            'NetworkBurstDetection.MIN_PEAK2_GAP';
    'NetworkBurstPlayers.GAP_WITH_PEAK',            'NetworkBurstDetection.GAP_WITH_PEAK'};
