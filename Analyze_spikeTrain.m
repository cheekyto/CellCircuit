function Analyze_spikeTrain(inputs)

% (optional) inputs:    structure to overwrite any variable of the sections INPUTS and PARAMETERS (e.g., inputs.folders, inputs.files_data); see "override_inputs"
%                       when it is not provided, the values defined in the code are used
% (revised October 2026: raster plots drawn with a single plot call; the dose loop no longer reuses "j"; see also "SpikeFeatures")

if ~parameter_set('active'), clc; end; close all;      % (revised October 2026) no "clc" while "run_pipeline" runs (it would erase its progress messages)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system = 3;                 % =1, Ubuntu; =2, Yang's PC; =3, Yang's MacBook; =0, do not modify the MATLAB path (e.g., test harness)
if nargin < 1, inputs = struct(); end                                       % no inputs provided -> use the values defined in the code
if isfield(inputs, 'system'), system = inputs.system; end                   % the system can also be provided as an input
switch system
    case 1      % Anlysis PCP - Ubuntu
        addpath(genpath('/home/cortex/Documents/Yang'));
    case 2      % Yang's PC - GitHub
        addpath(genpath('C:\Users\Yang\OneDrive - Harvard University\Documents\GitHub\MEA-analysis\'));
    case 3      % Yang's MacBook
        addpath(genpath('/Users/sungminyang/Documents/GitHub/CellCircuit'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folders = [
"/path/to/recording_folder"];
files_data = [
"recording_name_curated"];

condition_analysis = 1;                             % list of recording condition to analyze
                                                    % =-1, all the recording conditions 
window = 600;                                       % the analysis time, at the end of each recording condition, in sec (default: 300 sec)
                                                    % =-1, the entire recording time, for each recording condition                                                   
opt_rerun = 1;                                       % =1, it is a rerun, an output file already exists and the new data will be appended
                                                    % =0, it is the first time the analysis of this experiment is run 
opt_distance = 0;                                   % option to remove primary network burst events that are too close to the previous one

opt_clusterAgain = 0;                               % option to re-run the clustering of network bursting events                 

BURST_DETECT_METHOD = 1;                            % algorithm to be used for individual burst detection
                                                    %  =1       Maximal Interval 
NETWORK_BURST_DETECT_METHOD = 3;                    % algorithm to be used for newtork burst detection
                                                    %  =1       Population firing rate
                                                    %  =2       looking for peaks, combining the regular population firing rate and the normalized population firing; and using a progressive narrowing of the binSize
                                                    %  =3       option 2 + clustering
NETWORK_BURST_ANALYSIS_METHOD = 2;                  % algorithm to be used for network burst analysis
                                                    %  =1       single family of network bursting, with a single set of participating channels
                                                    %  =2       cluster analysis to identify different families of network bursts
opt_checking = 0;                             
print_Plots = 1;                                    % option to save the rater plot of all neurons
plot_option = 0;
plot_network_burst = 1;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% matlab file sprintf('%s\\%s\\%s_Analyzed.mat', main_folder, project_folder, file_data), containing the following variables
%  1) spikeTimestamp, 2) PerCondition, 3) NetworkBurstDetected
% PerCondition(cond) has the following fields (per each recording condition): 
%   1) spikeTrain, 2) SpikeFeature, 3) BurstFeature, 4) BurstDetected, 5) BurstFeature_Population, 6) TonicFeature_Population, 7) NetworkBurstFeature, 8) NetworkBurstStructure, 9) condition, 10) window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                        % factor to rescale from msec to sec
K_SEC_MIN = 1/60;                           % factor to rescale from sec to min
MIN_N_BURST = 5;                            % minimum number of burst events to be considered a bursting neuron (default: 5)

MIN_WINDOW = 0.4;                           

RATIO_MATCH = 0.60;                         % minimum matching between the treatment labels and the filename
TREATMENT{1} = 'CGP';
TREATMENT{2} = 'Gabazine';
TREATMENT{3} = 'DNQX';
TREATMENT{4} = 'D-AP5';
TREATMENT{5} = 'TTX';

apply_parameters('Analyze_spikeTrain');                                     % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")
override_inputs(inputs);                    % overwrite the INPUTS / PARAMETERS above with the ones provided as input (if any)
n_treatment = size(TREATMENT,2);                     % (revised October 2026: after the parameter set and the inputs, which can change TREATMENT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check the inputs (folder + files)

n_recordings = numel(files_data);                                               % total number of files to analyze
for k=1:n_recordings                                                            % across all files
    if ~isfile(fullfile(folders(k), sprintf('%s_QCmanual.mat', files_data(k)))) % check the exitence of the file
        error("The file %s (index: %g) can not be found!! ", files_data(k), k)
    end
end

for k=1:n_recordings                                                        % across all recording files
    folder = folders(k);                                                    % path to locate the input file k
    file_data = files_data(k)                                               % file under analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the pre-processed data

if opt_rerun
    file_output = fullfile(folder, sprintf('%s_Analyzed.mat',file_data));
    load(file_output, 'spikeTimestamp','PerCondition','n_clusterOriginal','n_cluster_Kilosort','n_cluster_3Brain');
    file_manual = fullfile(folder, sprintf('%s_QCmanual.mat',file_data));
    load(file_manual, 'rec_conditions');
    fprintf("The recording %s will be analyzed again and the results will be appended to the existing 'Analyzed' file \n", file_data)
else
    file_manual = fullfile(folder, sprintf('%s_QCmanual.mat',file_data));
    load(file_manual, 'spikeTimestamp','rec_conditions','n_clusterOriginal','n_cluster_Kilosort','n_cluster_3Brain');
    fprintf("The recording %s will be analyzed for the first time and a new result file will be created \n", file_data)
end
n_clusters = size(spikeTimestamp,2);
n_conditions = size(rec_conditions,1);
n_condition_analysis = numel(condition_analysis);

fprintf("The total number of final good clusters in %s is %g \n", file_data, n_clusters)

condition_analysis0 = condition_analysis;
if n_condition_analysis == 1 
    if condition_analysis == -1
        condition_analysis0 = 1:n_conditions;
        n_condition_analysis = n_conditions;
        fprintf("All %g recording conditions will be analyzed \n", n_condition_analysis)  
    else
        fprintf("% g recording condition will be analyzed (out of  %g  total) \n", n_condition_analysis ,n_conditions)
    end
elseif n_condition_analysis > 0
    if n_condition_analysis > n_conditions
        error('There is NOT that many recording conditions')
    else
        fprintf("% g recording conditions will be analyzed (out of  %g  total) \n", n_condition_analysis ,n_conditions)
    end
else
    error('Please check the list of recording conditions to be analyzed')
end

condition_file = strings(1,n_conditions);                                   % file name for the condition  
condition = strings(1,n_conditions);                                        % name of the condition    
dose = strings(1,n_conditions);

tic
for j=1:n_condition_analysis                                               % across all recording conditions to be analyzed

    cond = condition_analysis0(j);
    if window==-1                                                           % the entire recording window
        t_window = rec_conditions{cond,2};
    else
        t_window = rec_conditions{cond,2};
        if t_window(1)<t_window(2)-window                       
            t_window(1) = t_window(2) - window;                             % if there is enough room, it will analyze "window" (at the end of the condition)                          
        else   
            if t_window(2)-t_window(1) >= window*MIN_WINDOW
                t_window(1) = t_window(1);
            else
                fprintf('\n\n The analysis window is larger than the recording length \n')
                t_window(1) = nan;            
                t_window(2) = nan;            
            end
        end
    end

    spikeTrain = Reshape_spikeTimestamp(spikeTimestamp, t_window);          % reshape the structure spikeTimestamp for the specified time window
    aux_cluster = split(rec_conditions{cond,1},".b");                               % assuming the raw data file was a .brw or .bxr file
    condition_file(cond) = aux_cluster(1);                                          % file name for the condition  
    pat = caseInsensitivePattern("_" + digitsPattern + ("u"|"m"|"n") + "M");
    if cond==1                                                              % assuming that the first condition is the baseline
        condition(cond) = 'Baseline';                                       % name of the condition    
    else
        aux_condition0 = rec_conditions{1,1};                               % assuming that the first condition is the baseline
        aux0 = split(aux_condition0,".b");
        condition0 = aux0(1);                                               % main name for the experiment
        aux_condition = extractAfter(rec_conditions{cond,1},condition0);
        aux_label = split(aux_condition,".b");                              % specific label for the condition, contained in the file name
        condition_label = [];
        for k=1:n_treatment
            if contains(convertCharsToStrings(aux_label{1}),lavelVAR(TREATMENT{k},RATIO_MATCH),'IgnoreCase',true)
                condition_label = horzcat(condition_label,'+');
                condition_label = horzcat(condition_label,TREATMENT{k});
            end
        end
        if size(condition_label)==0
            condition(cond) = horzcat('unknown', aux_label{1});
        else
            condition(cond) = condition_label;
            aux = extractAfter(extract(aux_condition, pat),"_");
            if size(aux,1)>1
                dose_aux = [];
                for jd=1:size(aux,1)                                        % (revised October 2026: "jd" instead of "j", the index of the loop across conditions)
                    dose_aux = strcat(dose_aux, string(aux{jd}));
                end
                dose(cond) = dose_aux; 
            else
                if size(aux,1)>0
                    dose(cond) = string(aux);
                end
            end
        end
        clear condition_label 
    end

    fprintf("\nAnalyzing the condition %s \n", string(condition_file(cond)))
    fprintf("It will be saved as the condition:  %s (%s) \n", string(condition(cond)), dose(cond))
    fprintf("(the temporal window for analysis is %g to %g sec) \n", t_window(1), t_window(2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyze the spiking activity

SpikeFeature = SpikeFeatures(spikeTrain, t_window);                         % call the function "SpikeFeatures" to analyze the recording condition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyze the bursting activity (population assessment)

nBurst = nan(n_clusters,1);       BurstDetected = cell(n_clusters,1);   BurstFeature = cell(n_clusters,1);
for i=1:n_clusters                                                          % across all clusters
    BurstDetected{i} = BurstDetection(spikeTrain(i).timeStamp, BURST_DETECT_METHOD);    % detect the bursts for cluster i
    BurstFeature{i} = BurstFeatures(BurstDetected{i}, t_window);            % characterize the bursts for cluster i (revised October 2026: with the analysis window, to count the bursts in windows)
    nBurst(i) = BurstFeature{i}.nBurst;                                     % total number of burst events identified for cluster i
end

aux_cluster = find(nBurst>MIN_N_BURST);                                             % index of clusters with number of burst events above threshold
aux_mean = mean(nBurst(aux_cluster));                                               % mean number of burst events, only considering the clusters with at least a few events
aux_std = std(nBurst(aux_cluster));                                                 % standard deviation for number of burst events
nBurst_threshold = max(MIN_N_BURST, aux_mean-2*aux_std);                    % set the threshold for number of burst events to classify the firing activity of each cluster as tonic vs phasic
bursting_idx = find(nBurst>nBurst_threshold);                               % index of clusters with phasic firing activity
n_bursting = numel(bursting_idx);                                           % total number of clusters with phasic activity

if numel(n_bursting)>0                                                      % NOTE (October 2026): numel(n_bursting) is always 1, so this block always runs; when no cluster 
                                                                            % is bursting, the population metrics are NaN (and all the fields exist, as expected by "Compile_MEAproject")
                                                                            % variables to temporaly collect the data for individual clusters and then calculate population metrics
    aux_nSpikes = nan(n_bursting,1); aux_duration = nan(n_bursting,1); aux_amplitude = nan(n_bursting,1); aux_amplitudeFF = nan(n_bursting,1); aux_burstRate = nan(n_bursting,1);
    aux_burstRateFF = nan(n_bursting,1); aux_spikingRate_nonBurst = nan(n_bursting,1); aux_spikingRateFF_nonBurst = nan(n_bursting,1); aux_nSpikes_nonBurst = nan(n_bursting,1);
    % (revised October 2026) variability metrics: CV and Fano factor in two forms (CV^2 and variance/mean) of the intervals, and Fano factor of the burst counts
    aux_amplitudeCV = nan(n_bursting,1); aux_amplitudeFFvm = nan(n_bursting,1); aux_burstRateCV = nan(n_bursting,1); aux_burstRateFFvm = nan(n_bursting,1);
    aux_burstFFcount = nan(n_bursting,1); aux_burstMeanCount = nan(n_bursting,1); aux_spikingRateCV_nonBurst = nan(n_bursting,1); aux_spikingRateFFvm_nonBurst = nan(n_bursting,1);

    for i=1:n_bursting                                                      % across all clusters with phasic firing activity
        idx = bursting_idx(i);                                              % index of the cluster i in the original array; collect the data for cluster i
        aux_nSpikes(i) = BurstFeature{idx}.neuronMetric.nSpikes.norm;       % normalized number of spikes within bursts
        aux_duration(i) = BurstFeature{idx}.neuronMetric.duration.mean;     % mean burst duration
        aux_amplitude(i) = BurstFeature{idx}.neuronMetric.isi.mean;         % mean isi for spikes within bursts
        aux_amplitudeCV(i) = BurstFeature{idx}.neuronMetric.isi.cv;         % (revised October 2026) mean coefficient of variation for isi of spikes within bursts
        aux_amplitudeFF(i) = BurstFeature{idx}.neuronMetric.isi.FF_cv2;     % Fano factor for isi of spikes within bursts, as CV^2
        aux_amplitudeFFvm(i) = BurstFeature{idx}.neuronMetric.isi.FF_varmean;  % Fano factor for isi of spikes within bursts, as variance / mean (msec)
        aux_burstRate(i) = BurstFeature{idx}.neuronMetric.ibi.mean;         % mean inter-burst interval    
        aux_burstRateCV(i) = BurstFeature{idx}.neuronMetric.ibi.cv;         % (revised October 2026) coefficient of variation for ibi
        aux_burstRateFF(i) = BurstFeature{idx}.neuronMetric.ibi.FF_cv2;     % Fano factor for ibi, as CV^2
        aux_burstRateFFvm(i) = BurstFeature{idx}.neuronMetric.ibi.FF_varmean;   % Fano factor for ibi, as variance / mean (sec)
        aux_burstFFcount(i) = BurstFeature{idx}.neuronMetric.count.FFcount;     % Fano factor of the number of bursts per window
        aux_burstMeanCount(i) = BurstFeature{idx}.neuronMetric.count.meanCount; % mean number of bursts per window
        aux_spikingRate_nonBurst(i) = BurstFeature{idx}.neuronMetric.nonBursting.isi.mean;          % mean isi for spikes outside bursts
        aux_spikingRateCV_nonBurst(i) = BurstFeature{idx}.neuronMetric.nonBursting.isi.cv;          % (revised October 2026) coefficient of variation for isi of spikes outside bursts
        aux_spikingRateFF_nonBurst(i) = BurstFeature{idx}.neuronMetric.nonBursting.isi.FF_cv2;      % Fano factor for isi of spikes outside bursts, as CV^2
        aux_spikingRateFFvm_nonBurst(i) = BurstFeature{idx}.neuronMetric.nonBursting.isi.FF_varmean;    % Fano factor for isi of spikes outside bursts, as variance / mean (msec)
        aux_nSpikes_nonBurst(i) = BurstFeature{idx}.neuronMetric.nonBursting.nSpikes_norm ;         % normalized number of spikes outside bursts
    end
                                                                                                                % calculate the population metrics
    BurstFeature_Population.burst.neurons = n_bursting/n_clusters;                                              % normalized number of clusters with phasic firing activity
    BurstFeature_Population.burst.nSpikes.norm = mean(aux_nSpikes, 'omitnan');                                  % mean normalized number of spikes within bursts
    BurstFeature_Population.burst.duration.mean = mean(aux_duration, 'omitnan');                                % mean burst duration, in sec
    BurstFeature_Population.burst.duration.cv = std(aux_duration, 'omitnan')/mean(aux_duration, 'omitnan');     % coefficient of variation for burst duration
    BurstFeature_Population.burst.amplitude.mean = 1/(mean(aux_amplitude, 'omitnan')*K_MSEC_SEC);               % mean spiking frequency within bursts, in Hz
    BurstFeature_Population.burst.amplitude.cv = std(aux_amplitude, 'omitnan')/mean(aux_amplitude, 'omitnan');  % coefficient of variation (across clusters) for spiking frequency
    BurstFeature_Population.burst.amplitude.cvISI = mean(aux_amplitudeCV, 'omitnan');                           % (revised October 2026) mean CV of the isi within bursts (mean across clusters)
    BurstFeature_Population.burst.amplitude.FFisi_cv2 = mean(aux_amplitudeFF, 'omitnan');                       % mean Fano factor of the isi within bursts, as CV^2 (mean across clusters; before, "fanoFactor")
    BurstFeature_Population.burst.amplitude.FFisi_varmean = mean(aux_amplitudeFFvm, 'omitnan');                 % mean Fano factor of the isi within bursts, as variance / mean, in msec
    BurstFeature_Population.burst.rate.mean = 1/(mean(aux_burstRate, 'omitnan')*K_SEC_MIN);                     % mean bursting rate, in 1/min
    BurstFeature_Population.burst.rate.cv = std(aux_burstRate, 'omitnan')/mean(aux_burstRate, 'omitnan');       % coefficient of variation for bursting rate
    BurstFeature_Population.burst.rate.cvIBI = mean(aux_burstRateCV, 'omitnan');                                % (revised October 2026) mean CV of the ibi (mean across clusters; "cv" above is across clusters)
    BurstFeature_Population.burst.rate.FFibi_cv2 = mean(aux_burstRateFF, 'omitnan');                            % mean Fano factor of the ibi, as CV^2 (before, "fanoFactor")
    BurstFeature_Population.burst.rate.FFibi_varmean = mean(aux_burstRateFFvm, 'omitnan');                      % mean Fano factor of the ibi, as variance / mean, in sec
    BurstFeature_Population.burst.rate.FFcount = mean(aux_burstFFcount, 'omitnan');                             % mean Fano factor of the number of bursts per window
    BurstFeature_Population.burst.rate.meanCount = mean(aux_burstMeanCount, 'omitnan');                         % mean number of bursts per window
    
    BurstFeature_Population.nonburst.nSpikes.norm = mean(aux_nSpikes_nonBurst, 'omitnan');                      % mean normalized number of spikes outside bursts 
    BurstFeature_Population.nonburst.spikeRate.mean = 1/(mean(aux_spikingRate_nonBurst,'omitnan')*K_MSEC_SEC);  % mean spiking frequency outside bursts, in Hz
    BurstFeature_Population.nonburst.spikeRate.cv = std(aux_spikingRate_nonBurst,'omitnan')/mean(aux_spikingRate_nonBurst,'omitnan');   % coefficient of variation (across clusters) for spiking frequency
    BurstFeature_Population.nonburst.spikeRate.cvISI = mean(aux_spikingRateCV_nonBurst, 'omitnan');             % (revised October 2026) mean CV of the isi outside bursts (mean across clusters)
    BurstFeature_Population.nonburst.spikeRate.FFisi_cv2 = mean(aux_spikingRateFF_nonBurst, 'omitnan');         % mean Fano factor of the isi outside bursts, as CV^2 (before, "fanoFactor")
    BurstFeature_Population.nonburst.spikeRate.FFisi_varmean = mean(aux_spikingRateFFvm_nonBurst, 'omitnan');   % mean Fano factor of the isi outside bursts, as variance / mean, in msec

else
    BurstFeature_Population.burst.neurons = 0;                              % normalized number of clusters with phasic firing activity
end

if cond==1
    nBurst0 = nBurst;
    nBurst_threshold0 = nBurst_threshold;
    n_bursting0 = n_bursting;
    aux_burstRate0 = aux_burstRate;
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyze the activity of neurons with tonic firing (population assessment)

nonBursting_idx = find(nBurst<=nBurst_threshold);                           % index of clusters with tonic firing activity 
n_nonBursting = numel(nonBursting_idx);                                     % total number of clusters with tonic activity

if numel(n_nonBursting)>0                                                   % NOTE (October 2026): always true (see above); when no cluster is tonic, the metrics are NaN
    aux_spikingRate_nonBurst = nan(n_nonBursting,1);                        % variables to temporaly collect the data for individual clusters and then calculate population metrics
    aux_spikingRateFF_nonBurst = nan(n_nonBursting,1);
    aux_spikingRateCV_nonBurst = nan(n_nonBursting,1);    aux_spikingRateFFvm_nonBurst = nan(n_nonBursting,1);      % (revised October 2026)

    for i=1:n_nonBursting                                                   % across all clusters with tonic firing activity
        idx = nonBursting_idx(i);                                           % index of the cluster i in the original array
                                                                            % collect the data for cluster i
        aux_spikingRate_nonBurst(i) = BurstFeature{idx}.neuronMetric.nonBursting.isi.mean;          % mean isi for spikes outside bursts 
        aux_spikingRateCV_nonBurst(i) = BurstFeature{idx}.neuronMetric.nonBursting.isi.cv;          % (revised October 2026) coefficient of variation for isi of spikes outside bursts
        aux_spikingRateFF_nonBurst(i) = BurstFeature{idx}.neuronMetric.nonBursting.isi.FF_cv2;      % Fano factor for isi of spikes outside bursts, as CV^2
        aux_spikingRateFFvm_nonBurst(i) = BurstFeature{idx}.neuronMetric.nonBursting.isi.FF_varmean;    % Fano factor for isi of spikes outside bursts, as variance / mean (msec)
    end

    TonicFeature_Population.spikeRate.mean = 1/(mean(aux_spikingRate_nonBurst,'omitnan')*K_MSEC_SEC);                          % mean spiking frequency outside bursts, in Hz      
    TonicFeature_Population.spikeRate.cv = std(aux_spikingRate_nonBurst,'omitnan')/mean(aux_spikingRate_nonBurst,'omitnan');   % coefficient of variation (across clusters) for spiking frequency
    TonicFeature_Population.spikeRate.cvISI = mean(aux_spikingRateCV_nonBurst, 'omitnan');                                     % (revised October 2026) mean CV of the isi (mean across clusters)
    TonicFeature_Population.spikeRate.FFisi_cv2 = mean(aux_spikingRateFF_nonBurst, 'omitnan');                                 % mean Fano factor of the isi, as CV^2 (before, "fanoFactor")
    TonicFeature_Population.spikeRate.FFisi_varmean = mean(aux_spikingRateFFvm_nonBurst, 'omitnan');                           % mean Fano factor of the isi, as variance / mean, in msec
    TonicFeature_Population.neurons = n_nonBursting/n_clusters;             % normalized number of clusters with tonic firing activity
else
    TonicFeature_Population.neurons = 0;                                    % normalized number of clusters with tonic firing activity
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyze the network bursting

if opt_clusterAgain == 0
    NetworkBurstDetected = NetworkBurstDetection(spikeTrain, BurstDetected, t_window, NETWORK_BURST_DETECT_METHOD, opt_distance, plot_network_burst);
else
    file_input = fullfile(folder, sprintf('%s_Analyzed.mat',file_data));
    load(file_input, 'PerCondition');
    NetworkBurstDetected = PerCondition(cond).NetworkBurstDetected;
    flag2 = NetworkBurstDetected.flag2;
%    flag2 = 0;
    NetworkBurstDetected.clustering = clusterNetburst(NetworkBurstDetected, spikeTrain, BurstDetected, t_window, flag2, 1);
end

% (revised October 2026) report the number of families of network burst events and the quality of the clustering, to know how much weight to put on it
% (quality_clustering, saved as NetworkBurstDetected.clustering.clusteringQuality; see "cluster_neuronal_events")
if isstruct(NetworkBurstDetected.clustering) && isfield(NetworkBurstDetected.clustering, 'clusteringQuality')
    fprintf("\n Clustering of the network burst events: %g families, quality_clustering = %g \n", NetworkBurstDetected.clustering.n_cluster, NetworkBurstDetected.clustering.clusteringQuality);
    fprintf("   (2 = reliable; 1 = it seems to work, but the bootstrap and the permutation test do not give a clear answer; 0 = not significant, single family; NaN = not run) \n");
else
    fprintf("\n No clustering of the network burst events (no network bursting): quality_clustering = NaN \n");
end

nNetBurst1 = NetworkBurstDetected.nBurst;
nNetBurst2 = NetworkBurstDetected.Secondary.nBurst;
if ~isnan(nNetBurst1) && nNetBurst1 > 0
    switch NETWORK_BURST_ANALYSIS_METHOD
        
        case 1     % single family of network bursting, with a single set of participating channels
            netburst_member = ones(nNetBurst1,1);
            channels_part = NetworkBurstDetected.Channels(:,1);
            NetworkBurstFeature{1} = NetworkBurstFeatures(NetworkBurstDetected, spikeTrain, BurstDetected, t_window/K_MSEC_SEC, netburst_member, channels_part, nNetBurst1, 0);   % (revised October 2026) the argument nBurst1 was missing (crash)
            NetworkBurstFeature0 = NetworkBurstFeature{1};                  % (revised October 2026) all the primary events (single family); it was not defined -> crash when saving
            if nNetBurst2 > 0 
                netburst_member2 = ones(nNetBurst2,1);
                channels_part2 = NetworkBurstDetected.Secondary.Channels(:,1);
                NetworkBurstFeature{2} = NetworkBurstFeatures(NetworkBurstDetected.Secondary, spikeTrain, BurstDetected, t_window/K_MSEC_SEC, netburst_member2, channels_part2, NetworkBurstDetected.nBurst_level1, 1);   % (revised October 2026) the argument nBurst1 was missing (crash); level-1 primary events
            end

        case 2     % cluster analysis to identify different families of network bursts
            
            channels_part = NetworkBurstDetected.Channels(:,1);
            cluster1 = NetworkBurstDetected.clustering.label(1:nNetBurst1,1);
            aux_cluster = unique(cluster1);
            nFamilies = numel(aux_cluster);
            for k=1:nFamilies
                netburst_member_k = (cluster1 == aux_cluster(k));
                % (revised October 2026) participating units of the family k: the participating units of all the primary events that fire in
                % more than 20% of the family's events and in at least 2 of them (see "familyParticipants"); before, every family used all the 
                % participating units ("channels_part"), so "nNeuron.normBurst" was always 1. Their number is saved in "nNeuron.Burst"
                channels_part_k = familyParticipants(NetworkBurstDetected, BurstDetected, netburst_member_k, channels_part);
                NetworkBurstFeature{k} = NetworkBurstFeatures(NetworkBurstDetected, spikeTrain, BurstDetected, t_window/K_MSEC_SEC, netburst_member_k, channels_part_k, nNetBurst1, 0);
            end
            netburst_member = ones(nNetBurst1,1);
            NetworkBurstFeature0 = NetworkBurstFeatures(NetworkBurstDetected, spikeTrain, BurstDetected, t_window/K_MSEC_SEC, netburst_member, channels_part, nNetBurst1, 0);

            if nNetBurst2 > 0 
                netburst_member2 = ones(nNetBurst2,1);
                channels_part2 = NetworkBurstDetected.Secondary.Channels(:,1);
                NetworkBurstFeature{k+1} = NetworkBurstFeatures(NetworkBurstDetected.Secondary, spikeTrain, BurstDetected, t_window/K_MSEC_SEC, netburst_member2, channels_part2, NetworkBurstDetected.nBurst_level1, 1);  % (revised October 2026) number of LEVEL-1 primary events (the parents of the secondary events), instead of all primary events
            end
    end

    netburst_member = ones(nNetBurst1,1);                                   % to run "NetworkBurstStructure"
    units_participant = NetworkBurstDetected.Channels(:,1);                 % to run "NetworkBurstStructure": participating units of the primary events
                                                                            % (revised October 2026: only them, the same units as "NetworkBurstFeature0"; 
                                                                            %  before, also the participating units of the secondary events)
else
   netburst_member = [];
   channels_part = ones(n_clusters,1);
   NetworkBurstFeature{1} = NetworkBurstFeatures(NetworkBurstDetected, spikeTrain, BurstDetected, t_window/K_MSEC_SEC, netburst_member,channels_part, nNetBurst1, 0);
   units_participant = [];
   NetworkBurstFeature0 = [];
end

% Only assess the network burst structure for all primary events (together, disregarding the clustering)
NetworkBurstStruct = NetworkBurstStructure(NetworkBurstDetected, spikeTrain, netburst_member, units_participant, plot_network_burst);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare the output

PerCondition(cond).spikeTrain = spikeTrain;
PerCondition(cond).SpikeFeature = SpikeFeature;
PerCondition(cond).BurstFeature = BurstFeature;
PerCondition(cond).BurstDetected = BurstDetected;
PerCondition(cond).BurstFeature_Population = BurstFeature_Population;
PerCondition(cond).TonicFeature_Population = TonicFeature_Population;
PerCondition(cond).NetworkBurstFeature = NetworkBurstFeature;
PerCondition(cond).NetworkBurstStructure = NetworkBurstStruct;
PerCondition(cond).NetworkBurstFeature0 = NetworkBurstFeature0;
PerCondition(cond).NetworkBurstDetected = NetworkBurstDetected;
PerCondition(cond).condition = condition(cond);
PerCondition(cond).dose = dose(cond);
PerCondition(cond).window = t_window;
PerCondition(cond).file = rec_conditions{cond,1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the raster plot and population firing for all participating neurons

if print_Plots 
        
    figure('Name','Raster Plot and Population Firing','Color','w');
    colors_event = {'r', 'b', 'c', 'g', 'm'};

    if ~isnan(nNetBurst1) && nNetBurst1>0 && numel(NetworkBurstDetected.Channels(:,1))>0
    
        subplot(2,1,1);
        data = spikeTimestamp(units_participant);                           % raster plot of the participating units, with a single plot call (revised October 2026)
        all_rows = repelem((1:numel(data))', cellfun(@(t) size(t,1), {data.timeStamp})') + 1;   % row of each spike (same as "ones(n_unit,1)+i")
        plot(vertcat(data.timeStamp)*K_MSEC_SEC, all_rows,'|',"MarkerEdgeColor",[0 0 0]); hold on
        i = numel(units_participant);
        ylabel('Neurons');  title('Selected units - Participating in network burst');  xlim([t_window(1) t_window(2)]); ylim([1 i+2])
    
        subplot(2,1,2);
        bar(NetworkBurstDetected.time(:), NetworkBurstDetected.spikesBin(:,1)); hold on;
        xlabel('Time (sec)');   ylabel('Firing rate (norm.)');  xlim([t_window(1) t_window(2)])

        switch NETWORK_BURST_ANALYSIS_METHOD
            case 1
                for i=1:nNetBurst1
                    plot(NetworkBurstDetected.IndivBurst(i).peak, NetworkBurstDetected.IndivBurst(i).peakAmp*1.1,'o', 'MarkerFaceColor',colors_event{1}); hold on
                end
                if nNetBurst2 > 0
                    for i=1:nNetBurst2
                        plot(NetworkBurstDetected.Secondary.IndivBurst(i).peak, NetworkBurstDetected.Secondary.IndivBurst(i).peakAmp*1.5,'o', 'MarkerFaceColor',colors_event{2}); hold on
                    end
                end
                title(sprintf('%g primary events + %g secondary events', nNetBurst1, nNetBurst2));

            case 2
                    
                for k=1:nFamilies
                    label1 = NetworkBurstDetected.clustering.label(1:nNetBurst1,1);
                    member_k = find((label1 == aux_cluster(k)));
                    for i=1:numel(member_k)
                        idx = member_k(i);
                        plot(NetworkBurstDetected.IndivBurst(idx).peak, NetworkBurstDetected.IndivBurst(idx).peakAmp*1.1,'o', 'MarkerFaceColor',colors_event{k}); hold on
                    end
                end
                if nNetBurst2 > 0
                    for i=1:nNetBurst2
                        plot(NetworkBurstDetected.Secondary.IndivBurst(i).peak, NetworkBurstDetected.Secondary.IndivBurst(i).peakAmp*1.5,'o', 'MarkerFaceColor',colors_event{k+1}); hold on
                    end
                end
                title(sprintf('%g primary clusters + %g secondary events', nFamilies, nNetBurst2));

            end
    else
        subtitle('NO network bursting', 'Color','red');
    end  

    file_output = fullfile(folder, sprintf('%s_NetBurst_%s%s.fig', file_data, string(condition(cond)),dose(cond)));
    saveas(gcf, file_output)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the results

if opt_checking
    fprintf("\n Let me know when you are ready to move on... \n");                              
    pause();
end

clear NetworkBurstFeature NetworkBurstStruct NetworkBurstFeature0 NetworkBurstDetected;

end         % end of the for loop across recording conditions

file_output = fullfile(folder, sprintf('%s_Analyzed.mat',file_data));
if opt_rerun
    save(file_output, 'spikeTimestamp','PerCondition','n_clusterOriginal','n_conditions','n_cluster_Kilosort','n_cluster_3Brain',"-append");
else
    save(file_output, 'spikeTimestamp','PerCondition','n_clusterOriginal','n_conditions','n_cluster_Kilosort','n_cluster_3Brain');
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the results

if plot_option

    % Plot the number of burst events
    figure;
    h=histogram(nBurst0);        hold on
    plot([nBurst_threshold0 nBurst_threshold0], [0 max(h.Values)],':r','DisplayName','threshold');
    text(nBurst_threshold0, max(h.Values)/2,'threshold');
    xlabel('Number of burst events');       ylabel('Number of neurons');
    title(sprintf('Histogram: number of burst events during the first condition (%g bursting neurons out of %g neurons)',  n_bursting0, n_clusters));

    % Plot the inter-burst interval
    figure;
    histogram(1./aux_burstRate0*K_SEC_MIN);
    xlabel('Bursting Rate (#/min)');        ylabel('Number of neurons');
    title('Histogram: bursting rate during the first condition');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the raster plot and population firing for all neurons

LABEL_POS = -5;
if print_Plots && ~opt_rerun

    % ----- Raster plot of all active neurons -----------------
    figure;
    all_rows = repelem((1:n_clusters)', cellfun(@numel, {spikeTimestamp.timeStamp})');     % row of each spike (same as "ones(data.nUnits,1)+k-1"), single plot call (revised October 2026)
    plot(vertcat(spikeTimestamp.timeStamp)*K_MSEC_SEC, all_rows,'|',"MarkerEdgeColor",[0 0 0]);     hold on;
    k = n_clusters;
    ylim([0 k+1]);      xlabel('Time (sec)');       ylabel('Neuron')
    title(sprintf('Raster Plot (%g neurons)', n_clusters));

    for cond=1:n_conditions
        t_ini = rec_conditions{cond,2}(1);
        plot([t_ini t_ini], [0 n_clusters], '--r');     hold on
        text(t_ini, LABEL_POS, string(condition(cond)),'Color', 'red');     hold on
    end
    t_fin = rec_conditions{cond,2}(2);
    plot([t_fin t_fin], [0 n_clusters], '--r');

    file_output = fullfile(folder, sprintf('%s_RasterPlot.fig', file_data));
    saveas(gcf, file_output)

end
clear('condition_analysis0', 'spikeTimestamp','PerCondition','n_clusterOriginal','n_cluster_Kilosort','n_cluster_3Brain','rec_conditions');
end

end