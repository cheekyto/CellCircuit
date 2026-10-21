function Compile_MEAproject(inputs)

% function to compile the results from a project based on data recorded in a MEA platform.
% Assessment of the spiking activity, the bursting pattern, and the network bursts
% Each recording was previously analyzed using the function "Analyze_spikeTrain" and the files corresponding to a single group were collected in a single folder
% (optional) inputs:    structure to overwrite any variable of the sections INPUTS and PARAMETERS (e.g., inputs.folderPath, inputs.group); see "override_inputs"
%                       when it is not provided, the values defined in the code are used

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

folderPath = '/path/to/project_folder';
group = 'groupA';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matlab file and Excel file containing the following tables (with population metrics for each recording):
% 1) Firing;    2) Bursting;    3) NETbursting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MIN_NETBURST = 1/200;                      % minimum rate of network burst events, in # per sec (default: 1/150); it is applied to each family of events
                                           % (revised October 2026: 1/300 -> 1/200; before, this value was not used -only in the disabled block "if 0" below- and 
                                           %  "Compile_NetBurstResults" used its own value, 1/150). Below it, the metrics stay NaN, but the number of events is reported
MIN_SECONDARY_RATE = 0.5;                  % minimum average number of secondary peaks per primary event (default: 0.5)
MIN_PARENT = 0.25;                         % minimum proportion of primary peaks with secondary peaks (default: 0.25)

K_MSEC_SEC = 1/1000;                        % factor to rescale from msec to sec

RATIO_MATCH = 0.70;                         % minimum matching between the treatment labels and the filename
TREATMENT{1} = 'BASELINE';
TREATMENT{2} = 'Gabazine';
TREATMENT{3} = 'CGP';
TREATMENT{4} = 'DNQX';
TREATMENT{5} = 'D-AP5';
TREATMENT{6} = 'TTX';
TREATMENT{7} = 'unknown';

COL_COND = 4;

apply_parameters('Compile_MEAproject');                                     % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")
override_inputs(inputs);                    % overwrite the INPUTS / PARAMETERS above with the ones provided as input (if any)
n_treatment = size(TREATMENT,2);                     % (revised October 2026: after the parameter set and the inputs, which can change TREATMENT)
thresholds = struct('MIN_NETBURST', MIN_NETBURST, 'MIN_SECONDARY_RATE', MIN_SECONDARY_RATE, 'MIN_PARENT', MIN_PARENT);   % (revised October 2026) used by "Compile_NetBurstResults"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify all the files

folder = fullfile(folderPath, group);
listing = dir(folder);

n_recordings=0;
n_condition=0;
for i=1:size(listing,1)
    if ~listing(i).isdir && endsWith(listing(i).name, '_Analyzed.mat')    % look for all the result files of "Analyze_spikeTrain" in the folder (revised October 2026: 
                                                                            % before, any name containing "Analyzed" -e.g., an Excel or figure file- was compiled again or crashed)
        n_recordings = n_recordings + 1;
        Recording{n_recordings,1} = extractBefore(listing(i).name,'_Analyzed');   % cell array with the name of all recording files found in the folder
        file_input = fullfile(folderPath, group, sprintf('%s_Analyzed.mat',Recording{n_recordings}));
        load(file_input,'n_conditions');
        n_condition = n_condition + n_conditions;
        clear n_conditions
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compile the results on spiking activity (general), bursting activity and network bursting
% (revised October 2026) each file is loaded only once (before, once per table, plus once per condition and option in "Compile_NetBurstResults");
% the conditions that were not analyzed yet (see "opt_rerun" in "Analyze_spikeTrain") are skipped with a warning (before, it crashed)

Organoid = cell(n_condition,1);     File = cell(n_condition,1);     Condition = cell(n_condition,1);    dose = cell(n_condition,1);
ttx_qc = nan(n_condition,1);
T_start = nan(n_condition,1); 
T_end = nan(n_condition,1); 
Number_Good_Clusters = nan(n_condition,1);
Number_Kilosort = nan(n_condition,1);
Number_3Brain = nan(n_condition,1);
Firing_Rate = nan(n_condition,1);
CV_ISI = nan(n_condition,1);     FFisi_CV2 = nan(n_condition,1);     FFisi_VarMean = nan(n_condition,1);       % (revised October 2026: they replace "Fano_Factor")
FFcount_Spikes = nan(n_condition,1);     MeanCount_Spikes = nan(n_condition,1);
Entropy = nan(n_condition,1);
LVR = nan(n_condition,1); 
mutualInfo = nan(n_condition,1);
sttc =  nan(n_condition,1);

NumberNorm_Bursting_Neurons = nan(n_condition,1);                          % variables to temporaly collect the data from each recording
NumberNorm_Spikes_Burst = nan(n_condition,1);
Rate_Burst = nan(n_condition,1);
RateCV_Burst = nan(n_condition,1);
CV_IBI = nan(n_condition,1);     FFibi_CV2 = nan(n_condition,1);     FFibi_VarMean = nan(n_condition,1);       % (revised October 2026: they replace "RateFF_Burst")
FFcount_Burst = nan(n_condition,1);      MeanCount_Burst = nan(n_condition,1);
Duration_Burst = nan(n_condition,1);
DurationCV_Burst = nan(n_condition,1);
Amplitude_Burst = nan(n_condition,1);
AmplitudeCV_Burst = nan(n_condition,1);
CV_ISI_Burst = nan(n_condition,1);       FFisi_CV2_Burst = nan(n_condition,1);       FFisi_VarMean_Burst = nan(n_condition,1);      % (revised October 2026: they replace "AmplitudeFF_Burst")

NumberNorm_Spikes_nonBurst = nan(n_condition,1);
Rate_Spiking_nonBurst = nan(n_condition,1);
RateCV_Spiking_nonBurst = nan(n_condition,1);
CV_ISI_nonBurst = nan(n_condition,1);    FFisi_CV2_nonBurst = nan(n_condition,1);    FFisi_VarMean_nonBurst = nan(n_condition,1);   % (revised October 2026: they replace "RateFF_Spiking_nonBurst")

NumberNorm_Tonic_Neurons = nan(n_condition,1);
Rate_Spiking_Tonic = nan(n_condition,1);
RateCV_Spiking_Tonic = nan(n_condition,1);
CV_ISI_Tonic = nan(n_condition,1);       FFisi_CV2_Tonic = nan(n_condition,1);       FFisi_VarMean_Tonic = nan(n_condition,1);      % (revised October 2026: they replace "RateFF_Spiking_Tonic")

NETbursting0 = [];
NETbursting1 = [];
NETbursting2 = [];

idx = 0;
for i=1:n_recordings                                                        % across all experiments
    file_input = fullfile(folderPath, group, sprintf('%s_Analyzed.mat', Recording{i}));
    load(file_input, 'PerCondition','n_conditions','n_cluster_Kilosort','n_cluster_3Brain');
    info = whos('-file', file_input, 'spikeTimestamp');                     % size of "spikeTimestamp", without loading all the spikes
    n_good_clusters = info.size(2);                                         % total number of good clusters (as "size(spikeTimestamp,2)")

    analyzed = false(n_conditions,1);                                       % =true, if the condition k was analyzed
    for k=1:min(n_conditions, numel(PerCondition))
        analyzed(k) = ~isempty(PerCondition(k).SpikeFeature);
    end
    if ~all(analyzed)
        fprintf("\n WARNING: %s has %g recording conditions, but only %g were analyzed (conditions skipped: %s) \n", Recording{i}, n_conditions, nnz(analyzed), mat2str(find(~analyzed)'));
    end

    ttx = 0;                                                                % =1, if the recording has a TTX condition (among the analyzed conditions)
    for k=find(analyzed)'
        if strcmp(PerCondition(k).condition, '+TTX')
            ttx = 1;
        end
    end

    for k=find(analyzed)'                                                   % across all (analyzed) conditions in experiment i
        idx = idx+1;
        Organoid{idx,1} = Recording{i};                                     % name of the experiment
        File{idx,1} = PerCondition(k).file;                                 % name of the file
        Condition{idx,1} = PerCondition(k).condition;                       % label of the experimental condition
        dose{idx,1} = PerCondition(k).dose;                                 % dose/concentration of the treatment drug
        ttx_qc(idx) = ttx;

        % spiking activity (general)
        T_start(idx) = PerCondition(k).window(1);                           % start of the analysis window, in sec
        T_end(idx) = PerCondition(k).window(2);                             % end of the analysis window, in sec
        Number_Good_Clusters(idx) = n_good_clusters;                        % total number of good clusters
        Number_Kilosort(idx) = n_cluster_Kilosort;                          % total number of good clusters coming from Kilosort dataset
        Number_3Brain(idx) = n_cluster_3Brain;                              % total number of good clusters coming from 3Brain dataset
        SpikeFeature = PerCondition(k).SpikeFeature;
        Firing_Rate(idx) = SpikeFeature.Population.firingRateISI;           % mean spiking frequency
        CV_ISI(idx) = SpikeFeature.Population.cv;                           % (revised October 2026) mean coefficient of variation of the isi
        FFisi_CV2(idx) = SpikeFeature.Population.FFisi_cv2;                 % mean Fano factor of the isi, as CV^2
        FFisi_VarMean(idx) = SpikeFeature.Population.FFisi_varmean;         % mean Fano factor of the isi, as variance / mean, in msec
        FFcount_Spikes(idx) = SpikeFeature.Population.FFcount;              % mean Fano factor of the spike counts (windows of COUNT_WINDOW sec, see "SpikeFeatures")
        MeanCount_Spikes(idx) = SpikeFeature.Population.meanCount;          % mean number of spikes per window
        Entropy(idx) = SpikeFeature.Population.entropy;                     % mean entropy
        LVR(idx) = SpikeFeature.Population.lvr; 
        mutualInfo(idx) = SpikeFeature.Population.mutualInfo;
        sttc(idx) =  SpikeFeature.Population.sttc;

        % bursting activity
        BurstFeature_Population = PerCondition(k).BurstFeature_Population;
        TonicFeature_Population = PerCondition(k).TonicFeature_Population;

        NumberNorm_Bursting_Neurons(idx) = BurstFeature_Population.burst.neurons;             % normalized number of clusters with phasic firing activity 
        NumberNorm_Tonic_Neurons(idx) = TonicFeature_Population.neurons;                      % normalized number of clusters with tonic firing activity

        if BurstFeature_Population.burst.neurons>0                                           % if there is no phasic clusters, the metrics stay NaN
            NumberNorm_Spikes_Burst(idx) = BurstFeature_Population.burst.nSpikes.norm;        % normalized number of spikes within bursts
            Rate_Burst(idx) = BurstFeature_Population.burst.rate.mean;                        % bursting rate, in 1/min
            RateCV_Burst(idx) = BurstFeature_Population.burst.rate.cv;                        % coefficient of variation for bursting rate
            CV_IBI(idx) = BurstFeature_Population.burst.rate.cvIBI;                           % (revised October 2026) mean coefficient of variation of the ibi
            FFibi_CV2(idx) = BurstFeature_Population.burst.rate.FFibi_cv2;                    % mean Fano factor of the ibi, as CV^2
            FFibi_VarMean(idx) = BurstFeature_Population.burst.rate.FFibi_varmean;            % mean Fano factor of the ibi, as variance / mean, in sec
            FFcount_Burst(idx) = BurstFeature_Population.burst.rate.FFcount;                  % mean Fano factor of the burst counts (windows of COUNT_WINDOW sec, see "BurstFeatures")
            MeanCount_Burst(idx) = BurstFeature_Population.burst.rate.meanCount;              % mean number of bursts per window
            Duration_Burst(idx) = BurstFeature_Population.burst.duration.mean;                % burst duration, in sec
            DurationCV_Burst(idx) = BurstFeature_Population.burst.duration.cv;                % coefficient of variation for burst duration
            Amplitude_Burst(idx) = BurstFeature_Population.burst.amplitude.mean;              % spiking frequency within bursts, in Hz
            AmplitudeCV_Burst(idx) = BurstFeature_Population.burst.amplitude.cv;              % coefficient of variation (across clusters) for spiking frequency
            CV_ISI_Burst(idx) = BurstFeature_Population.burst.amplitude.cvISI;                % (revised October 2026) mean coefficient of variation of the isi within bursts
            FFisi_CV2_Burst(idx) = BurstFeature_Population.burst.amplitude.FFisi_cv2;         % mean Fano factor of the isi within bursts, as CV^2
            FFisi_VarMean_Burst(idx) = BurstFeature_Population.burst.amplitude.FFisi_varmean; % mean Fano factor of the isi within bursts, as variance / mean, in msec

            NumberNorm_Spikes_nonBurst(idx) = BurstFeature_Population.nonburst.nSpikes.norm;  % normalized number of spikes outside bursts 
            Rate_Spiking_nonBurst(idx) = BurstFeature_Population.nonburst.spikeRate.mean;     % spiking frequency outside bursts, in Hz
            RateCV_Spiking_nonBurst(idx) = BurstFeature_Population.nonburst.spikeRate.cv;     % coefficient of variation (across clusters) for spiking frequency
            CV_ISI_nonBurst(idx) = BurstFeature_Population.nonburst.spikeRate.cvISI;          % (revised October 2026) mean coefficient of variation of the isi outside bursts
            FFisi_CV2_nonBurst(idx) = BurstFeature_Population.nonburst.spikeRate.FFisi_cv2;   % mean Fano factor of the isi outside bursts, as CV^2
            FFisi_VarMean_nonBurst(idx) = BurstFeature_Population.nonburst.spikeRate.FFisi_varmean;   % mean Fano factor of the isi outside bursts, as variance / mean, in msec
        end

        if TonicFeature_Population.neurons>0                                                % if there is no tonic clusters, the metrics stay NaN
            Rate_Spiking_Tonic(idx) = TonicFeature_Population.spikeRate.mean;                 % spiking frequency, in Hz      
            RateCV_Spiking_Tonic(idx) = TonicFeature_Population.spikeRate.cv;                 % coefficient of variation (across clusters) for spiking frequency
            CV_ISI_Tonic(idx) = TonicFeature_Population.spikeRate.cvISI;                      % (revised October 2026) mean coefficient of variation of the isi
            FFisi_CV2_Tonic(idx) = TonicFeature_Population.spikeRate.FFisi_cv2;               % mean Fano factor of the isi, as CV^2
            FFisi_VarMean_Tonic(idx) = TonicFeature_Population.spikeRate.FFisi_varmean;       % mean Fano factor of the isi, as variance / mean, in msec
        end

        % network bursting (see "Compile_NetBurstResults"; it receives "PerCondition", already loaded)
        % option = 0
        Compile_NetBurstResult = Compile_NetBurstResults(1, k, Recording{i}, PerCondition, ttx, 0, thresholds);
        NETbursting0 = vertcat(NETbursting0, Compile_NetBurstResult);
        % option = 1
        Compile_NetBurstResult = Compile_NetBurstResults(1, k, Recording{i}, PerCondition, ttx, 1, thresholds);
        NETbursting1 = vertcat(NETbursting1, Compile_NetBurstResult);
        % option = 2
        if PerCondition(k).NetworkBurstDetected.nBurst > 0
            n_cluster = numel(unique(PerCondition(k).NetworkBurstDetected.clustering.label(1:PerCondition(k).NetworkBurstDetected.nBurst,1)));  % number of families of primary events
            ttxQC_array = ones(n_cluster,1)*ttx;                            % (revised October 2026: before, "clustering.n_cluster", which also counts the secondary-only families)
        else
            n_cluster = 1;
            ttxQC_array = ttx;
        end
        Compile_NetBurstResult = Compile_NetBurstResults(n_cluster, k, Recording{i}, PerCondition, ttxQC_array, 2, thresholds);
        NETbursting2 = vertcat(NETbursting2, Compile_NetBurstResult);

        clear SpikeFeature BurstFeature_Population TonicFeature_Population Compile_NetBurstResult
    end
    
    clear n_conditions PerCondition n_cluster_Kilosort n_cluster_3Brain
end

Firing = table(Organoid,ttx_qc,File,Condition,dose,T_start,T_end,Number_Good_Clusters,Number_Kilosort,Number_3Brain,Firing_Rate,CV_ISI,FFisi_CV2,FFisi_VarMean,FFcount_Spikes,MeanCount_Spikes,Entropy,LVR,mutualInfo,sttc); % build a table with the population metrics for spiking activity (general)
Firing = Firing(1:idx,:);                                                   % (only the analyzed conditions)

Bursting = table(Organoid,File, ttx_qc, Condition, dose, NumberNorm_Bursting_Neurons, NumberNorm_Spikes_Burst, Rate_Burst, RateCV_Burst, CV_IBI, FFibi_CV2, FFibi_VarMean, FFcount_Burst, MeanCount_Burst, ...
                Duration_Burst, DurationCV_Burst, Amplitude_Burst, AmplitudeCV_Burst, CV_ISI_Burst, FFisi_CV2_Burst, FFisi_VarMean_Burst ... ...
                ,NumberNorm_Spikes_nonBurst, Rate_Spiking_nonBurst, RateCV_Spiking_nonBurst, CV_ISI_nonBurst, FFisi_CV2_nonBurst, FFisi_VarMean_nonBurst, NumberNorm_Tonic_Neurons, Rate_Spiking_Tonic, RateCV_Spiking_Tonic, CV_ISI_Tonic, FFisi_CV2_Tonic, FFisi_VarMean_Tonic);
Bursting = Bursting(1:idx,:);                                               % (only the analyzed conditions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compile the results on network bursting 

if 0
Number_NETburst_Events = nan(n_condition,1);                                % variables to temporaly collect the data from each recording
Number_NETburst_Neurons = nan(n_condition,1);
NumberNorm_NETburst_Neurons = nan(n_condition,1);
Number_Spikes_NETburst = nan(n_condition,1);
NumberNorm_Spikes_NETburst = nan(n_condition,1);
NumberNorm_SpikesTotal_NETburst = nan(n_condition,1);
Participation_NETburst = nan(n_condition,1);
ParticipationCV_NETburst = nan(n_condition,1);
Rate_NETburst = nan(n_condition,1);
RateCV_NETburst = nan(n_condition,1);
Duration_NETburst = nan(n_condition,1);
DurationCV_NETburst = nan(n_condition,1);
Tail_NETburst = nan(n_condition,1);
TailCV_NETburst = nan(n_condition,1);
Amplitude_NETburst = nan(n_condition,1);
AmplitudeCV_NETburst = nan(n_condition,1);
Distance_NETburst = nan(n_condition,1);
DistanceCV_NETburst = nan(n_condition,1);
Onset_NETburst = nan(n_condition,1);
OnsetCV_NETburst = nan(n_condition,1);
Peak_NETburst = nan(n_condition,1);
PeakCV_NETburst = nan(n_condition,1);
Correlation_NETburst = nan(n_condition,1);
tempStructNeuron_NETburst = nan(n_condition,1);
PeakFire_NETburst = nan(n_condition,1);
PeakFireNorm_NETburst = nan(n_condition,1);

Number_clusters_Netburst = nan(n_condition,1);

Number_SecondaryBurst_Neurons = nan(n_condition,1);
NumberNorm_SecondaryBurst_Neurons = nan(n_condition,1);
Amplitude_SecondaryBurst = nan(n_condition,1);
ParentNorm_SecondaryBurst = nan(n_condition,1);
NumberNorm_SecondaryBurst = nan(n_condition,1);

idx = 0;
for i=1:n_recordings                                                                % across all recordings
    file_input = fullfile(folderPath, group, sprintf('%s_Analyzed.mat', Recording{i}));
    load(file_input, 'PerCondition', 'n_conditions');

    for k=1:n_conditions                                                    % across all conditions in experiment i
        idx = idx+1;
        Organoid{idx,1} = Recording{i};                                     % name of the experiment
        File{idx,1} = PerCondition(k).file;                                 % name of the file
        Condition{idx,1} = PerCondition(k).condition;                       % label of the experimental condition
        dose{idx,1} = PerCondition(k).dose;                                 % dose/concentration of the treatment drug

        NetworkBurstFeature = PerCondition(k).NetworkBurstFeature0;         % all the events (no clustering)
        NetworkBurstStructure = PerCondition(k).NetworkBurstStructure;
        NetworkBurstDetected = PerCondition(k).NetworkBurstDetected;
        Number_NETburst_Events(idx) = NetworkBurstDetected.nBurst;          % total number of network burst events

        nNetBurst1 = NetworkBurstDetected.nBurst;
        nNetBurst2 = NetworkBurstDetected.Secondary.nBurst;

        if nNetBurst1/(PerCondition(k).window(2)-PerCondition(k).window(1)) > MIN_NETBURST            % if there is not enough number of network burst events, the population metrics stay NaN

            Number_NETburst_Neurons(idx) = NetworkBurstFeature.nNeuron.Burst;             % total number of clusters participating in network bursting
            NumberNorm_NETburst_Neurons(idx) = NetworkBurstFeature.nNeuron.normBurst * NetworkBurstFeature.nNeuron.normRec;          % total number of clusters participating in network bursting, normalized to the total number of clusters
            Number_Spikes_NETburst(idx) = NetworkBurstFeature.nSpikes.mean;               % total number of spikes per network burst, per cluster
            NumberNorm_Spikes_NETburst(idx) = NetworkBurstFeature.nSpikes.norm;           % total number of spikes in network burst, normalized to the total number of spikes in bursting participating clusters
            NumberNorm_SpikesTotal_NETburst(idx) = NetworkBurstFeature.nSpikes.normTotal; % total number of spikes in network burst, normalized to the total number of spikes across all clusters
            Participation_NETburst(idx) = NetworkBurstFeature.neuronsAct.mean;            % mean number of active neurons within burst normalized to the total number of globally participating neurons
            ParticipationCV_NETburst(idx) = NetworkBurstFeature.neuronsAct.cv;
            Rate_NETburst(idx) = NetworkBurstFeature.rate.mean;                           % rate of network bursting, in 1/min
            RateCV_NETburst(idx) = NetworkBurstFeature.rate.cv;                           % coefficient of variation for rate of network bursting 
            Duration_NETburst(idx) = NetworkBurstFeature.duration.mean;                   % duration of network burst, in sec 
            DurationCV_NETburst(idx) = NetworkBurstFeature.duration.cv;                   % coefficient of variation for duration of network burst 
            Amplitude_NETburst(idx) = NetworkBurstFeature.amplitude.mean;                 % amplitude (spiking frequency) per cluster in network burst, in Hz
            AmplitudeCV_NETburst(idx) = NetworkBurstFeature.amplitude.cv;                 % coefficient of variation of burst amplitude
            Distance_NETburst(idx) = NetworkBurstFeature.distance.mean;                   % distance of the network burst to the individual spike trains, in msec (peak to midpoint) 
            DistanceCV_NETburst(idx) = NetworkBurstFeature.distance.cv;
            Onset_NETburst(idx) = mean(NetworkBurstStructure.structure(:,2),'omitnan');   % onset of the activity, in ms (relative to the peak of the burst)  
            OnsetCV_NETburst(idx) = std(NetworkBurstStructure.structure(:,2),'omitnan');  % standard deviation, in ms  
            Peak_NETburst(idx) = mean(NetworkBurstStructure.structure(:,3),'omitnan');    % peak of the activity, in ms (relative to the peak of the burst)   
            PeakCV_NETburst(idx) = std(NetworkBurstStructure.structure(:,3),'omitnan');   % standard deviation, in ms  
            Correlation_NETburst(idx) = NetworkBurstStructure.correlation;                % correlation coefficient  
            
            tempStructNeuron_NETburst(idx) = Number_NETburst_Neurons(idx) - sum(isnan(NetworkBurstStructure.structure(:,1))); % number of neurons with robust temporal structure

            popFire = nan(NetworkBurstFeature.nBurst,1);
            bin = (NetworkBurstStructure.populationFiring{1}(1,3)-NetworkBurstStructure.populationFiring{1}(1,2))*K_MSEC_SEC;
            for j = 1:NetworkBurstFeature.nBurst
                popFire(j,1) = max(NetworkBurstStructure.populationFiring{j}(:,1))/bin;
            end
            PeakFire_NETburst(idx) = mean(popFire(:,1),'omitnan')/Number_NETburst_Neurons(idx);
            PeakFireNorm_NETburst(idx) = PeakFire_NETburst(idx)/Amplitude_NETburst(idx);

            cluster1 = NetworkBurstDetected.clustering.label(1:nNetBurst1,1);
            aux_cluster = unique(cluster1);
            nFamilies = numel(aux_cluster);
            Number_clusters_Netburst(idx) = nFamilies;                      % number of families found from the primary peaks

            if nNetBurst2/nNetBurst1 > MIN_SECONDARY_RATE

               Feature_secondary = PerCondition(k).NetworkBurstFeature{Number_clusters_Netburst(idx)+1};         % all the events (no clustering)
           
                if Feature_secondary.PrimarywithSecondary > MIN_PARENT
                    
                    NumberNorm_SecondaryBurst(idx) = Feature_secondary.SecondaryperPrimary;         % number of secondary events, per primary peaks
                    ParentNorm_SecondaryBurst(idx) = Feature_secondary.PrimarywithSecondary;        % proportion of primary peaks with secondary events
                    Number_SecondaryBurst_Neurons(idx) = Feature_secondary.nNeuron.Burst;           % # of participating units in secondary events
                    NumberNorm_SecondaryBurst_Neurons(idx) = Feature_secondary.nNeuron.normRec;     % # of participating units, normalized to the total number of units in the recording
                    Amplitude_SecondaryBurst(idx) = Feature_secondary.amplitude.mean;               % mean amplitude of the secondary events
                    
                else
                    ParentNorm_SecondaryBurst(idx) = 0;
                    NumberNorm_SecondaryBurst(idx) = 0;
                end
            else
                ParentNorm_SecondaryBurst(idx) = 0;
                NumberNorm_SecondaryBurst(idx) = 0;
            end

        end

        clear NetworkBurstFeature NetworkBurstDetected Feature_secondary
    end
    clear PerCondition n_conditions
end

NETbursting = table(Organoid, File, ttx_qc, Condition, dose, Number_NETburst_Events, Number_clusters_Netburst, Number_NETburst_Neurons, NumberNorm_NETburst_Neurons, Number_Spikes_NETburst, NumberNorm_Spikes_NETburst, NumberNorm_SpikesTotal_NETburst, ...   ...
                    Participation_NETburst, ParticipationCV_NETburst, Rate_NETburst, RateCV_NETburst, Duration_NETburst, DurationCV_NETburst, Amplitude_NETburst, AmplitudeCV_NETburst, ...
                    Distance_NETburst, DistanceCV_NETburst, tempStructNeuron_NETburst, Onset_NETburst, OnsetCV_NETburst, Peak_NETburst, PeakCV_NETburst, Correlation_NETburst, ...
                    PeakFire_NETburst, PeakFireNorm_NETburst, Number_SecondaryBurst_Neurons, NumberNorm_SecondaryBurst_Neurons,NumberNorm_SecondaryBurst,ParentNorm_SecondaryBurst,Amplitude_SecondaryBurst);
end

% (revised October 2026) the network bursting results are compiled in the loop above ("Compile_NetBurstResults", options 0, 1 and 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the compiled results 

file_output = fullfile(folderPath, sprintf('%s_Compiled.mat', group));
save(file_output, 'Firing','Bursting','NETbursting0', 'NETbursting1','NETbursting2');

file_output_excel = fullfile(folderPath, sprintf('%s_Compiled.xlsx', group));

% units of each column, written in the first row of each sheet (revised October 2026: defined once; the rows of "Bursting" -one extra empty cell- 
% and of the network bursting tables -"#" on the 7th to 9th columns instead of the 6th, 7th and 9th- were shifted; "OnsetCV" and "PeakCV" are 
% standard deviations, in msec; revised October 2026: column "Clustering_Quality" after "Number_clusters_Netburst")
units_firing = {'','','','','','sec', 'sec', '#','#','#','Hz','ad','ad','msec','ad','#/window','ad','ad','ad','ad'};    % (revised October 2026: CV and Fano factors)
units_bursting = {'','','','','','norm', 'norm', '1/min','ad','ad','ad','sec','ad','#/window','sec','ad','1/sec','ad','ad','ad','msec','norm','1/sec','ad','ad','ad','msec','norm','1/sec','ad','ad','ad','msec'};
units_netbursting = {'','','','','','#','#','','','#','norm','#','norm','norm','norm','ad','1/min','ad','ad','sec','ad','#/window','sec','ad','ad','ad','1/sec','ad','msec','ad','#','msec','msec','msec','msec','ad','Hz','ad','#','norm','#/primary','ad','Hz'};

writetable(Firing, file_output_excel,'Sheet','Firing','Range','A2',"WriteMode",'replacefile')
writecell(units_firing,file_output_excel,'Sheet','Firing','Range','A1')

writetable(Bursting, file_output_excel,'Sheet','Bursting','Range','A2')
writecell(units_bursting,file_output_excel,'Sheet','Bursting','Range','A1')

writetable(NETbursting0, file_output_excel,'Sheet','Network Bursting Pooled','Range','A2')
writecell(units_netbursting,file_output_excel,'Sheet','Network Bursting Pooled','Range','A1')

writetable(NETbursting1, file_output_excel,'Sheet','Network Bursting Primary','Range','A2')
writecell(units_netbursting,file_output_excel,'Sheet','Network Bursting Primary','Range','A1')

writetable(NETbursting2, file_output_excel,'Sheet','Network Bursting ALL clusters','Range','A2')
writecell(units_netbursting,file_output_excel,'Sheet','Network Bursting ALL clusters','Range','A1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Group the files per condition and Save the results 

file_output_excel = fullfile(folderPath, sprintf('%s_Grouped&Compiled.xlsx', group));

% (revised October 2026) warning for the conditions that do not match any treatment label (they are missing from the grouped file) or that 
% match more than one (they appear more than once)
condition_all = lower(string(Firing.(COL_COND)));
n_match = zeros(height(Firing), 1);
for k=1:n_treatment
    for i=1:height(Firing)
        n_match(i) = n_match(i) + (stringSimilarity(condition_all(i), lower(TREATMENT{k})) >= RATIO_MATCH);
    end
end
if any(n_match == 0)
    fprintf("\n WARNING: conditions without a matching treatment label (not included in the grouped file): %s \n", strjoin(unique(condition_all(n_match == 0)), ', '));
end
if any(n_match > 1)
    fprintf("\n WARNING: conditions matching more than one treatment label (included more than once in the grouped file): %s \n", strjoin(unique(condition_all(n_match > 1)), ', '));
end

for k=1:n_treatment
    condition = lower(string(Firing.(COL_COND)));                                  % case-insensitive: convert to lowercase
    target = lower(TREATMENT{k});
    n = height(Firing);
    simScore = zeros(n, 1);
    for i = 1:n
        simScore(i) = stringSimilarity(condition(i), target);
    end
    group_cond = Firing(simScore >= RATIO_MATCH, :);

    if k==1
        writetable(group_cond, file_output_excel,'Sheet','Firing','Range','A2',"WriteMode",'replacefile')
%        writecell({'', '','sec', 'sec', '#','#','#','Hz','ad','ad','ad','ad','ad'},file_output_excel,'Sheet','Firing','Range','A1')
        writecell(units_firing,file_output_excel,'Sheet','Firing','Range','A1')
    else
        writetable(group_cond, file_output_excel,'Sheet','Firing',"WriteMode",'append')
    end

    clear simScore target group_cond
end

for k=1:n_treatment
    condition = lower(string(Bursting.(COL_COND)));                                  % case-insensitive: convert to lowercase
    target = lower(TREATMENT{k});
    n = height(Bursting);
    simScore = zeros(n, 1);
    for i = 1:n
        simScore(i) = stringSimilarity(condition(i), target);
    end
    group_cond = Bursting(simScore >= RATIO_MATCH, :);

    if k==1
        writetable(group_cond, file_output_excel,'Sheet','Bursting','Range','A2')
%        writecell({'', '','norm', 'norm', '1/min','ad','ad','sec','ad','1/sec','ad','ad','norm','1/sec','ad','ad','norm','1/sec','ad','ad'},file_output_excel,'Sheet','Bursting','Range','A1')
        writecell(units_bursting,file_output_excel,'Sheet','Bursting','Range','A1')
    else
        writetable(group_cond, file_output_excel,'Sheet','Bursting',"WriteMode",'append')
    end

    clear simScore target group_cond
end

for k=1:n_treatment
    condition = lower(string(NETbursting0.(COL_COND)));                                  % case-insensitive: convert to lowercase
    target = lower(TREATMENT{k});
    n = height(NETbursting0);
    simScore = zeros(n, 1);
    for i = 1:n
        simScore(i) = stringSimilarity(condition(i), target);
    end
    group_cond = NETbursting0(simScore >= RATIO_MATCH, :);

    if k==1
        writetable(group_cond, file_output_excel,'Sheet','Network Bursting Pooled','Range','A2')
%        writecell({'','', '#','#','norm','#','norm','norm','norm','ad','1/min','ad','sec','ad','1/sec','ad','msec','ad','#','msec','ad','msec','ad','ad','Hz','ad'},file_output_excel,'Sheet','Network Bursting','Range','A1')
        writecell(units_netbursting,file_output_excel,'Sheet','Network Bursting Pooled','Range','A1')
    else
        writetable(group_cond, file_output_excel,'Sheet','Network Bursting Pooled',"WriteMode",'append')
    end

    clear simScore target group_cond
end

end


