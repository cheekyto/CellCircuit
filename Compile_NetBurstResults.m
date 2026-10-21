function Compile_NetBurstResults = Compile_NetBurstResults(n_cluster, cond_idx, Recording, file_input, ttx_qc, option, thresholds)

%   n_cluster:      total number of clusters in the recording
%   Recording:      file to be assessed
%   file_input:     path of the "_Analyzed.mat" file, or (revised October 2026) the variable "PerCondition" of that file, already loaded
%   cond_idx:       index of the recording condition to be assessed
%   option:         =0, all primary clusters together
%                   =1, only the first primary cluster
%                   =2, all primary clusters, in individual rows
%                   ** all the different option include the secondary events
%   thresholds:     (optional, revised October 2026) structure with the fields MIN_NETBURST, MIN_SECONDARY_RATE and MIN_PARENT (set in "Compile_MEAproject");
%                   without it, the values below are used

MIN_NETBURST = 1/200;                      % minimum rate of network burst events, in # per sec (default: 1/150) (revised October 2026: 1/150 -> 1/200)
MIN_SECONDARY_RATE = 0.5;                  % minimum average number of secondary peaks per primary event (default: 0.5)
MIN_PARENT = 0.25;                         % minimum proportion of primary peaks with secondary peaks (default: 0.25)
if nargin >= 7 && ~isempty(thresholds)     % (revised October 2026) the values set in "Compile_MEAproject" (a single place for the thresholds)
    MIN_NETBURST = thresholds.MIN_NETBURST;
    MIN_SECONDARY_RATE = thresholds.MIN_SECONDARY_RATE;
    MIN_PARENT = thresholds.MIN_PARENT;
end
K_MSEC_SEC = 1/1000;                        % factor to rescale from msec to sec

switch option
    case 0  % all primary clusters together
        n_cluster = 1;
    case 1  % only the first primary cluster
        n_cluster = 1;
    case 2  % all primary clusters, in individual rows
end

%file_input = fullfile(folderPath, group, sprintf('%s_Analyzed.mat', Recording{i}));
if ischar(file_input) || isstring(file_input)                               % path of the "_Analyzed.mat" file
    load(file_input, 'PerCondition');
else                                                                        % (revised October 2026) the variable "PerCondition", already loaded by "Compile_MEAproject"
    PerCondition = file_input;                                              % (before, the file was loaded in every call: 3 times per condition)
end

NetworkBurstDetected = PerCondition(cond_idx).NetworkBurstDetected;
nNetBurst1 = NetworkBurstDetected.nBurst;
nNetBurst2 = NetworkBurstDetected.Secondary.nBurst;

% (revised October 2026) families of PRIMARY events: in "Analyze_spikeTrain", NetworkBurstFeature{k} is the family "families(k)" of the primary
% events (k = 1..n_families) and NetworkBurstFeature{n_families+1} are the secondary events. "clustering.n_cluster" also counts the families 
% with only secondary events, so it can be larger than n_families (before, it was used instead of n_families: the secondary events were 
% compiled as a primary family and the secondary features were read from a wrong index, which crashed)
families = 1;   n_families = 1;
if NetworkBurstDetected.nBurst > 0
    cluster1 = NetworkBurstDetected.clustering.label(1:nNetBurst1,1);       % family of each primary event
    families = unique(cluster1);                                            % families with primary events
    n_families = numel(families);
end

% Identify the primary cluster (cluster with the most number of events)
idx_cluster = 1;
if option == 1 && NetworkBurstDetected.nBurst > 0
    for i=2:n_families
        if PerCondition(cond_idx).NetworkBurstFeature{1,i}.nBurst > PerCondition(cond_idx).NetworkBurstFeature{1,idx_cluster}.nBurst
            idx_cluster = i;
        end
    end
end

%% Compile the results on network bursting 

Number_NETburst_Events = nan(n_cluster,1);                                % variables to temporaly collect the data from each recording
Number_NETburst_Neurons = nan(n_cluster,1);
NumberNorm_NETburst_Neurons = nan(n_cluster,1);
Number_Spikes_NETburst = nan(n_cluster,1);
NumberNorm_Spikes_NETburst = nan(n_cluster,1);
NumberNorm_SpikesTotal_NETburst = nan(n_cluster,1);
Participation_NETburst = nan(n_cluster,1);
ParticipationCV_NETburst = nan(n_cluster,1);
Rate_NETburst = nan(n_cluster,1);
RateCV_NETburst = nan(n_cluster,1);
FFibi_CV2_NETburst = nan(n_cluster,1);     FFibi_VarMean_NETburst = nan(n_cluster,1);     % (revised October 2026) Fano factors of the inter-burst intervals
FFcount_NETburst = nan(n_cluster,1);       MeanCount_NETburst = nan(n_cluster,1);         % (revised October 2026) Fano factor of the counts per window
Duration_NETburst = nan(n_cluster,1);
DurationCV_NETburst = nan(n_cluster,1);
Tail_NETburst = nan(n_cluster,1);
TailCV_NETburst = nan(n_cluster,1);
Amplitude_NETburst = nan(n_cluster,1);
AmplitudeCV_NETburst = nan(n_cluster,1);
Distance_NETburst = nan(n_cluster,1);
DistanceCV_NETburst = nan(n_cluster,1);
Onset_NETburst = nan(n_cluster,1);
OnsetCV_NETburst = nan(n_cluster,1);
Peak_NETburst = nan(n_cluster,1);
PeakCV_NETburst = nan(n_cluster,1);
Correlation_NETburst = nan(n_cluster,1);
tempStructNeuron_NETburst = nan(n_cluster,1);
PeakFire_NETburst = nan(n_cluster,1);
PeakFireNorm_NETburst = nan(n_cluster,1);
Number_clusters_Netburst = nan(n_cluster,1);                      % number of families found from the primary peaks
Clustering_Quality = nan(n_cluster,1);                            % (revised October 2026) quality of the clustering ("quality_clustering" in "cluster_neuronal_events"):
                                                                  % =2, reliable; =1, it seems to work (bootstrap / permutation not clear); =0, not significant; NaN, not run
label_cluster_Netburst = nan(n_cluster,1);          

Number_SecondaryBurst_Neurons = nan(n_cluster,1);
NumberNorm_SecondaryBurst_Neurons = nan(n_cluster,1);
Amplitude_SecondaryBurst = nan(n_cluster,1);
ParentNorm_SecondaryBurst = nan(n_cluster,1);
NumberNorm_SecondaryBurst = nan(n_cluster,1);

for k=1:n_cluster                                                           % across all network burst clusters    
    Organoid{k,1} = Recording;                                              % name of the experiment
    File{k,1} = PerCondition(cond_idx).file;                                % name of the file
    Condition{k,1} = PerCondition(cond_idx).condition;                      % label of the experimental condition
    dose{k,1} = PerCondition(cond_idx).dose;                                % dose/concentration of the treatment drug
        
    if NetworkBurstDetected.nBurst > 0 && isstruct(NetworkBurstDetected.clustering) && isfield(NetworkBurstDetected.clustering, 'clusteringQuality')
        Clustering_Quality(k) = NetworkBurstDetected.clustering.clusteringQuality;  % (revised October 2026) always reported, as the number of events (same for all the rows of the condition)
    end
    switch option
        case 0  % all primary clusters together
            NetworkBurstFeature = PerCondition(cond_idx).NetworkBurstFeature0;      % all the events (no clustering)
            Number_NETburst_Events(k) = NetworkBurstDetected.nBurst;                 % total number of network burst events
            NetworkBurstStructur = PerCondition(cond_idx).NetworkBurstStructure;
            label_cluster_Netburst(k) = nan;

        case 1  % only the first primary cluster
            NetworkBurstFeature = PerCondition(cond_idx).NetworkBurstFeature{1, idx_cluster}; % primary cluster (idx_cluster)
            Number_NETburst_Events(k) = NetworkBurstFeature.nBurst;                 % total number of network burst events
            label_cluster_Netburst(k) = families(idx_cluster);                     % (revised October 2026: label of the family, not its position)
            % prepare and run NetworkBurstStructure for the primary cluster
            if nNetBurst1 > 0
                netburst_member = (cluster1 == families(idx_cluster));                  % events corresponding to the cluster "idx_cluster"
                % (revised October 2026) the temporal structure of the family uses its own participating units (the same ones as its features, see
                % "familyParticipants"); before, all the participating units of primary and secondary events, so the units that do not participate in 
                % the family shifted its onset and peak, and "tempStructNeuron" and "PeakFire" mixed two different sets of units
                units_participant = familyParticipants(NetworkBurstDetected, PerCondition(cond_idx).BurstDetected, netburst_member, NetworkBurstDetected.Channels(:,1));
            else
                netburst_member = [];
                units_participant = [];
            end
            spikeTrain = PerCondition(cond_idx).spikeTrain;
            NetworkBurstStructur = NetworkBurstStructure(NetworkBurstDetected, spikeTrain, netburst_member, units_participant, 0);

        case 2  % all primary clusters, in individual rows
            NetworkBurstFeature = PerCondition(cond_idx).NetworkBurstFeature{1, k}; % information for cluster k
            Number_NETburst_Events(k) = NetworkBurstFeature.nBurst;                 % total number of network burst events
            label_cluster_Netburst(k) = families(k);                                % (revised October 2026: label of the family, not its position)

            % prepare and run NetworkBurstStructure for the cluster ksize
            if nNetBurst1 > 0
                netburst_member = (cluster1 == families(k));                  % events corresponding to the family k
                % (revised October 2026) the participating units of the family k, as for option 1 (see above)
                units_participant = familyParticipants(NetworkBurstDetected, PerCondition(cond_idx).BurstDetected, netburst_member, NetworkBurstDetected.Channels(:,1));
            else
                netburst_member = [];
                units_participant = [];
            end
            spikeTrain = PerCondition(cond_idx).spikeTrain;
            NetworkBurstStructur = NetworkBurstStructure(NetworkBurstDetected, spikeTrain, netburst_member, units_participant, 0);
            clear netburst_member
    end
            
    if Number_NETburst_Events(k)/(PerCondition(cond_idx).window(2)-PerCondition(cond_idx).window(1)) > MIN_NETBURST            % if there is not enough number of network burst events, the population metrics stay NaN
                                                                                                                                % (the number of events, "Number_NETburst_Events", is always reported)

            idx = k;
            Number_NETburst_Neurons(idx) = NetworkBurstFeature.nNeuron.Burst;             % total number of units participating in network bursting
            NumberNorm_NETburst_Neurons(idx) = NetworkBurstFeature.nNeuron.normBurst * NetworkBurstFeature.nNeuron.normRec;          % total number of units participating in network bursting, normalized to the total number of units
            Number_Spikes_NETburst(idx) = NetworkBurstFeature.nSpikes.mean;               % total number of spikes per network burst, per unit
            NumberNorm_Spikes_NETburst(idx) = NetworkBurstFeature.nSpikes.norm;           % total number of spikes in network burst, normalized to the total number of spikes in bursting participating units
            NumberNorm_SpikesTotal_NETburst(idx) = NetworkBurstFeature.nSpikes.normTotal; % total number of spikes in network burst, normalized to the total number of spikes across all units
            Participation_NETburst(idx) = NetworkBurstFeature.neuronsAct.mean;            % mean number of active neurons within burst normalized to the total number of globally participating neurons
            ParticipationCV_NETburst(idx) = NetworkBurstFeature.neuronsAct.cv;
            Rate_NETburst(idx) = NetworkBurstFeature.rate.mean;                           % rate of network bursting, in 1/min
            RateCV_NETburst(idx) = NetworkBurstFeature.rate.cv;                           % coefficient of variation for rate of network bursting 
            FFibi_CV2_NETburst(idx) = NetworkBurstFeature.rate.FF_cv2;                    % (revised October 2026) Fano factor of the inter-burst intervals, as CV^2
            FFibi_VarMean_NETburst(idx) = NetworkBurstFeature.rate.FF_varmean;            % Fano factor of the inter-burst intervals, as variance / mean, in sec
            FFcount_NETburst(idx) = NetworkBurstFeature.rate.FFcount;                     % Fano factor of the number of network bursts per window (COUNT_WINDOW sec, see "NetworkBurstFeatures")
            MeanCount_NETburst(idx) = NetworkBurstFeature.rate.meanCount;                 % mean number of network bursts per window
            Duration_NETburst(idx) = NetworkBurstFeature.duration.mean;                   % duration of network burst, in sec 
            DurationCV_NETburst(idx) = NetworkBurstFeature.duration.cv;                   % coefficient of variation for duration of network burst 
            Tail_NETburst(idx) = NetworkBurstFeature.tail.mean;
            TailCV_NETburst(idx) = NetworkBurstFeature.tail.cv;
            Amplitude_NETburst(idx) = NetworkBurstFeature.amplitude.mean;                 % amplitude (spiking frequency) per cluster in network burst, in Hz
            AmplitudeCV_NETburst(idx) = NetworkBurstFeature.amplitude.cv;                 % coefficient of variation of burst amplitude
            Distance_NETburst(idx) = NetworkBurstFeature.distance.mean;                   % distance of the network burst to the individual spike trains, in msec (peak to midpoint) 
            DistanceCV_NETburst(idx) = NetworkBurstFeature.distance.cv;
            if ~isempty(NetworkBurstStructur.structure)                                  % (revised October 2026) the metrics stay NaN if there is no structure (no participating units)
            Onset_NETburst(idx) = mean(NetworkBurstStructur.structure(:,2),'omitnan');   % onset of the activity, in ms (relative to the peak of the burst)  
            OnsetCV_NETburst(idx) = std(NetworkBurstStructur.structure(:,2),'omitnan');  % standard deviation, in ms (intended: the standard deviation, not the CV, because the mean is close to 0)
            Peak_NETburst(idx) = mean(NetworkBurstStructur.structure(:,3),'omitnan');    % peak of the activity, in ms (relative to the peak of the burst)   
            PeakCV_NETburst(idx) = std(NetworkBurstStructur.structure(:,3),'omitnan');   % standard deviation, in ms (intended, as above)
            Correlation_NETburst(idx) = NetworkBurstStructur.correlation;                % correlation coefficient  
            
            % (revised October 2026) the structure and the features use the same participating units (the family's units for options 1 and 2, the 
            % participating units of the primary events for option 0), so the number of units of the structure is "Number_NETburst_Neurons"
            tempStructNeuron_NETburst(idx) = Number_NETburst_Neurons(idx) - sum(isnan(NetworkBurstStructur.structure(:,1))); % number of neurons with robust temporal structure

            popFire = nan(Number_NETburst_Events(k),1);
            bin = (NetworkBurstStructur.populationFiring{1}(1,3)-NetworkBurstStructur.populationFiring{1}(1,2))*K_MSEC_SEC;
            for j = 1:Number_NETburst_Events(k)
                popFire(j,1) = max(NetworkBurstStructur.populationFiring{j}(:,1))/bin;
            end
            PeakFire_NETburst(idx) = mean(popFire(:,1),'omitnan')/Number_NETburst_Neurons(idx);
            PeakFireNorm_NETburst(idx) = PeakFire_NETburst(idx)/Amplitude_NETburst(idx);
            end

            Number_clusters_Netburst(idx) = n_families;                                 % number of families found from the primary peaks (revised October 2026: 
                                                                                        %  before, "clustering.n_cluster", which also counts the secondary-only families)

            %%% Assessment of secondary events  %%%
            if NetworkBurstDetected.nBurst_level1 > 0 && nNetBurst2/NetworkBurstDetected.nBurst_level1 > MIN_SECONDARY_RATE

               Feature_secondary = PerCondition(cond_idx).NetworkBurstFeature{n_families+1};    % all the secondary events (revised October 2026: after the n_families
                                                                                                % families of primary events)
           
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
        
        clear NetworkBurstFeature Feature_secondary NetworkBurstStructur 
end

NETbursting = table(Organoid, File, ttx_qc, Condition, dose, Number_NETburst_Events, Number_clusters_Netburst, Clustering_Quality, label_cluster_Netburst, Number_NETburst_Neurons, NumberNorm_NETburst_Neurons, Number_Spikes_NETburst, NumberNorm_Spikes_NETburst, NumberNorm_SpikesTotal_NETburst, ...
                    Participation_NETburst, ParticipationCV_NETburst, Rate_NETburst, RateCV_NETburst, FFibi_CV2_NETburst, FFibi_VarMean_NETburst, FFcount_NETburst, MeanCount_NETburst, Duration_NETburst, DurationCV_NETburst, Tail_NETburst, TailCV_NETburst, Amplitude_NETburst, AmplitudeCV_NETburst, ...
                    Distance_NETburst, DistanceCV_NETburst, tempStructNeuron_NETburst, Onset_NETburst, OnsetCV_NETburst, Peak_NETburst, PeakCV_NETburst, Correlation_NETburst, ...
                    PeakFire_NETburst, PeakFireNorm_NETburst, Number_SecondaryBurst_Neurons, NumberNorm_SecondaryBurst_Neurons,NumberNorm_SecondaryBurst,ParentNorm_SecondaryBurst,Amplitude_SecondaryBurst);

Compile_NetBurstResults = NETbursting;
end