function SpikeFeatures = SpikeFeatures(spikeTimestamp, tRange)

% function to analyze the spiking activity (bursting is analyzed by a different function) by calculating 22 different metrics on MEA-based recordings
% It is called by the function "Analyze_spikeTrain"
% (revised October 2026: the STTC uses the analysis window -tRange- as the recording interval; faster entropy, mutual information and STTC 
%  -see the functions "Entropy", "Mutual_Information" and "sttc_one_direction"-; new variability metrics: CV of the isi, Fano factor of the isi 
%  in two forms -CV^2 and variance/mean- and Fano factor of the spike counts in windows of COUNT_WINDOW sec -see "countFanoFactor"-; 
%  the previous "fanoFactor" -CV^2, and 2^CV before the review- was removed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   structure array spikeTimestamp (its length is the number of recording channels), with the following fields
% 	.recLength:	length of the recording, in seconds.
% 	.nUnits:	total number of identified units in the channel
% 	.timeStamp:	time stamp for all units (nUnits) in the channel, in msec
%   .window:    [start end] time window for the unit extraction, in sec
% 	.position:	position of the electrode
% 	.cluster:	cluster id
%   .t_step:    temporal resolution, in msec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SpikeFeatures.nSpikes(i)                     % total number of spikes in electrode i                    
% SpikeFeatures.nSpikes_norm(i)                % normalized number of spikes in electrode i 
% SpikeFeatures.entropy(i)                     % entropy in electrode i
% SpikeFeatures.lvr(i)                         % local variation ratio in electrode i
% SpikeFeatures.mutualInfo(i,j):               % mutual information for the pair of electrodes (i,j) 
% SpikeFeatures.sttc(i,j):                     % spike train tiling coefficient for the pair of electrodes (i,j)
% SpikeFeatures.isi(i).time(j)                 % time stamp for the spike j in electrode i, in msec
%                     .value(j)                % inter-spike interval associated to the spike j in electrode i, in msec [isi(j)=T(j+1)-T(j)]
% SpikeFeatures.isiMetric(i).mean              % mean isi for electrode i, in msec  
%                           .median            % median isi for electrode i, in msec
%                           .std               % isi standard deviation for electrode i
%                           .cv                % isi coefficient of variation for electrode i
%                           .FF_cv2            % isi Fano factor, as CV^2 (= the Fano factor of the counts for long windows, for a renewal process)
%                           .FF_varmean        % isi Fano factor, as variance / mean of the isi, in msec
% SpikeFeatures.countMetric(i).FFcount         % Fano factor of the number of spikes in windows of COUNT_WINDOW sec (variance / mean across windows)
%                            .meanCount        % mean number of spikes per window
%                            .nWindows         % number of windows
% SpikeFeatures.firingRate(i).isi              % firing rate, based on the isi, for electrode i, in 1/sec
%                            .counting         % firing rate, based on counting the total # of spikes, for electrode i
% SpikeFeatures.Population.firingRateISI       % population mean firing rate, based on the isi, in 1/sec 
%                         .firingRateCount     % population mean firing rate, based on counting the total # of spikes, in 1/sec 
%                         .cv                  % population mean isi coefficient of variation
%                         .FFisi_cv2           % population mean isi Fano factor (CV^2)
%                         .FFisi_varmean       % population mean isi Fano factor (variance / mean), in msec
%                         .FFcount             % population mean Fano factor of the spike counts
%                         .meanCount           % population mean number of spikes per counting window
%                         .entropy             % population mean entropy 
%                         .lvr                 % population mean lvr
%                         .mutualInfo          % population mean Mutual Information
%                         .sttc                % population mean STTC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_BIN_ENTROPY = 30;                            % temporal bin to calculate the entropy and the mutual information, in msec (default: 30)
DT_STTC = 15;                                   % dt for the spike train tiling coefficient, in msec (default: 15)
COUNT_WINDOW = 60;                              % length of the windows to count the spikes for the Fano factor of the counts, in sec (default: 60)
apply_parameters('SpikeFeatures');                                          % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

K_MSEC_SEC = 1/1000;                            % factor to rescale from msec to sec
n_electr = size(spikeTimestamp,2);              % total number of clusters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inter-spike interval and number of spikes per electrode

nSpikes = zeros(n_electr,1);    totalSpikes = 0;
cv_isi = nan(n_electr,1);   FF_cv2 = nan(n_electr,1);   FF_varmean = nan(n_electr,1);   FFcount = nan(n_electr,1);   meanCount = nan(n_electr,1);

for elect=1:n_electr                                                        % across all clusters
    nSpikes(elect) = spikeTimestamp(elect).nUnits;                          % total number of good spikes for cluster "elect"
    totalSpikes = totalSpikes + nSpikes(elect);                             % counting only the active electrodes and after removing the spikes that are too close
    isi_indv = diff(spikeTimestamp(elect).timeStamp);                       % inter-spike interval, in msec
    time_indiv = spikeTimestamp(elect).timeStamp(1:end-1);                  % start of each ISI, in msec 

    SpikeFeatures.isi(elect).time = time_indiv;                                             % time stamp for each isi (first spike of the pair), in msec
    SpikeFeatures.isi(elect).value = isi_indv;                                              % value of the isi, in msec
    SpikeFeatures.isiMetric(elect).mean = mean(isi_indv,'omitnan');                         % mean isi, in msec
    SpikeFeatures.isiMetric(elect).median = median(isi_indv,'omitnan');                     % median isi, in msec
    SpikeFeatures.isiMetric(elect).std = std(isi_indv,'omitnan');                           % isi standard deviation, in msec  
    SpikeFeatures.isiMetric(elect).cv = std(isi_indv,'omitnan')/mean(isi_indv,'omitnan');   % isi coefficient of variation
    SpikeFeatures.isiMetric(elect).FF_cv2 = SpikeFeatures.isiMetric(elect).cv^2;            % (revised October 2026) isi Fano factor, as CV^2
    SpikeFeatures.isiMetric(elect).FF_varmean = var(isi_indv,'omitnan')/mean(isi_indv,'omitnan');   % (revised October 2026) isi Fano factor, as variance / mean, in msec
    counting = countFanoFactor(spikeTimestamp(elect).timeStamp, tRange, COUNT_WINDOW);      % (revised October 2026) spikes counted in windows of COUNT_WINDOW sec
    SpikeFeatures.countMetric(elect).FFcount = counting.FFcount;                            % Fano factor of the spike counts
    SpikeFeatures.countMetric(elect).meanCount = counting.meanCount;                        % mean number of spikes per window
    SpikeFeatures.countMetric(elect).nWindows = counting.nWindows;                          % number of windows
    cv_isi(elect) = SpikeFeatures.isiMetric(elect).cv;                                      % (to calculate the population means)
    FF_cv2(elect) = SpikeFeatures.isiMetric(elect).FF_cv2;
    FF_varmean(elect) = SpikeFeatures.isiMetric(elect).FF_varmean;
    FFcount(elect) = counting.FFcount;
    meanCount(elect) = counting.meanCount;
end

SpikeFeatures.nSpikes = nSpikes;                                            % total number of good spikes for each cluster                                        
SpikeFeatures.nSpikes_norm = nSpikes / totalSpikes;                         % proportion of spikes coming from this cluster (it does not take in account potential differences in the lenght of recording between electrodes)
SpikeFeatures.Population.cv = mean(cv_isi,'omitnan');                       % (revised October 2026) mean isi coefficient of variation (across active clusters)
SpikeFeatures.Population.FFisi_cv2 = mean(FF_cv2,'omitnan');                % mean isi Fano factor, as CV^2 (across active clusters)
SpikeFeatures.Population.FFisi_varmean = mean(FF_varmean,'omitnan');        % mean isi Fano factor, as variance / mean, in msec
SpikeFeatures.Population.FFcount = mean(FFcount,'omitnan');                 % mean Fano factor of the spike counts
SpikeFeatures.Population.meanCount = mean(meanCount,'omitnan');             % mean number of spikes per counting window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Firing rate

firingRate_isi = nan(n_electr,1);   firingRate_counting  = nan(n_electr,1);

for elect=1:n_electr
    SpikeFeatures.firingRate(elect).isi = 1/(SpikeFeatures.isiMetric(elect).mean*K_MSEC_SEC);
    SpikeFeatures.firingRate(elect).counting = nSpikes(elect) / spikeTimestamp(elect).recLength;
    firingRate_isi(elect) = SpikeFeatures.firingRate(elect).isi;
    firingRate_counting(elect)  = SpikeFeatures.firingRate(elect).counting;
end

SpikeFeatures.Population.firingRateISI = mean(firingRate_isi,'omitnan');            % mean spiking frequency (calculated from isi), in 1/sec
SpikeFeatures.Population.firingRateCount = mean(firingRate_counting,'omitnan');     % mean spiking frequency (calculated from "counting" the spikes), in 1/sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Entropy

entropy = nan(n_electr,1);
signal = spikeTimestamp(1).timeStamp;                                   % time stamp time stamp for all units (nUnits) in the cluster, in msec    
aux = bin_spike_time(signal, N_BIN_ENTROPY, tRange);
n_bin = size(aux.binMidPoint,2);
spikeTrain_binned = nan(n_electr,n_bin);

for i=1:n_electr                                                            % across all clusters
    signal = spikeTimestamp(i).timeStamp;                                   % time stamp time stamp for all units (nUnits) in the cluster, in msec    
    aux = bin_spike_time(signal, N_BIN_ENTROPY, tRange);
    spikeTrain_binned(i,:) = aux.count;
    entropy(i) = Entropy(spikeTrain_binned(i,:));                          % entropy for cluster i    
end                                                                         % (revised October 2026: no "clear" inside the loop)

SpikeFeatures.entropy = entropy;                                            % entropy for cluster i
SpikeFeatures.Population.entropy = mean(entropy,'omitnan');                 % mean entropy for the population

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local variation ratio

localVarRatio = nan(n_electr,1);
for i=1:n_electr
    if spikeTimestamp(i).nUnits > 2
       isi = diff(spikeTimestamp(i).timeStamp);                             % isi for all spikes in cluster i (across all recording conditions)
       localVarRatio(i) = lvr(isi);
    end
end

SpikeFeatures.lvr = localVarRatio;                                          % local variation ratio for cluster i
SpikeFeatures.Population.lvr = mean(localVarRatio,'omitnan');               % mean lvr for the population

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mutual Information

mutualInfo = nan(n_electr, n_electr);
for i=1:n_electr
    for j=i+1:n_electr
        Mij = Mutual_Information(spikeTrain_binned(i,:), spikeTrain_binned(j,:));
        mutualInfo(i,j) = Mij;                                              % Since the Mutual information is symmetric
        mutualInfo(j,i) = Mij;
    end
end

SpikeFeatures.mutualInfo = mutualInfo;
SpikeFeatures.Population.mutualInfo = mean(mean(mutualInfo, 'omitnan'), 'omitnan');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spike Train Tiling Coefficient (STTC)

sttc_matrix = nan(n_electr, n_electr);
for i=1:n_electr
    for j=1+i:n_electr
        Cij = sttc(spikeTimestamp(i).timeStamp, spikeTimestamp(j).timeStamp, DT_STTC, tRange/K_MSEC_SEC);  % (revised October 2026) recording interval = analysis window, in msec
        sttc_matrix(i,j) = Cij;
        sttc_matrix(j,i) = Cij;
    end
end

SpikeFeatures.sttc = sttc_matrix;
SpikeFeatures.Population.sttc = mean(mean(abs(sttc_matrix), 'omitnan'), 'omitnan');

end