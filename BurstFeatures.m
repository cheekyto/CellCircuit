function BurstFeatures=BurstFeatures(BurstDetection, tRange)

% function to analyze the bursting activity (simple spiking is analyzed by a different function) by calculating different metrics on MEA-based recordings
% (revised October 2026: the number of spikes per burst is counted in a single pass, instead of searching the whole spike train once per burst;
%  new variability metrics: Fano factor of the isi / ibi in two forms -CV^2 and variance/mean- and Fano factor of the burst counts in windows
%  of COUNT_WINDOW sec -see "countFanoFactor"-; the previous "fanoFactor" -CV^2, and 2^CV before the review- was removed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% structure array BurstDetection.IndivBurst (its length is the number of burst), with the following fields
%   .start: time stamp of the first spike of the burst, in ms
%   .last:  time stamp of the last spike of the burst, in ms
%   .peak:  time stamp of the peak the burst, in ms
%   .isi(i):isi corresponding to the spike i of the burst
%
% BurstDetection.spikes(i,1):   time stamp of spike i, in ms
%               .spikes(i,2):   isi corresponding to the spike i, in ms     [isi(j)=T(j+1)-T(j)]     
%               .spikes(i,3):   # of burst which the spike i belongs to (=0 if spike i is NOT part of a burst)
% tRange:       (optional, revised October 2026) [start end] of the analysis window, in sec, to count the bursts in windows (Fano factor of 
%               the counts); without it, the count metrics are NaN (e.g., when it is called to shuffle the spike trains)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BurstFeatures.nBurst                                      % total number of burst 
% BurstFeatures.indivMetric(i).nSpikes                      % total number of spikes in burst i
% BurstFeatures.indivMetric(i).duration                     % duration of the burst i, in sec
% BurstFeatures.indivMetric(i).isi.mean                     % mean isi for the burst i, in ms
% BurstFeatures.indivMetric(i).isi.median                   % median isi for the burst i, in ms
% BurstFeatures.indivMetric(i).isi.cv                       % coefficient of variation for the isi in burst i
% BurstFeatures.indivMetric(i).isi.FF_cv2                   % Fano factor for the isi in burst i, as CV^2
% BurstFeatures.indivMetric(i).isi.FF_varmean               % Fano factor for the isi in burst i, as variance / mean, in ms

% BurstFeatures.neuronMetric.nSpikes.mean                   % mean number of spikes per burst
% BurstFeatures.neuronMetric.nSpikes.median                 % median number of spikes per burst
% BurstFeatures.neuronMetric.nSpikes.cv                     % coefficient of variation for the number of spikes per burst
% BurstFeatures.neuronMetric.nSpikes.norm                   % total number of spikes within burst normalized by the total number of spikes in the recording (ratio of spikes participating in bursting)
% BurstFeatures.neuronMetric.duration.mean                  % mean burst duration, in sec
% BurstFeatures.neuronMetric.duration.median                % median (across burst events) burst duration, in sec
% BurstFeatures.neuronMetric.duration.cv                    % coefficient of variation for the burst duration
% BurstFeatures.neuronMetric.isi.mean                       % mean (across burst events) of the mean isi within bursts, in ms
% BurstFeatures.neuronMetric.isi.median                     % median (across burst events) of the median isi within bursts, in ms
% BurstFeatures.neuronMetric.isi.cv                         % mean coefficient of variation for the isi
% BurstFeatures.neuronMetric.isi.FF_cv2                     % mean (across burst events) Fano factor for the isi, as CV^2
% BurstFeatures.neuronMetric.isi.FF_varmean                 % mean (across burst events) Fano factor for the isi, as variance / mean, in ms
% BurstFeatures.neuronMetric.ibi.mean                       % mean inter-burst interval, in sec
% BurstFeatures.neuronMetric.ibi.median                     % median inter-burst interval, in sec
% BurstFeatures.neuronMetric.ibi.cv                         % coefficient of variation for the IBI
% BurstFeatures.neuronMetric.ibi.FF_cv2                     % Fano factor for the IBI, as CV^2
% BurstFeatures.neuronMetric.ibi.FF_varmean                 % Fano factor for the IBI, as variance / mean, in sec
% BurstFeatures.neuronMetric.count.FFcount                  % Fano factor of the number of bursts in windows of COUNT_WINDOW sec (NaN without tRange)
% BurstFeatures.neuronMetric.count.meanCount                % mean number of bursts per window
% BurstFeatures.neuronMetric.count.nWindows                 % number of windows

% BurstFeatures.neuronMetric.nonBursting.isi.mean           % mean inter-burst interval outside of bursts, in ms
% BurstFeatures.neuronMetric.nonBursting.isi.median         % median inter-burst interval outside of bursts, in ms
% BurstFeatures.neuronMetric.nonBursting.isi.cv             % coefficient of variation for the ISI outside of bursts
% BurstFeatures.neuronMetric.nonBursting.isi.FF_cv2         % Fano factor for the ISI outside of bursts, as CV^2
% BurstFeatures.neuronMetric.nonBursting.isi.FF_varmean     % Fano factor for the ISI outside of bursts, as variance / mean, in ms
% BurstFeatures.neuronMetric.nonBursting.nSpikes_norm       % proportion of non-bursting spikes (ratio of spikes non-participating in bursting)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                                                        % factor to rescale from msec to sec
COUNT_WINDOW = 60;                                                          % length of the windows to count the bursts for the Fano factor of the counts, in sec (default: 60)
apply_parameters('BurstFeatures');                                          % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")
if nargin < 2
    tRange = [];
end

n_burst = BurstDetection.nBurst;                                            % total number of burst events
n_spike = size(BurstDetection.spikes,1);                                    % total number of spikes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyze the non-bursting spikes (spikes outside of bursts)

burst_label = BurstDetection.spikes(:,3);                                                                           % burst to which each spike belongs (=0, outside of bursts)
aux = burst_label==0;                                                                                               % identify all the spikes that do not belong to a burst
n_spikes_noBurst = sum(aux);                                                                                        % total number of spikes outside of bursts
isi_noBurst = BurstDetection.spikes(aux,2);                                                                         % isi for spikes outside of bursts

BurstFeatures.neuronMetric.nonBursting.nSpikes_norm = n_spikes_noBurst/n_spike;                                     % normalized number of spikes outside of bursts
if n_spikes_noBurst>0
    BurstFeatures.neuronMetric.nonBursting.isi.mean = mean(isi_noBurst,'omitnan');                                  % mean isi for spikes outside of bursts
    BurstFeatures.neuronMetric.nonBursting.isi.median = median(isi_noBurst,'omitnan');                              % median isi for spikes outside of bursts
    BurstFeatures.neuronMetric.nonBursting.isi.cv = std(isi_noBurst,'omitnan')/mean(isi_noBurst,'omitnan');         % coefficient of variation for isi of spikes outside of bursts
    BurstFeatures.neuronMetric.nonBursting.isi.FF_cv2 = BurstFeatures.neuronMetric.nonBursting.isi.cv^2;            % Fano factor for isi of spikes outside of bursts, as CV^2
    BurstFeatures.neuronMetric.nonBursting.isi.FF_varmean = var(isi_noBurst,'omitnan')/mean(isi_noBurst,'omitnan'); % Fano factor for isi of spikes outside of bursts, as variance / mean, in ms
else                        
    BurstFeatures.neuronMetric.nonBursting.isi.mean = nan;                                                          % if there are no spikes outside of bursts, nan all metrics
    BurstFeatures.neuronMetric.nonBursting.isi.median = nan;
    BurstFeatures.neuronMetric.nonBursting.isi.cv = nan;
    BurstFeatures.neuronMetric.nonBursting.isi.FF_cv2 = nan;
    BurstFeatures.neuronMetric.nonBursting.isi.FF_varmean = nan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyze the inter-burst interval (measured from start-to-start of each burst event)

BurstFeatures.nBurst = n_burst;                                                         % total number of burst events
aux = nan(n_burst,1);

for i=1:n_burst                                                                         % across all burst events
    aux(i) = BurstDetection.IndivBurst(i).start;                                        % time stamp of the first spike in burst i, in msec
end
ibi = diff(aux)*K_MSEC_SEC;                                                             % inter-burst interval, in sec (even if there are other spikes in between bursts)
BurstFeatures.neuronMetric.ibi.mean = mean(ibi, 'omitnan');                             % mean inter-burst interval, in sec
BurstFeatures.neuronMetric.ibi.median = median(ibi, 'omitnan');                         % median inter-burst interval, in sec
BurstFeatures.neuronMetric.ibi.cv = std(ibi,'omitnan')/mean(ibi, 'omitnan');            % coefficient of variation for the IBI
BurstFeatures.neuronMetric.ibi.FF_cv2 = BurstFeatures.neuronMetric.ibi.cv^2;            % Fano factor for the IBI, as CV^2
BurstFeatures.neuronMetric.ibi.FF_varmean = var(ibi,'omitnan')/mean(ibi, 'omitnan');   % Fano factor for the IBI, as variance / mean, in sec
if numel(tRange) == 2                                                                   % bursts counted in windows of COUNT_WINDOW sec (by the time of their first spike)
    counting = countFanoFactor(aux, tRange, COUNT_WINDOW);
else
    counting = struct('FFcount', nan, 'meanCount', nan, 'nWindows', nan);
end
BurstFeatures.neuronMetric.count.FFcount = counting.FFcount;                            % Fano factor of the burst counts
BurstFeatures.neuronMetric.count.meanCount = counting.meanCount;                        % mean number of bursts per window
BurstFeatures.neuronMetric.count.nWindows = counting.nWindows;                          % number of windows

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyze individual bursts

duration = nan(n_burst,1); isi_mean = nan(n_burst,1); isi_median = nan(n_burst,1); isi_cv = nan(n_burst,1); isi_FF = nan(n_burst,1); isi_FFvm = nan(n_burst,1);
aux = burst_label(burst_label>0);                                                                       % burst number of all the spikes within bursts
nspikes_burst = accumarray(aux(:), 1, [n_burst 1]);                                                     % total number of spikes within each burst (counted in a single pass across all spikes)
if n_burst > 0                                                                                          % pre-allocate the structure array (same fields and order as before)
    BurstFeatures.indivMetric = struct('duration', cell(1,n_burst), 'isi', [], 'nSpikes', []);
end

for i=1:n_burst                                                                                         % across all burst events
    duration(i) = (BurstDetection.IndivBurst(i).last-BurstDetection.IndivBurst(i).start)*K_MSEC_SEC;    % duration of the burst i, in sec
    isi_mean(i) = mean(BurstDetection.IndivBurst(i).isi,'omitnan');                                     % mean isi for the burst i, in msec
    isi_median(i) = median(BurstDetection.IndivBurst(i).isi,'omitnan');                                 % median isi for the burst i, in msec
    isi_cv(i) = std(BurstDetection.IndivBurst(i).isi,'omitnan')/isi_mean(i);                            % coefficient of variation for the isi
    isi_FF(i) = isi_cv(i)^2;                                                                            % Fano factor for the isi, as CV^2
    isi_FFvm(i) = var(BurstDetection.IndivBurst(i).isi,'omitnan')/isi_mean(i);                          % Fano factor for the isi, as variance / mean, in ms

    BurstFeatures.indivMetric(i).duration = duration(i);                                                % duration of the burst i, in sec
    BurstFeatures.indivMetric(i).isi.mean = isi_mean(i);                                                % mean isi for the burst i, in ms
    BurstFeatures.indivMetric(i).isi.median = isi_median(i);                                            % median isi for the burst i, in ms
    BurstFeatures.indivMetric(i).isi.cv = isi_cv(i);                                                    % coefficient of variation for the isi in burst i
    BurstFeatures.indivMetric(i).isi.FF_cv2 = isi_FF(i);                                                % Fano factor for the isi in burst i, as CV^2
    BurstFeatures.indivMetric(i).isi.FF_varmean = isi_FFvm(i);                                          % Fano factor for the isi in burst i, as variance / mean, in ms
    BurstFeatures.indivMetric(i).nSpikes = nspikes_burst(i);                                            % total number of spikes in burst i
end

BurstFeatures.neuronMetric.nSpikes.mean = mean(nspikes_burst,'omitnan');                                % mean number of spikes per burst
BurstFeatures.neuronMetric.nSpikes.median = median(nspikes_burst,'omitnan');                            % median number of spikes per burst
BurstFeatures.neuronMetric.nSpikes.cv = std(nspikes_burst,'omitnan')/mean(nspikes_burst,'omitnan');     % coefficient of variation for the number of spikes per burst
BurstFeatures.neuronMetric.nSpikes.norm = sum(nspikes_burst)/n_spike;                                   % total number of spikes within burst normalized by the total number of spikes in the recording (ratio of spikes participating in bursting)
BurstFeatures.neuronMetric.duration.mean = mean(duration,'omitnan');                                    % mean (across burst events) burst duration, in sec
BurstFeatures.neuronMetric.duration.median = median(duration,'omitnan');                                % median (across burst events) burst duration, in sec
BurstFeatures.neuronMetric.duration.cv = std(duration,'omitnan')/mean(duration,'omitnan');              % coefficient of variation (across burst events) for the burst duration 
BurstFeatures.neuronMetric.isi.mean = mean(isi_mean,'omitnan');                                         % mean (across burst events) of the mean isi within bursts, in ms
BurstFeatures.neuronMetric.isi.median = median(isi_median,'omitnan');                                   % median (across burst events) of the median isi within bursts, in ms
BurstFeatures.neuronMetric.isi.cv = mean(isi_cv,'omitnan');                                             % mean (across burst events) coefficient of variation for the isi
BurstFeatures.neuronMetric.isi.FF_cv2 = mean(isi_FF,'omitnan');                                         % mean (across burst events) Fano factor for the isi, as CV^2
BurstFeatures.neuronMetric.isi.FF_varmean = mean(isi_FFvm,'omitnan');                                   % mean (across burst events) Fano factor for the isi, as variance / mean, in ms

end