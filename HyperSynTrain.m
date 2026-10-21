function QC = HyperSynTrain(signal1, signal2, t_step, range, n_iteration, method, refractory_period, BurstDetected2, BurstFeature2, ShuffleTrain)

% function to assess the statistical significance on hypersynchronicity - Analysis by comparing two channels
% It calls the functions "HyperSyncPair" and "ShuffleSpikeTrain"
% (revised October 2026: the shuffled spike trains can be provided as input -to reuse them across pairs-; and only the hypersynchronicity
%  level of each shuffling iteration is kept, not the full spike-level data of every iteration)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% signal1(i,1):         time stamp of the spike i recorded in channel 1, in msec (reference channel)
% signal1(i,2):         indicates if the spike i belongs to a burst (=0, No; =1, Yes)

% signal2(i,1):         time stamp of the spike i recorded in channel 1, in msec (study channel)
% signal2(i,2):         indicates if the spike i belongs to a burst (=0, No; =1, Yes)

% t_step:               temporal resolution - sampling rate, in ms
% n_iteration:          number of iteration for the shuffling
% method:               method for the shuffling of the spike train
% range:                range to generated the new spike train (range(1)=start; range(2)=end), in msec
% refractory_period:    absolute refractory period, in msec
% BurstDetected2:       output from the function "BurstDetection" corresponding to signal2
% BurstFeature2:        output from the function "BurstFeatures" corresponding to signal2
% ShuffleTrain:         (optional) output of "ShuffleSpikeTrain" for signal2, with at least "n_iteration" shuffled spike trains
%                       (the same for every pair with the same signal2, so it can be calculated only once); if it is not provided, it is calculated here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QC.burst.nspikes_hs                   proportion of spikes displaying hypersynchronicity within bursts
% QC.burst.percentile                   percentile/chance for the measured hypersynchronicity within burst (%) <the higher the percentile, the higher the chances to be statistically significant>
% QC.burst.nspikes_total                total number of spikes within bursts, ch1 + ch2
% QC.nonburst.nspikes_hs                proportion of spikes displaying hypersynchronicity outside of bursts
% QC.nonburst.percentile                percentile/chance for the measured hypersynchronicity outside of bursts (%)
% QC.nonburst.nspikes_total             total number of spikes within bursts, ch1 + ch2
% QC.total.nspikes_hs                   proportion of spikes displaying hypersynchronicity (entire recording)
% QC.total.percentile                   percentile/chance for the measured hypersynchronicity, ch1 + ch2 (entire recording) (%)
% QC.total.nspikes_total                total number of spikes within bursts, ch1 + ch2 (entire recording)
% QC.threshold                          threshold used to identify hypersynchronized spikes, in msec
% QC.HS_spikes_Indv{1}                  result for hypersynchronicity assessment done by "HyperSyncPair" for the original data set
%                                       (revised October 2026: the results for each shuffling iteration -HS_spikes_Indv{k+1}- are no longer saved, since they
%                                        contained two full copies of the spike trains per iteration; only their hypersynchronicity levels are kept, see below)
% QC.shuffled.burst(k)                  proportion of spikes displaying hypersynchronicity within bursts, for the shuffling iteration k
% QC.shuffled.nonburst(k)               proportion of spikes displaying hypersynchronicity outside of bursts, for the shuffling iteration k
% QC.shuffled.total(k)                  proportion of spikes displaying hypersynchronicity (entire recording), for the shuffling iteration k
% QC.shuffled.ch1(k) / ch2(k)           proportion of spikes displaying hypersynchronicity in channel 1 / 2, for the shuffling iteration k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Measure the level of hyper-synchronicity in the recorded data

HS_spikes0 = HyperSyncPair(signal1, signal2, t_step);                       % hypersynchronicity in the original dataset
QC.HS_spikes_Indv{1} = HS_spikes0;                                          % save the result for hypersynchronicity assessment corresponding to the original dataset
QC.burst.nspikes_hs = HS_spikes0.burst.nspikes_norm;                        % proportion of spikes displaying hypersynchronicity within bursts
QC.burst.nspikes_total = HS_spikes0.burst.nspikes_total;                    % total number of spikes within bursts, ch1 + ch2
QC.nonburst.nspikes_hs = HS_spikes0.nonburst.nspikes_norm;                  % proportion of spikes displaying hypersynchronicity outside of bursts
QC.nonburst.nspikes_total = HS_spikes0.nonburst.nspikes_total;              % total number of spikes outside of bursts, ch1 + ch2
QC.total.nspikes_hs = HS_spikes0.total.nspikes_norm;                        % proportion of spikes displaying hypersynchronicity (entire recording)
QC.total.nspikes_total = HS_spikes0.total.nspikes_total;                    % total number of spikes, ch1 + ch2 (entire recording)
QC.total.ch_1.nspikes_hs = HS_spikes0.total.ch_1.nspikes_norm;
QC.total.ch_1.nspikes_total = HS_spikes0.total.ch_1.nspikes_total;
QC.total.ch_2.nspikes_hs = HS_spikes0.total.ch_2.nspikes_norm;
QC.total.ch_2.nspikes_total = HS_spikes0.total.ch_2.nspikes_total;
QC.threshold = HS_spikes0.threshold;                                        % threshold used to identify hypersynchronized spikes, in msec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shuffle the spike trains and re-calculate the level of hypersynchronicity

if n_iteration>0                                                            % if iterations across shuffled data was requested

    if nargin < 10 || isempty(ShuffleTrain)                                 % the shuffled spike trains were not provided
        ShuffleTrain = ShuffleSpikeTrain(signal2(:,1), range, t_step, method, n_iteration, refractory_period,BurstDetected2,BurstFeature2);       % shuffle the signal2 "n_iteration" times
    end
    hs_burst = nan(n_iteration,1);    hs_nonburst = nan(n_iteration,1);    hs_total = nan(n_iteration,1);                                % initialize the variables to load the result of the iterations
    hs_ch1 = nan(n_iteration,1);    hs_ch2 = nan(n_iteration,1);

    for k=1:n_iteration                                                     % across all iterations
        HS_spikes = HyperSyncPair(signal1, ShuffleTrain{k}, t_step);        % hypersynchronicity analysis between signal1 and k-shuffled signal2
        hs_burst(k) = HS_spikes.burst.nspikes_norm;                         % proportion of spikes displaying hypersynchronicity within bursts, using suffling iteration k
        hs_nonburst(k) = HS_spikes.nonburst.nspikes_norm;                   % proportion of spikes displaying hypersynchronicity outside of bursts, using suffling iteration k
        hs_total(k) = HS_spikes.total.nspikes_norm;                         % proportion of spikes displaying hypersynchronicity (entire recording), using suffling iteration k
        hs_ch1(k) = HS_spikes.total.ch_1.nspikes_norm;
        hs_ch2(k) = HS_spikes.total.ch_2.nspikes_norm;
    end                                                                     % (revised October 2026: no "clear" inside the loop)
    QC.shuffled = struct('burst', hs_burst, 'nonburst', hs_nonburst, 'total', hs_total, 'ch1', hs_ch1, 'ch2', hs_ch2);    % hypersynchronicity level for each shuffling iteration

    QC.burst.percentile = numel(find(hs_burst<QC.burst.nspikes_hs))/n_iteration*100;            % percentile/chance for the measured hypersynchronicity within burst (%)
    QC.nonburst.percentile = numel(find(hs_nonburst<QC.nonburst.nspikes_hs))/n_iteration*100;   % percentile/chance for the measured hypersynchronicity outside of bursts (%)
    QC.total.percentile = numel(find(hs_total<QC.total.nspikes_hs))/n_iteration*100;            % percentile/chance for the measured hypersynchronicity, ch1 + ch2 (entire recording) (%)
    QC.ch1.percentile = numel(find(hs_ch1<QC.total.ch_1.nspikes_hs))/n_iteration*100;
    QC.ch2.percentile = numel(find(hs_ch2<QC.total.ch_2.nspikes_hs))/n_iteration*100;
%    QC.ch1.percentile = numel(find(hs_ch1<QC.total.ch_2.nspikes_norm))/n_iteration*100;
%    QC.ch2.percentile = numel(find(hs_ch2<QC.total.ch_2.nspikes_norm))/n_iteration*100;
else
    QC.burst.percentile = nan;                                              % percentile=NaN because n_iteration=0
    QC.nonburst.percentile = nan;                                           % percentile=NaN because n_iteration=0
    QC.total.percentile = nan;                                              % percentile=NaN because n_iteration=0
    QC.ch1.percentile = nan;
    QC.ch2.percentile = nan;
end

end
