function FiringPattern = FiringPatternFeatures(spikeTrain, tRange, netBursts, reference)

% function to measure the firing pattern of each unit outside the network bursts, for the classification of the neurons (new, November 2026)
% It is meant for the baseline condition: one value per unit for the relative firing rate, the local variation (LV), the coefficient of
% variation of the inter-spike interval (CV) and the burst index. Each is measured only on the spikes OUTSIDE the primary network bursts
% (each burst from its first to its last spike, widened by GUARD sec on each side): in the organoids about half of the spikes of a unit fall
% inside network bursts, so the same measures on the whole recording mostly restate the role of the unit in the network bursts
% (participation, onset order), which are the outcomes to compare between the classes. Measured outside the bursts, LV, CV and the burst
% index are independent of those roles (7,963 units of 58 organoids, November 2026)
% The consecutive intervals of a pair, and the pairs of spikes of the autocorrelogram, must lie within the same period between two bursts:
% no interval spans a burst
% Choice of the measures (November 2026): LV (Shinomoto et al. 2003; = the "lvr" of "SpikeFeatures" on the whole train), CV of the ISI,
% burst index from the autocorrelogram (Royer et al. 2012: spikes at 3-5 msec against 200-300 msec) and the firing rate relative to the
% other units of the same organoid (the overall activity of each organoid is removed). Left out: the CV2 and LvR (same information as LV),
% the fit of the autocorrelogram of CellExplorer (unreliable at the rates of these recordings), the entropy (follows the rate), the Fano
% factor of the counts and the metrics of the single-unit bursts (follow the rate and the CV)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   spikeTrain:     structure array with one element per unit (as PerCondition(cond).spikeTrain of "Analyze_spikeTrain"), with the field
%                   .timeStamp: time of each spike, in msec
%                   (new, November 2026) or a firing pattern already measured (output of this function): only its relativeRate is computed
%                   again, against the units of "reference" (with the MIN_REFERENCE it was measured with); tRange and netBursts are not used
%   tRange:         [start end] of the analysis window, in sec (as PerCondition(cond).window)
%   netBursts:      primary network bursts of the same window: the structure NetworkBurstDetected (its IndivBurst(j).start and .last,
%                   in sec) or a matrix with one row per burst [start last], in sec; [] if there are no network bursts (the whole window
%                   is used). (revised November 2026; user's choice) The verdict of the pipeline is followed: when NetworkBurstDetected.nBurst
%                   is NaN (the candidate bursts were rejected by the final check of "NetworkBurstDetection") or 0, the recording has no
%                   network bursts and the whole window is used, even if IndivBurst still holds the candidates
%   reference:      (optional) logical vector, one value per unit: the units that set the median rate of the organoid for "relativeRate"
%                   (e.g., the somatic units); all units if omitted or []

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The 4 features (one value per unit; NaN when the unit has less than MIN_SPIKES spikes outside the bursts):
% FiringPattern.relativeRate(i)     log10 of the firing rate outside the bursts minus the median of the same value across the reference
%                                   units with spikes outside the bursts (NaN for every unit if there are less than MIN_REFERENCE of them)
% FiringPattern.lv(i)               local variation of the consecutive intervals: 3/n * sum[((I(k)-I(k+1)) / (I(k)+I(k+1)))^2]
%                                   (=1 for a Poisson train, <1 regular, >1 bursty)
% FiringPattern.cv(i)               coefficient of variation of the inter-spike intervals (std / mean)
% FiringPattern.burstIndex(i)       log10 of the ratio of the autocorrelogram at 3-5 msec to the autocorrelogram at 200-300 msec (mean count
%                                   per 1-msec bin of each range, plus PSEUDO_COUNT to avoid log(0)); >0 when the unit fires bursts of
%                                   spikes a few msec apart
% Information about each unit:
% FiringPattern.rateOutside(i)      firing rate outside the bursts, in 1/sec
% FiringPattern.nSpikes(i)          number of spikes in the analysis window
% FiringPattern.nSpikesOutside(i)   number of spikes outside the bursts
% FiringPattern.enoughSpikes(i)     =true if nSpikesOutside >= MIN_SPIKES (the 4 features are measured)
% FiringPattern.timeOutside         time outside the bursts, in sec (the same for all the units)
% FiringPattern.periods             [start end] of each period outside the bursts, in msec
% FiringPattern.features            names of the 4 features (fields above), in this order
% FiringPattern.parameters          parameters used

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GUARD = 0.5;                        % margin added before and after each network burst, in sec (the rise and the tail of the burst are also
                                    % left out)
MIN_SPIKES = 50;                    % minimum number of spikes outside the bursts to measure the features (65% of the units in the test of
                                    % November 2026; the features of the two halves of the baseline agree with Spearman 0.74-0.93)
SHORT_LAGS = [3 5];                 % lags of the autocorrelogram for the bursts, in msec (Royer et al. 2012)
LONG_LAGS = [200 300];              % lags of the autocorrelogram for the reference, in msec (Royer et al. 2012)
PSEUDO_COUNT = 0.5;                 % added to the mean count per bin of each range of the autocorrelogram, to avoid log(0)
MIN_REFERENCE = 5;                  % minimum number of reference units with spikes outside the bursts to compute relativeRate
apply_parameters('FiringPatternFeatures');                                  % values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

FEATURES = {'relativeRate', 'lv', 'cv', 'burstIndex'};
K_MSEC_SEC = 1/1000;                % factor to rescale from msec to sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (new, November 2026) Firing pattern already measured: relative rate against other reference units

if isstruct(spikeTrain) && isscalar(spikeTrain) && isfield(spikeTrain, 'rateOutside') && isfield(spikeTrain, 'enoughSpikes')
    FiringPattern = spikeTrain;
    n_units = numel(FiringPattern.rateOutside);
    if nargin < 4 || isempty(reference)
        reference = true(n_units, 1);
    end
    reference = logical(reference(:));
    if numel(reference) ~= n_units
        error("FiringPatternFeatures: 'reference' must have one value per unit (%d units, %d values)", n_units, numel(reference));
    end
    if isfield(FiringPattern, 'parameters') && isfield(FiringPattern.parameters, 'MIN_REFERENCE')
        MIN_REFERENCE = FiringPattern.parameters.MIN_REFERENCE;
    end
    FiringPattern.relativeRate = relative_rate(FiringPattern.rateOutside(:), reference, logical(FiringPattern.enoughSpikes(:)), MIN_REFERENCE);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Periods outside the network bursts

n_units = numel(spikeTrain);
if nargin < 4 || isempty(reference)
    reference = true(n_units, 1);
end
reference = logical(reference(:));
if numel(reference) ~= n_units
    error("FiringPatternFeatures: 'reference' must have one value per unit (%d units, %d values)", n_units, numel(reference));
end

if isstruct(netBursts)                                                      % structure NetworkBurstDetected
    % (revised November 2026; user's choice) nBurst = NaN or 0: no network bursts, as for the rest of the pipeline (before, the windows of the
    % candidates rejected by the final check of "NetworkBurstDetection" were left out too: 9-70% of the baseline in 8 of 59 recordings)
    rejected = isfield(netBursts, 'nBurst') && isscalar(netBursts.nBurst) && (isnan(netBursts.nBurst) || netBursts.nBurst == 0);
    if ~rejected && isfield(netBursts, 'IndivBurst') && ~isempty(netBursts.IndivBurst) && isfield(netBursts.IndivBurst, 'start')
        bursts = [[netBursts.IndivBurst.start]' [netBursts.IndivBurst.last]'];
    else
        bursts = zeros(0, 2);                                               % no network bursts
    end
elseif isempty(netBursts)
    bursts = zeros(0, 2);
else
    bursts = netBursts;
end

w0 = tRange(1) / K_MSEC_SEC;    w1 = tRange(2) / K_MSEC_SEC;                % analysis window, in msec
bursts = sortrows(bursts / K_MSEC_SEC);                                     % in msec
periods = zeros(0, 2);  current = w0;
for j=1:size(bursts, 1)
    a = max(w0, bursts(j,1) - GUARD / K_MSEC_SEC);                          % start of the widened burst
    b = min(w1, bursts(j,2) + GUARD / K_MSEC_SEC);                          % end of the widened burst
    if a > current
        periods(end+1,:) = [current a];                                     %#ok<AGROW> period before this burst
    end
    current = max(current, b);
end
if current < w1
    periods(end+1,:) = [current w1];                                        % period after the last burst
end
time_outside = sum(diff(periods, 1, 2)) * K_MSEC_SEC;                       % in sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Features of each unit

n_spikes = zeros(n_units, 1);    n_outside = zeros(n_units, 1);
lv = nan(n_units, 1);   cv = nan(n_units, 1);   burst_index = nan(n_units, 1);
edges = 0:1:LONG_LAGS(2);                                                   % 1-msec bins of the autocorrelogram
short_bins = (edges(1:end-1) >= SHORT_LAGS(1)) & (edges(1:end-1) < SHORT_LAGS(2));
long_bins = (edges(1:end-1) >= LONG_LAGS(1)) & (edges(1:end-1) < LONG_LAGS(2));

for i=1:n_units
    t = sort(spikeTrain(i).timeStamp(:));
    t = t(t >= w0 & t <= w1);
    n_spikes(i) = numel(t);
    isi = [];   pairs = zeros(0, 2);    acg = zeros(1, numel(edges) - 1);
    for p=1:size(periods, 1)
        tp = t(t >= periods(p,1) & t < periods(p,2));                       % spikes of this period
        n_outside(i) = n_outside(i) + numel(tp);
        d = diff(tp);
        isi = [isi; d];                                                     %#ok<AGROW>
        if numel(d) > 1
            pairs = [pairs; d(1:end-1) d(2:end)];                           %#ok<AGROW> consecutive intervals within the period
        end
        for k=1:numel(tp)-1                                                 % autocorrelogram: all the pairs within LONG_LAGS(2) msec
            lag = tp(k+1:end) - tp(1:end-k);
            lag = lag(lag <= LONG_LAGS(2));
            if isempty(lag)
                break
            end
            acg = acg + histcounts(lag, edges);
        end
    end
    if n_outside(i) < MIN_SPIKES
        continue
    end
    if size(pairs, 1) >= 2
        lv(i) = 3 * mean(((pairs(:,1) - pairs(:,2)) ./ (pairs(:,1) + pairs(:,2))).^2);
    end
    if numel(isi) >= 3
        cv(i) = std(isi) / mean(isi);
    end
    burst_index(i) = log10((mean(acg(short_bins)) + PSEUDO_COUNT) / (mean(acg(long_bins)) + PSEUDO_COUNT));
end

if time_outside > 0
    rate_outside = n_outside / time_outside;
else
    rate_outside = nan(n_units, 1);                                         % the bursts cover the whole window
end
enough = n_outside >= MIN_SPIKES;
rel_rate = relative_rate(rate_outside, reference, enough, MIN_REFERENCE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output

FiringPattern.relativeRate = rel_rate;
FiringPattern.lv = lv;
FiringPattern.cv = cv;
FiringPattern.burstIndex = burst_index;
FiringPattern.rateOutside = rate_outside;
FiringPattern.nSpikes = n_spikes;
FiringPattern.nSpikesOutside = n_outside;
FiringPattern.enoughSpikes = enough;
FiringPattern.timeOutside = time_outside;
FiringPattern.periods = periods;
FiringPattern.features = FEATURES;
FiringPattern.parameters = struct('GUARD', GUARD, 'MIN_SPIKES', MIN_SPIKES, 'SHORT_LAGS', SHORT_LAGS, 'LONG_LAGS', LONG_LAGS, ...
    'PSEUDO_COUNT', PSEUDO_COUNT, 'MIN_REFERENCE', MIN_REFERENCE);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rr = relative_rate(rate_outside, reference, enough, min_reference)
% log10 of the rate outside the bursts minus its median across the reference units with spikes outside the bursts (NaN for every unit with
% less than min_reference of them, and for the units without enough spikes) (moved here, November 2026)
log_rate = log10(rate_outside);
log_rate(rate_outside == 0) = NaN;
use = reference & rate_outside > 0;
if sum(use) >= min_reference
    rr = log_rate - median(log_rate(use));
else
    rr = nan(numel(rate_outside), 1);
end
rr(~enough) = NaN;
end
