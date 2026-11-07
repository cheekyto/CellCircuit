function participationSurrogate = participationSurrogate(spikeData, range, bin, smoothWindow, n_surrogate, seed)

% function to build the null distributions of the population participation signal (fraction of units with at least one spike in each
% bin), with surrogates in which the spike train of each unit is circularly shifted by an independent random offset within "range" (each
% unit keeps its spikes and its inter-spike intervals; the timing between units is broken)
% It is used by "PopulationActivityDetection" (option 4 of "NetworkBurstDetection", opt_norm = 2), with two null distributions:
%   - the maximum of each surrogate (max statistic, for the large events and the secondary events: a peak above the (1-alpha) quantile
%     of the maxima has a probability alpha of happening by chance anywhere in "range")
%   - all the bins of all the surrogates, pooled (for the small events, as Bocchio et al. 2020, Nat Commun 11:4559)
% (new, November 2026)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   spikeData:      structure array with the selected units (field .timeStamp, in msec)
%   range:          [start end] of the analysis window, in sec
%   bin:            size of the bins, in sec (the same bins as "bin_spike_time")
%   smoothWindow:   length of the Gaussian smoothing window, in number of bins (the same smoothing as the real signal)
%   n_surrogate:    number of surrogates
%   seed:           seed for the random offsets (the same seed gives the same surrogates)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   participationSurrogate.maxima(b):       maximum of the smoothed participation signal of the surrogate b (column)
%   participationSurrogate.pooledCount(j):  number of bins (all the surrogates together) with a value in [(j-1), j)*HIST_WIDTH (column)
%   participationSurrogate.histWidth:       width of the bins of "pooledCount" (HIST_WIDTH)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                            % factor to rescale from msec to sec
HIST_WIDTH = 1e-4;                              % resolution of the pooled distribution (the participation signal goes from 0 to 1)
MAX_CACHE = 200;                                % maximum number of null distributions kept in memory (the screen of bin sizes asks for the same one more than once)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spikes of the selected units within the range

persistent cache_keys cache_values
if isempty(cache_keys)
    cache_keys = {};    cache_values = {};
end

n_units = numel(spikeData);
r0 = range(1)/K_MSEC_SEC;                                                   % start of the range, in msec
T = (range(2) - range(1))/K_MSEC_SEC;                                       % length of the range, in msec
bin_ms = bin/K_MSEC_SEC;                                                    % bin size, in msec
if T > bin_ms                                                               % bins as in "bin_spike_time"
    edges = r0 : bin_ms : range(2)/K_MSEC_SEC;
else
    edges = [r0, r0 + bin_ms];
end
n_bins = numel(edges) - 1;

times = cell(n_units,1);    unit = cell(n_units,1);
for i=1:n_units
    t = spikeData(i).timeStamp(:);
    t = t(t >= r0 & t <= r0 + T);
    times{i} = t - r0;                                                      % time from the start of the range, in msec
    unit{i} = repmat(i, numel(t), 1);
end
times = vertcat(times{:});      unit = vertcat(unit{:});

key = sprintf('%d|%d|%.10g|%.10g|%.10g|%d|%d|%d|%.12g', n_units, numel(times), r0, T, bin_ms, smoothWindow, n_surrogate, seed, sum(times));
hit = find(strcmp(cache_keys, key), 1);
if ~isempty(hit)                                                            % already computed in this run
    participationSurrogate = cache_values{hit};
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Surrogates: each unit is circularly shifted by its own random offset

stream = RandStream('twister', 'Seed', seed);                               % local random stream (the global generator is not modified)
offsets = rand(stream, n_units, n_surrogate) * T;                           % offset of each unit (row) in each surrogate (column), in msec
last_edge = edges(end) - r0;                                                % end of the last bin, from the start of the range
n_hist = round(1/HIST_WIDTH);
maxima = nan(n_surrogate,1);
pooled = zeros(n_hist,1);
for b=1:n_surrogate
    shifted = mod(times + offsets(unit,b), T);                              % shifted spike times, in [0, T)
    inside = shifted <= last_edge;                                          % (the incomplete bin at the end of the range is not used, as in the real signal)
    bin_idx = min(floor(shifted(inside)/bin_ms) + 1, n_bins);               % bin of each spike (the last bin includes its right edge, as "histcounts")
    pairs = unique((unit(inside) - 1)*n_bins + bin_idx);                    % bins with at least one spike, for each unit
    active = accumarray(mod(pairs - 1, n_bins) + 1, 1, [n_bins, 1]);        % number of active units in each bin
    signal = smoothdata(active/max(n_units,1), 'gaussian', smoothWindow);   % participation signal, smoothed as the real one
    maxima(b) = max(signal);
    hist_idx = min(max(floor(signal/HIST_WIDTH) + 1, 1), n_hist);           % pooled distribution of all the bins
    pooled = pooled + accumarray(hist_idx, 1, [n_hist, 1]);
end
participationSurrogate.maxima = maxima;
participationSurrogate.pooledCount = pooled;
participationSurrogate.histWidth = HIST_WIDTH;

if numel(cache_keys) >= MAX_CACHE
    cache_keys(1) = [];     cache_values(1) = [];
end
cache_keys{end+1} = key;      cache_values{end+1} = participationSurrogate;

end
