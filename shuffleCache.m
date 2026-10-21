function entry = shuffleCache(action, ch, data, range, t_step, n_shuffle, BurstDetected0, seed)

% function to keep the shuffled spike trains of each cluster (method 3 of "ShuffleSpikeTrain") during one run of "NetworkBurstDetection"
% (revised October 2026)
% "findChannel_in_Peak" assesses the participation of each cluster with the same N_SHUFFLING shuffles every time (same seed for each
% cluster), and it is called many times per run (17 times on the reference recording). With the cache, the shuffles of each cluster are
% generated only once per run. Same results.
% The shuffles are kept as the bin index of each spike (uint32: 4 bytes per spike and shuffle); above MAX_CACHE_MB, the new clusters are
% not kept (their shuffles are generated every time, as before)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% action:               ='reset', empty the cache (at the beginning and at the end of "NetworkBurstDetection")
%                       ='get', shuffles of the cluster ch (generated now, or taken from the cache)
% ch:                   index of the cluster
% data, range, t_step, n_shuffle, BurstDetected0, seed:   same inputs as "ShuffleSpikeTrain" (method 3), range in msec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% entry.keys(i):        key of the spike i of all the shuffles: (shuffle-1)*entry.M + bin index (uint32, ascending: shuffle 1, then 2, ...)
% entry.M:              distance between the keys of consecutive shuffles (n_bin + 2, so the shuffles never mix)
% entry.bin, entry.ini: bin size and start of the range, in msec: the time of the bin b is (b-0.5)*bin + ini, as in "ShuffleSpikeTrain"
% entry.n_bin:          total number of bins;       entry.n_shuffle: number of shuffles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MAX_CACHE_MB = 1000;                            % maximum memory for the cache, in MB (~175 MB for the reference recording, 87,000 spikes)
apply_parameters('shuffleCache');                                           % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent cache total_bytes
if isempty(total_bytes)
    cache = {};     total_bytes = 0;
end

switch action
    case 'reset'
        cache = {};     total_bytes = 0;
        entry = [];

    case 'get'
        % inputs that determine the shuffles (to make sure that a kept entry is never used with different data)
        key = [numel(data), sum(data), data(1), data(end), range(1), range(2), t_step, n_shuffle, seed, ...
               size(BurstDetected0.spikes,1), sum(BurstDetected0.spikes(:,3))];
        if ch <= numel(cache) && ~isempty(cache{ch}) && isequal(cache{ch}.key, key)
            entry = cache{ch};                                              % already generated in this run
            return
        end
        BurstFeature0 = BurstFeatures(BurstDetected0);                      % characterize the bursts of the cluster (as in "findChannel_in_Peak")
        [~, sb] = ShuffleSpikeTrain(data, range, t_step, 3, n_shuffle, 0, BurstDetected0, BurstFeature0, seed);
        M = sb.n_bin + 2;
        n_spikes = cellfun(@numel, sb.bins(:));                             % number of spikes of each shuffle
        if n_shuffle*M < double(intmax('uint32'))                           % usual case: the keys fit in uint32 (4 bytes each)
            key_class = 'uint32';
        else                                                                % (very long recordings with very small bins) exact integers as double
            key_class = 'double';
        end
        shuffle_id = repelem(cast(0:n_shuffle-1, key_class)', n_spikes);
        entry.keys = shuffle_id*cast(M, key_class) + cast(vertcat(sb.bins{:}, zeros(0,1,'uint32')), key_class);
        entry.M = M;    entry.bin = sb.bin;     entry.ini = sb.ini;     entry.n_bin = sb.n_bin;     entry.n_shuffle = n_shuffle;
        entry.key = key;
        bytes = numel(entry.keys)*(4 + 4*strcmp(key_class,'double'));
        if total_bytes + bytes <= MAX_CACHE_MB*1e6                          % keep it, if there is enough memory left
            cache{ch} = entry;
            total_bytes = total_bytes + bytes;
        end
end

end
