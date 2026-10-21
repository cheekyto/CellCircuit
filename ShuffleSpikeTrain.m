function [ShuffleSpikeTrain, shuffle_bins] = ShuffleSpikeTrain(signal, range, t_step, method, n_shuffle, refractory_period, BurstDetected0, BurstFeature0, seed)

% function to shuffle a spike train
% (revised October 2026: method 2 places the spikes directly, instead of drawing one random number per bin -same random process,
%  statistically equivalent results, ~100x faster-; no "clear" inside the loops; burst labels assigned in one pass (method 2);
%  method 3 no longer resets the burst label of the first spike of iteration 1 -bug-)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% signal(i,1):                  time stamp of the spike i, in msec
% range:                        range to generated the new spike train (range(1)=start; range(2)=end), in msec
% t_step:                       temporal resolution - sampling rate, in msec
% n_shuffle:                    number of iterations       
% refractory_period:            absolute refractory period, in msec (usually =0 because I do not care)             
% BurstDetected0:               output from the function "BurstDetection" corresponding to the signal to be shuffled 
% BurstFeature0:                output from the function "BurstFeatures" corresponding to the signal to be shuffled
% seed:                         (optional) seed for the random number generator (default: DEFAULT_SEED)
%                               same seed + same inputs = exactly the same shuffled spike trains

% method:                       method for spike distribution:    
%                               = 1,    uniform distribution
%                               = 2,    biphasic (burst and non-bust); keep the temporal position of the bursts and shuffle the spikes within them 
%                               = 3,    shuffling of individual bursts (keep the number and duration, but with uniform distribuion within the burst) and spikes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ShuffleSpikeTrain{j}(i,1):     time stamp of the spike i for newly shuffled spike train, iteration j, in msec
% ShuffleSpikeTrain{j}(i,2):     indicates if the spike i belongs to a burst (=0, No; =1, Yes)
% shuffle_bins:                  (optional, only method 3; revised October 2026) the same shuffled trains, as bin indexes (used by "shuffleCache")
%   shuffle_bins.bins{j}:        bin index of each spike of the iteration j (uint32); its time is (bin_index-0.5)*bin + ini, in msec
%   shuffle_bins.bin, .ini:      bin size and start of the range, in msec;     shuffle_bins.n_bin: total number of bins

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BIN = 50;                        % maximum bin size, for option #3, in ms (default: 50)     
BURST_METHOD = 1;                % burst detection method (=1      Maximal Interval)
DEFAULT_SEED = 0;                % seed for the random number generator when the input "seed" is not provided
apply_parameters('ShuffleSpikeTrain');                                      % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate the newly shuffled spike train

if nargin < 9 || isempty(seed)
    seed = DEFAULT_SEED;
end
shuffle_bins = [];                                                          % (only filled by method 3, when it is requested)
stream = RandStream("twister", "Seed", seed);                               % private random number generator (reproducible, and not affected by other rand calls)
ini = range(1);                                                             % start point of the time window, in msec
fin = range(2);                                                             % last point of the time window, in msec
n_bin = round((fin-ini)/t_step);                                            % total number of bins or temporal positions    

if size(signal,2)>1 || size(signal,1)<1                                     % if the variable signal is not a column array or it does not contain at least 1 element
    error("The input variable for *ShuffleSpikeTrain* has the wrong dimensions (must be a column array")
end

switch method                                                               % Select the shuffling method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1          % UNIFORM distribution

        n_spikes = numel(signal);                                           % total number of spikes in the train
        ShuffleSpikeTrain = cell(1, n_shuffle);                             % initialize the output variable

        for it=1:n_shuffle                                                  % across all iterations
            aux_train = n_bin.*rand(stream, n_spikes,1);                    % array of "n_spikes" random numbers in the range [0 n_bin]
            aux_sorted = sort(aux_train)*t_step;                            % spike time stamp (in ascending order), in msec
            aux_diff = diff(aux_sorted);                                    % inter-spike interval, in msec
            aux = find(aux_diff<refractory_period);                         % identify the isi that violates the refractory period
            if numel(aux)                                                   % there is at least one spike violating the absolute refractory period 
               aux_sorted(aux) = [];                                        % eliminate the faulty spikes 
            end
            ShuffleSpikeTrain{it}(:,1) = aux_sorted+ini;                    % the first column is the time stamp
            BurstDetected = BurstDetection(aux_sorted, BURST_METHOD);       % run the burst detection algorithm
            ShuffleSpikeTrain{it}(:,2) = zeros(size(aux_sorted));           % in principle, label all spikes as non-bursting    
            burst_idx = find(BurstDetected.spikes(:,3)>0);                  % identify all the spikes within burst
            if numel(burst_idx)>0                                           % if there is at least one spike within burst
                ShuffleSpikeTrain{it}(burst_idx,2) = 1;                     % the second column indicates that the spike belongs to a burst
            end                                                             % (revised October 2026: no "clear" inside the loop)
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2          % biphasic (bust and non-burst)

        % (revised October 2026) Same random process as before: each temporal bin (sampling step) gets a spike with probability t_step/isi
        % (isi within bursts inside the windows of the original bursts; isi outside bursts elsewhere). Before, one random number was drawn
        % for EVERY bin of the recording in every iteration (e.g., 53 million for 45 min at 20 kHz); now the spikes are placed directly,
        % drawing the number of bins between consecutive spikes (geometric distribution), which is ~300x faster. The shuffled spike trains
        % are statistically equivalent, but not identical to the ones generated before (different random numbers)
        isi_burst = BurstFeature0.neuronMetric.isi.mean;                    % average isi for spikes within burst, in msec
        isi_nonburst = BurstFeature0.neuronMetric.nonBursting.isi.mean;     % average isi for spikes outside of bursts, in msec
        p_nonburst = t_step/isi_nonburst;                                   % probability of a spike in a (regular) time bin
        p_burst = t_step/isi_burst;                                         % probability of a spike in a time bin within the window corresponding to a burst
        n_burst0 = BurstFeature0.nBurst;                                    % total number of bursts in the original data
        ini_b = zeros(n_burst0,1);  fin_b = zeros(n_burst0,1);              % start and end of the window corresponding to each burst, in msec
        ini_idx = zeros(n_burst0,1);  fin_idx = zeros(n_burst0,1);          % start and end of the window corresponding to each burst, in bins
        for k=1:n_burst0                                                    % across all bursts in the original data
            ini_b(k) = BurstDetected0.IndivBurst(k).start - isi_burst/2;    % start point for the window corresponding to burst k, in msec
            fin_b(k) = BurstDetected0.IndivBurst(k).last + isi_burst;       % last point for the window corresponding to burst k, in msec
            ini_idx(k) = max(round((ini_b(k)-ini)/t_step), 1);              % index for ini_b(k)
            fin_idx(k) = min(round((fin_b(k)-ini)/t_step), n_bin);          % index for fin_b(k)
        end
        [burst_seg, nonburst_seg] = bin_segments(ini_idx, fin_idx, n_bin);  % segments of consecutive bins within / outside the burst windows ([first last] bin of each segment)

        ShuffleSpikeTrain = cell(1, n_shuffle);                             % initialize the output variable
        for it=1:n_shuffle                                                  % across all iterations
            aux_positive = sort([random_bins(stream, burst_seg, p_burst); random_bins(stream, nonburst_seg, p_nonburst)]);   % all the temporal points with 1 spike (Poisson-like spike train)
            aux_train = nan(numel(aux_positive),1);
            aux_train(:,1) = aux_positive(:)*t_step+ini;                    % new shuffled spike train, in msec
            aux_diff = diff(aux_train);                                     % inter-spike interval for the new shuffled spike train, in msec
            aux = find(aux_diff<refractory_period);                         % identify the events that violate the refractory period    
            if numel(aux)                                                   % there is at least one spike violating the absolute refractory period 
               aux_train(aux) = [];                                         % eliminate the faulty spikes 
            end
            
            ShuffleSpikeTrain{it} = zeros(numel(aux_train),2);              % initiate the output variable
            ShuffleSpikeTrain{it}(:,1) = aux_train;                         % the first column is the time stamp for iteration it, in msec
            if BurstFeature0.nBurst > 0                                     % the second column indicates if the unit is within the temporal window of a burst from the original data
                n_windows = count_below(ini_b, aux_train, 0) - count_below(fin_b, aux_train, 1);   % number of burst windows that contain each spike: ini_b(k) < t < fin_b(k)
                ShuffleSpikeTrain{it}(n_windows>0,2) = 1;                   % (we are not running "BurstDetection", so it might not be an actual burst in the shuffled spike train)
            end                                                             % (revised October 2026: one pass across all spikes, instead of one search per burst; same result)
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3          % re-shuffling of individual bursts and spikes
    
        isi_burst = BurstFeature0.neuronMetric.isi.mean;                    % average isi for spikes within burst, in msec
        bin = min(BIN, isi_burst);                                          % bin size, in msec
        n_bin = round((fin-ini)/bin);                                       % total number of bins or temporal positions    
        ShuffleSpikeTrain = cell(1, n_shuffle);                             % initialize the output variable
        keep_bins = nargout > 1;                                            % (revised October 2026) also return the bin indexes
        if keep_bins
            shuffle_bins = struct('bins', {cell(1, n_shuffle)}, 'bin', bin, 'ini', ini, 'n_bin', n_bin);
        end
        n_burst0 = BurstFeature0.nBurst;                                    % total number of bursts in the original spike train
        bin_burst = zeros(n_burst0,1);                                      % duration of each burst, in number of bins (the same in every iteration)
        for k=1:n_burst0
            duration = BurstDetected0.IndivBurst(k).last - BurstDetected0.IndivBurst(k).start + bin/2;  % duration of the burst k, in msec
            bin_burst(k) = ceil(duration/bin);                              % duration of the burst, in number of bins 
        end
        n_single = sum(BurstDetected0.spikes(:,3)==0);                      % total number of individual spikes (not within a burst) in the original spike train

        % (revised October 2026) same random numbers and same placement rules as before, but: each iteration uses its own row "spikes_it" 
        % (before, the row "it" of a n_shuffle x n_bin matrix); the individual spikes are placed all at once in the first free bins of the 
        % random order (before, one by one in a while loop); and only the beginning of the random order is sorted (see below)
        burst_bins = sum(min(bin_burst+1, n_bin));                          % number of bins occupied by the bursts
        free_fraction = max(0.05, 1 - (burst_bins + n_single)/n_bin);       % approximate proportion of free bins
        n_needed = ceil(1.5*(n_burst0 + n_single)/free_fraction) + 100;     % expected number of positions of the random order that will be used (+ margin)
        cut_value = n_needed/n_bin;                                         % the bins with random number < cut_value are the first ~n_needed positions of the random order
        for it=1:n_shuffle                                                  % across all shuffling iterations
            random = rand(stream, n_bin,1);                                 % array of "n_bin" random numbers (uniformly distributed random numbers in the interval [0,1])
            if cut_value < 0.5                                              % only the bins with the smallest random numbers are sorted: exactly the first positions of
                candidates = find(random < cut_value);                      % "sort(random)" (also with ties: equal values are all below or all above the cut), much faster
                [~, order] = sort(random(candidates));
                [spikes_it, complete] = place_events(candidates(order), n_bin, bin_burst, n_single);
            else
                complete = false;
            end
            if ~complete                                                    % the first positions were not enough (or most of the bins are needed): full random order, as before
                [~, bin_order] = sort(random);                              % index of the temporal bins, in a random order
                spikes_it = place_events(bin_order, n_bin, bin_burst, n_single);
            end
            
            assigned_bins = find(spikes_it > 0);                            % bins with a spike assigned (1 = single spike; 2 = burst)
            ShuffleSpikeTrain{it} = [(assigned_bins(:)-0.5)*bin+ini, ...    % first column: locate the spike at the center of the bin, in msec
                                     double(spikes_it(assigned_bins)' == 2)];   % second column: =1, the spike belongs to a burst; =0, it is an individual spike
            if keep_bins
                shuffle_bins.bins{it} = uint32(assigned_bins(:));           % the same spikes, as bin indexes
            end
        end                                                                 % end of the shuffling iterations
end                                                                         % end of the switch to select the shuffling method
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [burst_seg, nonburst_seg] = bin_segments(ini_idx, fin_idx, n_bin)

% function to split the bins 1..n_bin into segments of consecutive bins within the burst windows (union of all windows) and outside of them
% (revised October 2026)
% burst_seg(s,:) / nonburst_seg(s,:):   [first last] bin of the segment s

valid = ini_idx <= fin_idx;                                                 % windows with at least one bin (as before, with "threshold(ini_idx:fin_idx)")
ini_idx = ini_idx(valid);   fin_idx = fin_idx(valid);
[ini_idx, order] = sort(ini_idx);   fin_idx = fin_idx(order);
burst_seg = zeros(0,2);
for k=1:numel(ini_idx)                                                      % across all windows, merging the ones that overlap (or touch each other)
    if ~isempty(burst_seg) && ini_idx(k) <= burst_seg(end,2)+1
        burst_seg(end,2) = max(burst_seg(end,2), fin_idx(k));
    else
        burst_seg(end+1,:) = [ini_idx(k) fin_idx(k)];                       %#ok<AGROW>
    end
end
nonburst_seg = [[1; burst_seg(:,2)+1], [burst_seg(:,1)-1; n_bin]];          % gaps between the burst windows
nonburst_seg = nonburst_seg(nonburst_seg(:,1) <= nonburst_seg(:,2), :);     % (only the non-empty ones)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bins = random_bins(stream, seg, p)

% function to draw the bins (within the segments "seg") that get a spike, when every bin independently gets a spike with probability p
% The segments are concatenated in a single sequence of bins; in such a sequence of independent trials, the number of bins between 
% consecutive spikes follows a geometric distribution, so the spikes can be placed directly (without drawing one random number per bin)
% (revised October 2026)

L_seg = seg(:,2) - seg(:,1) + 1;                                            % number of bins in each segment
L = sum(L_seg);                                                             % total number of bins
if L == 0 || ~(p > 0)                                                       % no bins, or probability zero (or undefined, e.g., no spikes outside bursts) -> no spikes
    bins = zeros(0,1);
    return
end
if p >= 1                                                                   % every bin gets a spike
    pos = (1:L)';
else
    log_q = log1p(-p);                                                      % log(1-p)
    pos = zeros(0,1);   last = 0;                                           % position of the spikes in the sequence of bins, and position of the last one
    while true
        n_draw = ceil((L-last)*p + 5*sqrt((L-last)*p) + 10);                % expected number of spikes in the remaining bins (+ margin)
        new_pos = last + cumsum(floor(log(rand(stream, n_draw, 1))/log_q) + 1);    % geometric distance (in bins) between consecutive spikes
        if new_pos(end) >= L                                                % the end of the sequence was reached
            pos = [pos; new_pos(new_pos <= L)];                             %#ok<AGROW>
            break
        end
        pos = [pos; new_pos];                                               %#ok<AGROW> (rarely needed: more spikes than expected)
        last = new_pos(end);
    end
end
cum_L = cumsum(L_seg);                                                      % last position (in the sequence) of each segment
s = count_below(cum_L, pos, 0) + 1;                                         % segment of each spike: number of segments that end before it, + 1
bins = seg(s,1) + pos - (cum_L(s) - L_seg(s)) - 1;                          % bin of each spike

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spikes_it, complete] = place_events(bin_order, n_bin, bin_burst, n_single)

% function to place the bursts and the individual spikes of one shuffling iteration (method 3), following the random order of the bins
% (revised October 2026: same rules as the original code)
% bin_order:            random order of the bins (it can be only the first positions of the full random order)
% spikes_it(j):         =0, empty bin; =1, individual spike; =2, burst
% complete:             =false if "bin_order" was not long enough (only possible when it is not the full random order)

n_order = numel(bin_order);
is_full = n_order == n_bin;                                                 % full random order (as in the original code)
spikes_it = zeros(1, n_bin);
complete = true;

% first, assigned the bursts to specific locations (all the bin with bursts will be labeled "2")
assigned = 1;                                                               % index that points out to the bin under analysis (position in the random order)
for k=1:numel(bin_burst)                                                    % across all burst (from the original spike train)
    while assigned <= n_order                                               % until the assignment is done (or the full length of the recordings was covered)
        ini_b = bin_order(assigned);                                        % index of the first bin for the burst k
        fin_b = min(n_bin, ini_b+bin_burst(k));                             % index of the last bin for the burst k
        if ~any(spikes_it(ini_b:fin_b))                                     % there is no spike within the window ini:fin
           spikes_it(ini_b:fin_b) = 2;                                      % assign spikes to all bins in that window (2 means it is a burst)
           break                                                            % terminate the assignement for burst k
        end
        assigned = assigned + 1;                                            % there are already spikes in the window ini:fin -> keep looking into the next bin (random order)
    end
    if assigned > n_order && ~is_full                                       % ran out of positions before the end of the full random order
        complete = false;
        return
    end
end

% second, assigned the individual spikes: each one to the next free bin in the random order (starting from "assigned")
if assigned <= n_order && n_single > 0
    candidates = bin_order(assigned:end);                                   % bins not tried yet, in the random order
    free_bins = candidates(spikes_it(candidates) == 0);                     % free bins, in the random order
    if numel(free_bins) < n_single && ~is_full                              % not enough free bins in the first positions of the random order
        complete = false;
        return
    end
    spikes_it(free_bins(1:min(n_single, numel(free_bins)))) = 1;            % assign the individual spikes (1 means it is a single spike)
elseif n_single > 0 && ~is_full                                             % no positions left for the individual spikes
    complete = false;
end

end
