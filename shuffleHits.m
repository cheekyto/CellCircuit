function hits = shuffleHits(entry, w_ini, w_fin)

% function to count, for each shuffled spike train, the number of temporal windows with at least one spike (revised October 2026)
% Same result as calling "windowCounts" with each shuffled train, one after the other (see "findChannel_in_Peak"), but all the shuffles at
% once: each spike is identified by an integer key (shuffle and bin index), so the comparisons with the windows are exact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% entry:                shuffled spike trains of a cluster, from "shuffleCache"
% w_ini(k), w_fin(k):   start and end of the window k, in msec; a spike t is within the window if w_ini(k) <= t < w_fin(k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hits(i):              number of windows in which the shuffled train i has at least one spike (column array)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_shuffle = entry.n_shuffle;
n_win = numel(w_ini);
lo = first_bin(entry, w_ini(:));                                            % first bin whose spike time is >= w_ini(k)
hi = first_bin(entry, w_fin(:));                                            % first bin whose spike time is >= w_fin(k)
offsets = (0:n_shuffle-1)*entry.M;                                          % key of the bin 0 of each shuffle
edges = [0; double(entry.keys); n_shuffle*entry.M];                         % all keys are integers in (0, n_shuffle*M)
below_ini = discretize(lo + offsets - 0.5, edges) - 1;                      % number of keys below each query (n_win x n_shuffle), i.e., number of spikes
below_fin = discretize(hi + offsets - 0.5, edges) - 1;                      % of the previous shuffles + spikes of this shuffle before the window

% same as "windowCounts": the search for the window k starts at the first spike found for the window k-1 ("nextStart"); the windows are
% taken shuffle after shuffle, and the keys of a shuffle are always above the ones of the previous shuffles, so "cummax" never mixes them
first = cummax(below_ini(:) + 1);
after = max(first, below_fin(:) + 1);
n_events = reshape(after - first, n_win, n_shuffle);                        % number of spikes of each shuffle within each window
hits = sum(n_events > 0, 1)';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = first_bin(entry, w)

% first bin b (1..n_bin+1) whose spike time, (b-0.5)*bin + ini (same expression as in "ShuffleSpikeTrain"), is >= w
bin = entry.bin;    ini = entry.ini;    n_bin = entry.n_bin;
b = min(max(ceil((w - ini)/bin + 0.5), 1), n_bin + 1);                      % estimate
for it = 1:10                                                               % correct the estimate with the exact expression (rounding errors)
    dec = b > 1 & ((b-1)-0.5)*bin + ini >= w;                               % the previous bin is also >= w
    inc = b <= n_bin & (b-0.5)*bin + ini < w;                               % this bin is still < w
    if ~any(dec) && ~any(inc)
        return
    end
    b(dec) = b(dec) - 1;
    b(inc) = b(inc) + 1;
end
error('shuffleHits: the first bin of a window could not be found');

end
