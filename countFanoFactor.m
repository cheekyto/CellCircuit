function countFanoFactor = countFanoFactor(eventTime, tRange, window)

% function to calculate the Fano factor of the number of events (spikes, bursts or network bursts) counted in consecutive time windows
% (revised October 2026)
% FFcount = variance / mean of the number of events per window, across windows (sample variance, n-1), as in Rajdl et al. 2020
% (Front. Comput. Neurosci. 14:569049, equation 20). The windows are consecutive and non-overlapping, starting at the beginning of the
% analysis window; the last, incomplete window is not used. For a renewal process FFcount tends to CV^2 of the intervals for long windows,
% and it is biased toward 1 when a window has few events (less than ~5-10 on average): "meanCount" helps to judge it

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   eventTime:      time of each event, in msec
%   tRange:         [start end] of the analysis window, in sec
%   window:         length of the counting windows, in sec (e.g., 60)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   countFanoFactor.FFcount:      Fano factor of the counts (NaN with less than 2 windows, or without events)
%   countFanoFactor.meanCount:    mean number of events per window
%   countFanoFactor.nWindows:     number of windows used
%   countFanoFactor.counts(j):    number of events in the window j (column)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                                                        % factor to rescale from msec to sec

n_windows = floor((tRange(2) - tRange(1))/window + 1e-9);                   % number of complete windows (the tolerance avoids losing a window to round-off)
countFanoFactor.FFcount = nan;
countFanoFactor.meanCount = nan;
countFanoFactor.nWindows = max(n_windows, 0);
countFanoFactor.counts = zeros(0,1);
if n_windows < 1
    return
end

edges = tRange(1) + (0:n_windows)'*window;                                  % edges of the windows, in sec
t = eventTime(:)*K_MSEC_SEC;                                                % time of the events, in sec
counts = count_below(t, edges(2:end), 0) - count_below(t, edges(1:end-1), 0);  % number of events within [edges(j), edges(j+1)) of each window j
countFanoFactor.counts = counts;
countFanoFactor.meanCount = mean(counts);
if n_windows >= 2 && countFanoFactor.meanCount > 0
    countFanoFactor.FFcount = var(counts)/mean(counts);                     % variance (n-1) to mean ratio
end

end
