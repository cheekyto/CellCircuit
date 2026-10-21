function [Tfrac, Pfrac] = sttc_one_direction(refTrain, otherTrain, dt, tStart, tEnd, T)

% Function to calculate the Spike Train Tiling Coefficient between two spike trains, in one specfic direction
% (revised October 2026: vectorized, without loops across the spikes -it uses "count_below"-; and fixed the search of coincident spikes, which
%  counted as coincident every spike of refTrain after the last spike of otherTrain -the search index got stuck at the last spike-)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% refTrain:         spikes whose TA / PA we compute
% otherTrain:       spikes used for coincidence / tiling relative to refTrain
% dt:               coincidence window (±dt)
% tStart, tEnd:     recording interval [tStart, tEnd]
% T:                total length of the recording

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tfrac:            fraction of total time within ±dt of ANY spike in refTrain
% Pfrac:            fraction of spikes in refTrain that have at least one spike in otherTrain within ±dt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

refTrain = sort(refTrain(:));       otherTrain = otherTrain(:);
nRef = numel(refTrain);

% ----- Pfrac: fraction of spikes in refTrain that have at least one spike in otherTrain within ±dt -----
if isempty(otherTrain) || nRef == 0
    Pfrac = 0;
else
    n_close = count_below(otherTrain, refTrain + dt, 1) - count_below(otherTrain, refTrain - dt, 0);   % number of spikes of otherTrain in [t-dt, t+dt], for each spike t of refTrain
    Pfrac = sum(n_close > 0) / nRef;
end

% ----- Tfrac: fraction of total time within ±dt of ANY spike in refTrain -----
if nRef == 0
    Tfrac = 0;
    return;
end
starts = max(refTrain - dt, tStart);                                        % [start,end] intervals around each spike, clipped to [tStart,tEnd]
ends   = min(refTrain + dt, tEnd);                                          % (sorted, like the spikes, since all intervals have the same width before clipping)
covered = min(ends(1:end-1), starts(2:end)) - starts(1:end-1);              % part of each interval that is not covered by the next one (merging the overlapping intervals)
total = sum(max(covered, 0)) + (ends(end) - starts(end));                   % total length of the union of all intervals
Tfrac = total / T;

end
