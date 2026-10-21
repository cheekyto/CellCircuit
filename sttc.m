function sttc = sttc(spikeTimes1, spikeTimes2, dt, tRange)

% Function to calculate the Spike Train Tiling Coefficient between two spike trains.
% (revised October 2026: the recording interval can be provided -tRange-; before, it was always taken from the first and the last spike
%  of the two trains, so the tiling fractions depended on when the spike trains happened to start and end)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   spikeTimes1, spikeTimes2 :  spike times (vectors, same units)
%   dt                       :  coincidence window (±dt)
%   tRange                   :  (optional) [tStart tEnd] recording interval, same units as the spike times
%                               if it is not provided, it is taken from the first and the last spike of the two trains (as before)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   sttc:    Spike Train Tiling Coefficient (range approx -1 to 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensure column vectors
xx = spikeTimes1(:);
yy = spikeTimes2(:);

if isempty(xx) && isempty(yy)
    sttc = NaN;                                             % undefined if no spikes in both trains
    return;
end
if nargin >= 4 && ~isempty(tRange)                          % recording interval provided
    tStart = tRange(1);
    tEnd   = tRange(2);
else                                                        % automatically determine tStart / tEnd from data
    allSpikes = [xx; yy];
    tStart = min(allSpikes);
    tEnd   = max(allSpikes);
end

% Restrict spikes to [tStart, tEnd]
xx = xx(xx >= tStart & xx <= tEnd);
yy = yy(yy >= tStart & yy <= tEnd);

T = tEnd - tStart;
if T < 0
   error('tEnd must be > tStart');
else
    % If either train has no spikes, STTC is 0 by definition (Cutts & Eglen)
    if numel(xx) == 0 || numel(yy) == 0
        sttc = 0;
        return;
    end

    % Compute TA, PA (A relative to B)
    [TA, PA] = sttc_one_direction(xx, yy, dt, tStart, tEnd, T);

    % Compute TB, PB (B relative to A)
    [TB, PB] = sttc_one_direction(yy, xx, dt, tStart, tEnd, T);

    % STTC formula
    termA = (PA - TB) / (1 - PA * TB);
    termB = (PB - TA) / (1 - PB * TA);

    sttc = 0.5 * (termA + termB);
end

end
