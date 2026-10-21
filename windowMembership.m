function [nEvents, index, nextStart] = windowMembership(data, ini, fin, start, option)

% Events of the sorted vector data0 that fall inside the window [ini, fin).
% The search begins at index "start". nextStart is the index of the first
% element >= ini, so it can be passed as "start" for the next window, as long
% as the windows are processed in increasing order of "ini" (they may overlap).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   data:                       % array with timeStamp for a spike train
%   ini:                        % left side of the temporal window, in the same unit that data (e.g., msec)
%   fin:                        % right side of the temporal window, in the same unit that data (e.g., msec)
%   start:                      % index within data to initiate the search
%   option = 1                  % retrieve the timeStamp of the spikes and their indexes
%                               % only retrieve the number of events 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   nEvents:                    % total number of spikes within the window [ini, fin)
%   timeStamp:                  % timeStamp of the spikes within the window [ini, fin); only when option==1
%   index:                      % indexes of the spikes within the window [ini, fin); only when option==1
%   nextStart:                  % index to start the search of the next window (first element >= ini)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = numel(data);

% first index >= start with data0(i) >= ini   (binary search)
lo = start;
hi = n + 1;
while lo < hi
    mid = floor((lo + hi) / 2);
    if data(mid) >= ini
        hi = mid;
    else
        lo = mid + 1;
    end
end
first = lo;

% first index >= first with data0(i) >= fin   (binary search)
hi = n + 1;
while lo < hi
    mid = floor((lo + hi) / 2);
    if data(mid) >= fin
        hi = mid;
    else
        lo = mid + 1;
    end
end
nextStart = first;
last = lo - 1;

nEvents = last - first + 1;          % number of units within the window [ini, fin)
if option == 1 && nEvents > 0
%    timeStamp = data(first:last);
    index = first:1:last;
else
%    timeStamp = [];
    index = [];
end

end
