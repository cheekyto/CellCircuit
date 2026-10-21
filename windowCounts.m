function [n_events, first] = windowCounts(data, ini, fin)

% function to count the spikes within several temporal windows at once (revised October 2026)
% Same result as calling "windowMembership" for the windows one after the other, passing "nextStart" from one call to the next
% (i.e., the search for the window k starts at the first spike found for the window k-1), but without a loop across the windows

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data:                 time stamp of the spikes (sorted in ascending order, as required by "windowMembership")
% ini(k), fin(k):       start and end of the window k; a spike t is within the window if ini(k) <= t < fin(k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% n_events(k):          number of spikes within the window k (column array)
% first(k):             index of the first spike within the window k (the spikes are first(k) : first(k)+n_events(k)-1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

first = cummax(count_below(data, ini, 0) + 1);                              % first spike >= ini(k), but never before the first spike found for the previous window ("nextStart")
after = max(first, count_below(data, fin, 0) + 1);                          % first spike >= fin(k), searching from first(k) (as "windowMembership")
n_events = after - first;                                                   % number of spikes in [ini(k), fin(k))

end
