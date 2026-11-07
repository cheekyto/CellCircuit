function event = characterizeEvent(units, spikeTimestamp, time_binned, bin, start_idx, end_idx, best, start_level, end_level, opt_norm, opt_plot, spike_merged)

% (new, November 2026) characterization of one event of options 1-3 of "NetworkBurstDetection" from its candidates (one, or several merged
% because their edges overlap). It is the code of "PopulationActivityDetection", moved here so that the events added by option 5 (large events
% of option 4, see "add_large_option4" in "NetworkBurstDetection") are characterized by exactly the same code: a change here applies to both.
% Each candidate, from its edges in bins (from "QuickNetBurstLimits"), is refined with "refinePeakDetection"; the peak of the event is the one
% of its best candidate, and its edges (20% level) and half-width edges (50% level) are the outermost ones of its candidates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   units(j):           =0, if cluster j should be ignored; =1, if cluster j is used for the population activity
%   spikeTimestamp:     structure array with the field .timeStamp (in msec), one element per cluster
%   time_binned(i):     midpoint of the bin i of the population signal, in sec
%   bin:                size of the bins of the population signal, in sec
%   start_idx, end_idx: edges of each candidate of the event, as bin indices of the population signal (one element per candidate)
%   best:               the candidate that gives the peak of the event
%   start_level, end_level:     proportion of the peak amplitude that defines the start and the end of the event (START_LEVEL, END_LEVEL)
%   opt_norm:           =0, no normalization; =1, normalized firing rate; =2, participation
%   opt_plot:           plots of "refinePeakDetection"
%   spike_merged:       (optional) output of "mergeSpikeTrains(units, spikeTimestamp)", to build it once for many events

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   event.peak, event.peakAmp:      time (sec) and amplitude of the peak (best candidate, refined)
%   event.start, event.end:         edges of the event (20% level), in sec (outermost of its candidates)
%   event.start20, event.end20:     the same
%   event.start50, event.end50:     edges at 50% of the height of the peak (half-width), in sec (outermost of its candidates)
%   (and the other fields of the output of "refinePeakDetection" for the best candidate)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 12
    spike_merged = mergeSpikeTrains(units, spikeTimestamp);
end
cand_start = time_binned(start_idx) + bin/2;                                % edges of the candidates, in sec
cand_end = time_binned(end_idx) - bin/2;
n_cand = numel(start_idx);
s_c = nan(n_cand,1);    e_c = nan(n_cand,1);    s50_c = nan(n_cand,1);  e50_c = nan(n_cand,1);
for c=1:n_cand
    refined = refinePeakDetection(units, spikeTimestamp, cand_start(c), cand_end(c), bin, start_level, end_level, opt_norm, opt_plot, spike_merged);
    s_c(c) = refined.start;     e_c(c) = refined.end;   s50_c(c) = refined.start50;     e50_c(c) = refined.end50;
    if c == best
        event = refined;                                                    % peak of the event: its best candidate
    end
end
event.start = min(s_c);     event.end = max(e_c);                           % edges of the event: the outermost ones (20% level)
event.start20 = event.start;    event.end20 = event.end;
event.start50 = min(s50_c);     event.end50 = max(e50_c);                   % half-width (50% level), reported as a metric

end
