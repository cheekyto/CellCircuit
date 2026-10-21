function bin_spike_times = bin_spike_time(spike_data, binSize, tRange)

% Function to convert an array of spike time to a bin-count vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   spike_data:     numeric vector of spike time, in msec   
%   binSize:        bin width, in msec
%   tRange:         [t_start, t_end], in sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bin_spike_times.count(i):                     % total number of spikes in bin i                    
% bin_spike_times.normCount(i):                 % normalized number of spikes in bin i
% bin_spike_times.binStart(i):                  % start of the bin i, in msec
% bin_spike_times.binMidPoint(i):               % midpoint of the bin i, in msec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                            % factor to rescale from msec to sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (tRange(2) - tRange(1))/K_MSEC_SEC > binSize 
    edges = tRange(1)/K_MSEC_SEC : binSize : tRange(2)/K_MSEC_SEC;          % edges for the bins, in msec
else
    edges = [tRange(1)/K_MSEC_SEC (tRange(1)/K_MSEC_SEC+binSize)];
end
nBins = numel(edges) - 1;                                                   % total number of bins
midpoint = edges + binSize/2.0;                                             % midpoint for each bin, in msec
midpoint(end) = [];

if nBins ~= numel(midpoint)
    error('Number of bins does not match the number of midpoints!!!');
end
bin_spike_times.binMidPoint = midpoint;

if numel(spike_data)>0
    counts = histcounts(spike_data, edges);
    bin_spike_times.normCount = counts/numel(spike_data);
else
    bin_spike_times.normCount = nan;
    counts = midpoint; counts(:) = 0;
end
edges(end) = [];
bin_spike_times.count = counts;
bin_spike_times.binStart = edges;

end