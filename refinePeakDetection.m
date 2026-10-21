function refinePeakDetection = refinePeakDetection(units, spikeTimestamp, ini0, fin0, bin0, start_level, end_level, opt_norm, opt_plot, spike_merged)

test = 0; 

% function to refine the characterization of detected population activation in a spike train data, using the average population firing rate 
% it analyses single peaks at a time
% (revised October 2026: all the clusters are binned at once, from the merged spike trains -see "mergeSpikeTrains"-, instead of one call of 
%  "bin_spike_time" per cluster and bin size, followed by a "clear"; same results)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   structure array spikeTimestamp (its length is the number of recording channels), with the following fields
% 	.recLength:	length of the recording, in seconds.
% 	.nUnits:	total number of identified units in the channel
% 	.timeStamp:	time stamp for all units (nUnits) in the channel, in msec

%   units(j):   =0, if channel j should be ignored;     =1, if channel j should be incoporated for population activity
%   ini0:       start of the population spike, in sec
%   fin0:       end of the population spike, in sec
%   bin0:       original size of the temporal bin to calculate the population firing rate, in sec        
%   start_level:proportion of the peak amplitude of the network burst to identifity the starting point    
%   end_level:  proportion of the peak amplitude of the network burst to identifity the ending point
%   opt_norm:   =0, No normalization; =1: normalization of the firing rate
%   spike_merged:   (optional) output of "mergeSpikeTrains(units, spikeTimestamp)"; PopulationActivityDetection builds it once for all the peaks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% refinePeakDetection.peak:       peak point of network burst, in sec       
% refinePeakDetection.peakAmp:   Amplitude of the peak of network burst      
% refinePeakDetection.start:      start point of network burst, in sec
% refinePeakDetection.end:        last point of network burst, in sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                            % factor to rescale from msec to sec

smoothWindow = 3;                               % smoothing window length, in number of bins (default: 3; odd for symmetric smoothing)
WINDOW_FACTOR = 0.5;                            % in terms of the width of the peak (default: 0.5)

BIN_FACTOR = [0.25; 0.40; 0.70; 1.0; 1.20];     % in terms of the original bin (default: 0.4, 0.7, 1.0)

MIN_SAMPLES = 3;
apply_parameters('refinePeakDetection');                                    % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Prepare the data 

if nargin < 10 || isempty(spike_merged)                                    % (revised October 2026) spike trains of the selected clusters, merged and sorted
    spike_merged = mergeSpikeTrains(units, spikeTimestamp);
end
n_clusters = spike_merged.n_clusters;                                       % total number of clusters

width = fin0-ini0;
ini = ini0-width*WINDOW_FACTOR;
fin = fin0+width*WINDOW_FACTOR;

peakTime = nan(numel(BIN_FACTOR),1);
startTime = nan(numel(BIN_FACTOR),1);
endTime = nan(numel(BIN_FACTOR),1); 
peakAmp = nan;

if test && opt_plot
    figure();
end

for k=1:numel(BIN_FACTOR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Population Firing rate 

bin = bin0*BIN_FACTOR(k);
binSize = bin/K_MSEC_SEC;                                                       % bin size, in msec
% (revised October 2026) the same bins as "bin_spike_time(spikes, bin/K_MSEC_SEC, [ini fin])" (same expressions)
if (fin - ini)/K_MSEC_SEC > binSize 
    edges = ini/K_MSEC_SEC : binSize : fin/K_MSEC_SEC;                          % edges for the bins, in msec
else
    edges = [ini/K_MSEC_SEC (ini/K_MSEC_SEC+binSize)];
end
binMidPoint = edges + binSize/2.0;                                              % midpoint for each bin, in msec
binMidPoint(end) = [];

if numel(binMidPoint)>MIN_SAMPLES

time_binned = binMidPoint*K_MSEC_SEC;                                           % time array with the midpoint for each bin, in sec
% (revised October 2026) all the clusters at once: the spikes within [edges(1), edges(end)] (as "histcounts", the last bin includes its right edge)
first_spike = first_at_or_above(spike_merged.time, edges(1));                   % first spike >= edges(1)
last_spike = first_above(spike_merged.time, edges(end)) - 1;                    % last spike <= edges(end)
bin_idx = discretize(spike_merged.time(first_spike:last_spike), edges);         % bin of each spike (same bin rules as "histcounts")
counts = accumarray([bin_idx(:), spike_merged.unit(first_spike:last_spike)], 1, [numel(binMidPoint), n_clusters]);  % spike count per bin (row) and cluster (column)
if opt_norm
    spikes_indiv = counts ./ spike_merged.nSpikes(:)';                          % normalized spike count (a cluster without spikes gives NaN, as "bin_spike_time")
else
    spikes_indiv = counts;                                                      % total spike count per bin
end
spikes = mean(spikes_indiv,2, 'omitnan');                                       % mean Population spike count
smoothedSpikes = smoothdata(spikes, 'gaussian', smoothWindow);                  % Gaussian smoothing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Identify the supra-threshold time points

[pks, locs, widths, proms] = findpeaks(smoothedSpikes);

if numel(pks)>0
    if numel(pks)>1
       [~, aux] = max(pks);
       idx_peak = locs(aux);
       peak = pks(aux);
       amplitude = peak;                                                    % (revised October 2026) the height of the peak, as when there is a single peak 
                                                                            % (before, its prominence: "amplitude = proms(aux)"; two different criteria)
    else
        idx_peak = locs;
        peak = pks;
        amplitude = proms; 
        amplitude = peak; 
    end
else
    [peak, idx_peak] = max(smoothedSpikes);
    amplitude = peak;
end

peakTime(k) = time_binned(idx_peak);
if BIN_FACTOR(k)==1
%    peakAmp = amplitude;
    peakAmp = peak;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Identifity the starting and ending point of each potential network burst 

threshold_ini = peak-(1-start_level)*amplitude;                             % threshold to identify the start of the burst event k    
threshold_fin = peak-(1-end_level)*amplitude;                               % threshold to identify the end of the burst event k
            
aux = find(smoothedSpikes(1:idx_peak,1)<threshold_ini);                     % index of all bins below the threshold
[~, idx_min] = min(smoothedSpikes(1:idx_peak,1));                           % lowest point in between the peak k and the beginning
shift = 0;                                                                  % to adjust the index
if numel(aux)>0                                                             % at least one point crossing the threshold level
   start_Nburst = aux(end)+shift;                                           % first point below the threshold, to the left of the peak k
   % (revised October 2026) start = time at which the smoothed signal crosses the threshold, by linear interpolation between the last bin 
   % below the threshold and the next one (before, the midpoint of the bin below the threshold, which depends on the bin size)
   startTime(k) = time_binned(start_Nburst) + (threshold_ini - smoothedSpikes(start_Nburst)) / ...
                  (smoothedSpikes(start_Nburst+1) - smoothedSpikes(start_Nburst)) * (time_binned(start_Nburst+1) - time_binned(start_Nburst));
else
   start_Nburst = idx_min+shift;                                            % lowest point to the left of the peak k        
   startTime(k) = time_binned(start_Nburst);                                % (no crossing: midpoint of the lowest bin, as before)
end
clear aux shift idx_min
        
aux = find(smoothedSpikes(idx_peak:end,1)<threshold_fin);                   % index of all bins below the threshold    
[~, idx_min] = min(smoothedSpikes(idx_peak:numel(smoothedSpikes)));         % lowest point in between the peak k and the end
shift = idx_peak-1;
if numel(aux)>0                                                             % at least one point crossing the threshold level
    end_Nburst = aux(1)+shift;                                              % first point below the threshold, to the right of the peak
    % (revised October 2026) end = threshold crossing, by linear interpolation between the last bin above the threshold and the first one below it
    endTime(k) = time_binned(end_Nburst-1) + (smoothedSpikes(end_Nburst-1) - threshold_fin) / ...
                 (smoothedSpikes(end_Nburst-1) - smoothedSpikes(end_Nburst)) * (time_binned(end_Nburst) - time_binned(end_Nburst-1));
else
    end_Nburst = idx_min + shift;                                           % lowest point to the right of the peak k
    endTime(k) = time_binned(end_Nburst);                                   % (no crossing: midpoint of the lowest bin, as before)
end
clear aux shift idx_min

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Plot during the code developing phase

if test && opt_plot   
   subplot(numel(BIN_FACTOR),1,k)
   bar(time_binned, spikes); hold on;
   plot(time_binned, smoothedSpikes,'m--' ); hold on;
   plot(peakTime(k), smoothedSpikes(idx_peak)*1.1,'ro', 'MarkerFaceColor','r'); hold on
   plot(startTime(k), smoothedSpikes(idx_peak)*start_level*1.1,'k*'); hold on
   plot(endTime(k), smoothedSpikes(idx_peak)*end_level*1.1,'g*'); hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end % end of the for-loop across binning factors

refinePeakDetection.peak = median(peakTime, 'omitnan');
refinePeakDetection.peakAmp = peakAmp;
refinePeakDetection.start = median(startTime, 'omitnan');
refinePeakDetection.end = median(endTime, 'omitnan');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idx = first_at_or_above(t, x)

% binary search in the sorted array t: index of the first element >= x (numel(t)+1 if there is none) (revised October 2026)
lo = 1;     hi = numel(t) + 1;
while lo < hi
    mid = floor((lo + hi)/2);
    if t(mid) < x
        lo = mid + 1;
    else
        hi = mid;
    end
end
idx = lo;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idx = first_above(t, x)

% binary search in the sorted array t: index of the first element > x (numel(t)+1 if there is none) (revised October 2026)
lo = 1;     hi = numel(t) + 1;
while lo < hi
    mid = floor((lo + hi)/2);
    if t(mid) <= x
        lo = mid + 1;
    else
        hi = mid;
    end
end
idx = lo;

end
