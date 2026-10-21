function tailPopulationPeaks = tailPopulationPeaks(units, spikeTimestamp, PopulationActivityDetection, bin,  range, opt_norm, opt_secondary)

% function to measure the tail 
% take the binned population firing rate, measure the AUC right after the
% end of a population peak and normalized by the AUC during the peak

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   structure array spikeTimestamp (its length is the number of recording channels), with the following fields
% 	.recLength:	length of the recording, in seconds.
% 	.nUnits:	total number of identified units in the channel
% 	.timeStamp:	time stamp for all units (nUnits) in the channel, in msec

%   units(j):   =0, if channel j should be ignored;     =1, if channel j should be incoporated for population activity
%   bin:       original size of the temporal bin to calculate the population firing rate, in sec        
%   opt_norm:   =0, No normalization; =1: normalization of the firing rate

%   PopulationActivityDetection.IndivBurst(k).peak:     temporal location of the peak of event k, in sec 
%   PopulationActivityDetection.IndivBurst(k).start:    beginning of the event k, in sec
%   PopulationActivityDetection.IndivBurst(k).last:     end of the event k, in sec 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tailPopulationPeaks(k,1):       normalized tail of the peak k (measured using individual windows)     
% tailPopulationPeaks(k,2):       normalized tail of the peak k (measured using a common window for all peaks)      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                            % factor to rescale from msec to sec

WINDOW_FACTOR = 2.0;                            % in terms of the width of the peak (default: 2.0)
FACTOR_SECONDARY = 0.5;                         % factor to rescale WINDOW_FACTOR when analyzing secondary peaks (default: 0.5)

MIN_SAMPLES = 3;
apply_parameters('tailPopulationPeaks');                                    % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Prepare the data 

spikeData = spikeTimestamp;
for k=size(spikeData,2):-1:1
    if ~units(k)
        spikeData(k)= [];
    end
end
n_clusters = size(spikeData,2);                                            % total number of clusters

if opt_secondary
    window = WINDOW_FACTOR*FACTOR_SECONDARY;
else
    window = WINDOW_FACTOR;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Population Firing rate 

aux = bin_spike_time(spikeData(1).timeStamp, bin/K_MSEC_SEC, [range(1) range(2)]);
tail_indiv = nan(numel(PopulationActivityDetection.IndivBurst), 2);     % (revised October 2026) NaN if there are too few bins (before, the output was
                                                                        %  not assigned and the function crashed)

if numel(aux.binMidPoint)>MIN_SAMPLES

time_binned = aux.binMidPoint*K_MSEC_SEC;                                       % time array with the midpoint for each bin, in sec
spikes_indiv = nan(numel(aux.binMidPoint),n_clusters);                          % matrix to save the binned data
for i=1:n_clusters                                                              % across all clusters
    aux = bin_spike_time(spikeData(i).timeStamp, bin/K_MSEC_SEC, [range(1)  range(2)]);    
    if opt_norm
        spikes_indiv(:,i) = aux.normCount;                                      % normalized spike count
    else
        spikes_indiv(:,i) = aux.count;                                          % total spike count per bin
    end
    clear aux
end
spikes = mean(spikes_indiv,2, 'omitnan');                                       % mean Population spike count

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Measure the tail of the peaks, using individual windows

n_peaks = numel(PopulationActivityDetection.IndivBurst);
peakWidth = nan(n_peaks,1);
tail_indiv = nan(n_peaks,2);

for i=1:n_peaks
    peakWidth(i) = PopulationActivityDetection.IndivBurst(i).last - PopulationActivityDetection.IndivBurst(i).start;
    tailStart = PopulationActivityDetection.IndivBurst(i).last;
    tailEnd = min(tailStart + window * peakWidth(i), range(2));
    peakMask = time_binned >= PopulationActivityDetection.IndivBurst(i).start & time_binned <= PopulationActivityDetection.IndivBurst(i).last;
    tailMask = time_binned > tailStart & time_binned <= tailEnd;
    
    tail_rate = trapz(time_binned(tailMask), spikes(tailMask),1) / max(tailEnd - tailStart, eps);   % mean activity during the tail, per unit of time
    peak_rate = trapz(time_binned(peakMask), spikes(peakMask),1) / peakWidth(i);                      % mean activity during the peak, per unit of time
    if peak_rate > 0, tail_indiv(i, 1) = tail_rate / peak_rate; end     % (revised October 2026: NaN instead of Inf/NaN when there is no activity in the peak)
                                                                            % (revised October 2026) ratio of activity per unit of time; before, it was
                                                                            % AUC_tail/(AUC_peak/window) = window^2 times this ratio (x4 for primary, x1 for secondary peaks)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Measure the tail of the peaks, using a common window for all peaks

mean_peakwidth = mean(peakWidth, 'omitnan');

for i=1:n_peaks
    tailStart = PopulationActivityDetection.IndivBurst(i).last;
    tailEnd = min(tailStart + window * mean_peakwidth, range(2));
    peakMask = time_binned >= PopulationActivityDetection.IndivBurst(i).start & time_binned <= PopulationActivityDetection.IndivBurst(i).last;
    tailMask = time_binned > tailStart & time_binned <= tailEnd;
    
    tail_rate = trapz(time_binned(tailMask), spikes(tailMask),1) / max(tailEnd - tailStart, eps);   % mean activity during the tail, per unit of time
    peak_rate = trapz(time_binned(peakMask), spikes(peakMask),1) / peakWidth(i);                      % mean activity during the peak, per unit of time
    if peak_rate > 0, tail_indiv(i, 2) = tail_rate / peak_rate; end     % (revised October 2026: NaN instead of Inf/NaN when there is no activity in the peak)
end

end

tailPopulationPeaks = tail_indiv;                                       % (revised October 2026: after the "if", see above)
