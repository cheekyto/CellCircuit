function PopulationActivityDetection = PopulationActivityDetection(units, spikeTimestamp, range, bin, threshold, smoothWindowAdd, opt_norm, opt_secondary, opt_plot)

% function to detect population activation in a spike train data, using the average population firing rate 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   structure array spikeTimestamp (its length is the number of recording channels), with the following fields
% 	.recLength:	length of the recording, in seconds.
% 	.nUnits:	total number of identified units in the channel
% 	.timeStamp:	time stamp for all units (nUnits) in the channel, in msec

%   units(j):   =0, if channel j should be ignored;     =1, if channel j should be incoporated for population activity
%   range:      range of the recording (range(1)=start; range(2)=end), in sec
%   bin:        size of the temporal bin to calculate the population firing rate, in sec        
%   threshold:  in terms of the standard deviation 
%   opt_norm:   =0, No normalization; =1: normalization of the firing rate
%   opt_plot:   =0, no plot; =1: plot the data and the identified network bursts
%   smoothWindow:       smoothing window length, in number of bins (default: 3; odd for symmetric smoothing)
%   opt_secondary:      =1, detection of secondary peaks; =0, detection of primary peaks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PopulationActivityDetection.spikesBin(i,1):           total number of spikes in bin i
% PopulationActivityDetection.time(i):                  midpoint of the bin i, in sec 
% PopulationActivityDetection.nBurst:                   total number of network bursts 
% PopulationActivityDetection.IndivBurst(j).peak:       peak point of network burst j, in sec (after refinement)       
% PopulationActivityDetection.IndivBurst(j).start:      start point of network burst j, in sec (after refinement)
% PopulationActivityDetection.IndivBurst(j).last:       last point of network burst j, in sec (after refinement)
% PopulationActivityDetection.IndivBurst(j).peakIdx:    bin index of the peak point of network burst j
% PopulationActivityDetection.IndivBurst(j).tail:       normalized tail of the network burst j

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

START_LEVEL = 0.20;                             % proportion of the peak amplitude of the network burst to identifity the starting point (default: 0.10)
END_LEVEL = 0.20;                               % proportion of the peak amplitude of the network burst to identifity the ending point (default: 0.10)
K_MSEC_SEC = 1/1000;                            % factor to rescale from msec to sec

SMOOTH_FACTOR = [3; 5; 7];                      % standard smoothing window length, in number of bins (default: 3; odd for symmetric smoothing)
MIN_POINTS = 3;

OUTLIER_TOP = 99.5;                             % upper percentile to identify outliers, before calculating the noise (dafault: 99)

PERCENTILE_TOP = 95.0;                          % upper percentile to identify the "big mountains" (default: 95)
CHAIN = 5;                                      % total number of consecutive high-value to identify "big mountains" (default: 3)
FACTOR_MOUNTAIN1 = 2.0;                         % (default: 1.8) (in 100 bin and 5 high-values, the chances of having 3 or more consecutives high-values is 1/170)
FACTOR_MOUNTAIN2 = 2.5;                         % (default: 2.5) (in 100 bin and 5 high-values, the chances of having 3 or more consecutives high-values is 1/170)
BIN_REDUCT = 0.30;                              % (default: 0.30)
MIN_THRESHOLD = 3.0;                            % (default: 3.0)

ACTIVITY_THRESHOLD = 0.03;                      % minimum level of spiking frequency for an individual electrode to be considered active, in Hz (default: 0.05)
                                                % (revised October 2026: 0.05 -> 0.03 Hz, user's choice)

PEAK_THRESHOLD = 20;                            % threshold to run the second round of searching, relative to median firing rate (default: 20)
apply_parameters('PopulationActivityDetection');                            % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Population Firing rate 

spikeData = spikeTimestamp;
for k=size(spikeData,2):-1:1
    if ~units(k)
        spikeData(k)= [];
    else
        if ~opt_secondary && spikeData(k).nUnits/(range(2)-range(1)) < ACTIVITY_THRESHOLD
            spikeData(k)= [];
            units(k) = 0;                                                       % (revised October 2026) the cluster is also removed from "units", so the refinement 
                                                                                % of the peaks, the size check and the tail (below) use the same clusters as the
                                                                                % detection (before, they also used the clusters below ACTIVITY_THRESHOLD)
        end
    end
end
n_clusters = size(spikeData,2);                                                 % total number of clusters

if n_clusters > 0

aux = bin_spike_time(spikeData(1).timeStamp, bin/K_MSEC_SEC, range);
time_binned = aux.binMidPoint*K_MSEC_SEC;                                       % time array with the midpoint for each bin, in sec
spikes_indiv = nan(numel(aux.binMidPoint),n_clusters);                          % matrix to save the binned data
for i=1:n_clusters                                                              % across all clusters
    aux = bin_spike_time(spikeData(i).timeStamp, bin/K_MSEC_SEC, range);    
    if opt_norm
        spikes_indiv(:,i) = aux.normCount;                                      % normalized spike count
    else
        spikes_indiv(:,i) = aux.count;                                          % total spike count per bin
    end
    clear aux
end
spikes = mean(spikes_indiv,2, 'omitnan');                                       % mean Population spike count
PopulationActivityDetection.spikesBin = spikes;                                 % mean Population spike count (normalized or nominal, depending on "opt_norm"
PopulationActivityDetection.time = time_binned;                                 % array with the midpoint of the bins, in sec
PopulationActivityDetection.bin = bin;                                          % bin size

if opt_secondary
    sm_factor = SMOOTH_FACTOR(1);
else
    sm_factor = SMOOTH_FACTOR;
end    
n_smoothing = numel(sm_factor);
n_peaks_smooth = nan(n_smoothing,1);
n_peaksLevel1_smooth = nan(n_smoothing,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Check if there are big mountains that are driving the noise

threshold0 = threshold;
if ~opt_secondary
    signal = spikes;
    p95 = prctile(signal, PERCENTILE_TOP);
    isHigh = signal >= p95;
    d = diff([false; isHigh; false]);      % pad with false at both ends
    startIdx = find(d == 1);            % run starts
    endIdx   = find(d == -1) - 1;       % run ends
    runLengths = endIdx - startIdx + 1;
    validRuns = find(runLengths >= CHAIN);

    idxAllInRuns = [];
    for k = validRuns'
        idxAllInRuns = [idxAllInRuns, startIdx(k):endIdx(k)];
    end
    numHighInRuns = numel(idxAllInRuns);
    mountain = 100*numHighInRuns/numel(signal);
    if mountain > FACTOR_MOUNTAIN1
        threshold = threshold * (1-BIN_REDUCT*(1-exp(-1*mountain/FACTOR_MOUNTAIN2)));
    end
end

for sm=1:n_smoothing

smoothWindow = sm_factor(sm) + smoothWindowAdd;
smoothedSpikes0 = smoothdata(spikes, 'gaussian', smoothWindow);                  % Gaussian smoothing
smoothedSpikes = filloutliers(smoothedSpikes0(:,1),"clip","percentiles",[0 OUTLIER_TOP]);
smoothed_spikesBin{sm} = smoothedSpikes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Identify the supra-threshold time points

threshold = max(threshold, MIN_THRESHOLD);
threshold_peak = std(smoothedSpikes)*threshold;                                % threshold to identify the peak of potential network burst events
if numel(smoothedSpikes) >= MIN_POINTS
    [pks, locs, widths, proms] = findpeaks(smoothedSpikes, MinPeakProminence=threshold_peak);
else
    pks = [];
end
contador = numel(pks);
n_peaksLevel1_smooth(sm) = contador;

if contador>0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Identifity the starting and ending point of each potential network burst 

    netBurstLimits = QuickNetBurstLimits(pks, locs, proms, smoothedSpikes, START_LEVEL, END_LEVEL);
    start_Nburst = netBurstLimits.start;
    end_Nburst = netBurstLimits.end;

    locs_final = locs;
    proms_final = proms;
    start_Nburst_final = start_Nburst;
    end_Nburst_final = end_Nburst;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Remove the strikingly large peaks and recalculate the threshold for peak detection 
    
    tamanio = mean(pks)/median(smoothedSpikes);
    smoothedSpikes_cropped = smoothedSpikes;

    if tamanio > PEAK_THRESHOLD && ~opt_secondary
        for j=1:contador
            ini = start_Nburst(j);
            fin = min(size(smoothedSpikes,1), end_Nburst(j) + end_Nburst(j) - ini);
            smoothedSpikes_cropped(ini:fin) = nan;
        end
%        threshold0 = max(threshold0, MIN_THRESHOLD);
        threshold0 = max(threshold, MIN_THRESHOLD);
        threshold_peak_low = std(smoothedSpikes_cropped,'omitnan')*threshold0;
        threshold_second = min(threshold_peak_low, threshold_peak);

        %fprintf("Running 2nd round of Population Peak Detection with an absolute threshold reduced by a factor of %f \n", threshold_peak/threshold_second)
    
        [pks2, locs2, widths2, proms2] = findpeaks(smoothedSpikes, MinPeakProminence=threshold_second);

        netBurstLimits = QuickNetBurstLimits(pks2, locs2, proms2, smoothedSpikes, START_LEVEL, END_LEVEL);
        start_Nburst2 = netBurstLimits.start;
        end_Nburst2 = netBurstLimits.end;

        locs_final = locs2;
        proms_final = proms2;
        start_Nburst_final = start_Nburst2;
        end_Nburst_final = end_Nburst2;                
    end
            
    peak_level = zeros(numel(locs_final),1);
    peak_level(find(proms_final > threshold_peak)) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Integrate the events that overlap and generate the final list of network bursts

    peak_final = locs_final;
    contador_final = numel(peak_final);
    start_burst = time_binned(start_Nburst_final)+bin/2;
    end_burst = time_binned(end_Nburst_final)-bin/2;
    
    flagged = zeros(contador_final,1);
    if contador_final>1
        for k=2:contador_final
           if start_Nburst_final(k) <= end_Nburst_final(k-1)                % there is overlap between consecutive events
               if peak_level(k) < peak_level(k-1)                           % only when the first event is higher level than the second one
                   flagged(k) = 1;                                          % flagged the second event (k) to be discarded
               else                                                         % if they have the same level or if the second event is higher level than the first one
                   flagged(k-1) = 1;                                        % flagged the first event (k-1) to be discarded
               end
           end
        end
        if sum(flagged) > 0                                                 % at least one of the events were flagged to be discarded
            for k=contador_final:-1:1
                if flagged(k)                                               % remove the events that were flagged
                    peak_final(k) = []; 
                    start_burst(k) = [];
                    end_burst(k) = [];
                    peak_level(k) = [];
                end
            end
        end
    end

    IndBurst(sm).peak = peak_final;
    IndBurst(sm).start = start_burst;
    IndBurst(sm).end = end_burst;
    IndBurst(sm).level = peak_level;
    n_peaks_smooth(sm) = numel(peak_final); 
    clear peak_final start_burst end_burst
end
end % end for the for loop across smoothing factors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Look for the best smoothing factor

[n_peaks_final, sm_optimal] = max(n_peaks_smooth);
PopulationActivityDetection.smoothed_spikesBin = smoothed_spikesBin{sm_optimal};                % smoothed mean Population spike count

if n_peaks_final>0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Count only the peaks level 1

    count_level1 = 0; 
    for k=1:n_peaks_final 
        if IndBurst(sm_optimal).level(k)
            count_level1 = count_level1 + 1;
        end
    end
    if count_level1 > n_peaksLevel1_smooth(sm_optimal)
        error('\n There was a problem in PopulationActivityDetection: a mismatch in the number of level1 peaks that were identified \n')
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Refine the characterization of the peaks
    
    peak=nan(n_peaks_final,1);    peakAmp=nan(n_peaks_final,1);     peakIdx=nan(n_peaks_final,1);   start=nan(n_peaks_final,1);
    last=nan(n_peaks_final,1);    level=nan(n_peaks_final,1);
    spike_merged = mergeSpikeTrains(units, spikeTimestamp);                 % (revised October 2026) spike trains of the selected clusters, merged once for all the peaks 
    for k=1:n_peaks_final       
        peakCharact = refinePeakDetection(units, spikeTimestamp, IndBurst(sm_optimal).start(k), IndBurst(sm_optimal).end(k), bin, START_LEVEL, END_LEVEL, opt_norm, opt_plot, spike_merged);
        peak(k) = peakCharact.peak;                                 % temporal location of the peak of event k, in sec 
        peakAmp(k) = peakCharact.peakAmp;                           % Amplitude of the peak of event k 
        peakIdx(k) = IndBurst(sm_optimal).peak(k);                  % location of the peak of the event k
        start(k) = peakCharact.start;                               % beginning of the event k, in sec
        last(k) = peakCharact.end;                                  % end of the event k, in sec 
        level(k) = IndBurst(sm_optimal).level(k);                   % level of the event k (=0, small event; =1, large event)  
    end
    label_size = PopulationActivitySize(units, spikeTimestamp, start/K_MSEC_SEC, last/K_MSEC_SEC, 1);

%    PopulationActivityDetection.nBurst = n_peaks_final;
    PopulationActivityDetection.nBurst = nnz(label_size);
    idx = 0;
    for k=1:n_peaks_final       
        if label_size(k)
            idx = idx + 1;
            PopulationActivityDetection.IndivBurst(idx).peak = peak(k);             % temporal location of the peak of event k, in sec 
            PopulationActivityDetection.IndivBurst(idx).peakAmp = peakAmp(k);       % Amplitude of the peak of event k 
            PopulationActivityDetection.IndivBurst(idx).peakIdx = peakIdx(k);       % location of the peak of the event k
            PopulationActivityDetection.IndivBurst(idx).start = start(k);           % beginning of the event k, in sec
            PopulationActivityDetection.IndivBurst(idx).last = last(k);             % end of the event k, in sec 
            PopulationActivityDetection.IndivBurst(idx).level = level(k);           % level of the event k (=0, small event; =1, large event)  
       else
            if level(k) == 1
                count_level1 = count_level1 - 1;
            end
        end
    end
    PopulationActivityDetection.nBurst_level1 = count_level1;
    if PopulationActivityDetection.nBurst ~= idx
        error('\n \n There was a mismatch while removing population peaks and counting the final number of peaks (in PopulationActivityDetection) \n \n')
    end
    if PopulationActivityDetection.nBurst > 0
        tailPopulationPeak = tailPopulationPeaks(units, spikeTimestamp, PopulationActivityDetection, bin,  range, opt_norm, opt_secondary);
        for k=1:PopulationActivityDetection.nBurst
            PopulationActivityDetection.IndivBurst(k).tail = tailPopulationPeak(k,1);
        end
    else
        PopulationActivityDetection.IndivBurst.tail = [];                       
    end

else
    PopulationActivityDetection.nBurst = 0;
    PopulationActivityDetection.nBurst_level1 = 0;
    PopulationActivityDetection.IndivBurst.peak = [];                       % temporal location of the peak of event k, in sec (center of the peak bin)
    PopulationActivityDetection.IndivBurst.peakIdx = [];                    % location of the peak of the event k
    PopulationActivityDetection.IndivBurst.start = [];                      % beginning of the event k, in sec
    PopulationActivityDetection.IndivBurst.last = [];                       % end of the event k, in sec 
    PopulationActivityDetection.IndivBurst.tail = [];                       % normalized tail of the event k
    PopulationActivityDetection.IndivBurst.peakAmp = [];
    PopulationActivityDetection.IndivBurst.level = [];
end

else
    PopulationActivityDetection.spikesBin = [];                             % (revised October 2026) the same fields as when there are clusters (no data)
    PopulationActivityDetection.time = [];
    PopulationActivityDetection.bin = bin;
    PopulationActivityDetection.smoothed_spikesBin = [];
    PopulationActivityDetection.nBurst = 0;
    PopulationActivityDetection.nBurst_level1 = 0;
    PopulationActivityDetection.IndivBurst.peak = [];                       % temporal location of the peak of event k, in sec (center of the peak bin)
    PopulationActivityDetection.IndivBurst.peakIdx = [];                    % location of the peak of the event k
    PopulationActivityDetection.IndivBurst.start = [];                      % beginning of the event k, in sec
    PopulationActivityDetection.IndivBurst.last = [];                       % end of the event k, in sec 
    PopulationActivityDetection.IndivBurst.tail = [];                       % normalized tail of the event k
    PopulationActivityDetection.IndivBurst.peakAmp = [];
    PopulationActivityDetection.IndivBurst.level = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     options to plot the process

if opt_plot && n_clusters > 0                                               % (revised October 2026: nothing to plot without clusters; before, it crashed)
    figure()
    bar(time_binned, spikes); hold on;
    plot(time_binned, smoothed_spikesBin{sm_optimal},'m--' ); hold on;
    if PopulationActivityDetection.nBurst > 0
        for i=1:PopulationActivityDetection.nBurst
            plot(PopulationActivityDetection.IndivBurst(i).peak, PopulationActivityDetection.IndivBurst(i).peakAmp*1.1,'ro', 'MarkerFaceColor','r'); hold on
            plot(PopulationActivityDetection.IndivBurst(i).start, PopulationActivityDetection.IndivBurst(i).peakAmp/2,'k*'); hold on
            plot(PopulationActivityDetection.IndivBurst(i).last, PopulationActivityDetection.IndivBurst(i).peakAmp/2,'g*'); hold on
        end
    end
end
end

