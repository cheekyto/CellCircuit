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
%   threshold:  in terms of the standard deviation (not used with opt_norm = 2, except =0: SD threshold, MIN_THRESHOLD, as in the scan of
%               "Optimal_secondaryPeak")
%   opt_norm:   =0, No normalization; =1: normalization of the firing rate
%               =2 (new, November 2026; option 4 of "NetworkBurstDetection"): participation (fraction of the units with at least one spike in
%               each bin), with a surrogate threshold (see "participationSurrogate"): large events (level 1), peak above the (1-ALPHA_LARGE)
%               quantile of the maxima of the surrogates; small events (level 0), peak above the PERCENTILE_SMALL percentile of all the
%               surrogate bins pooled; secondary events (opt_secondary = 1), peak above the (1-ALPHA_SECONDARY) quantile of the maxima of the
%               surrogates within the window. Only absolute values (no prominence): one candidate per stretch above the threshold, with the
%               edges at 50% of the height of its peak above the median of the signal (and 20%, reported); candidates whose 50% edges
%               overlap are merged
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
% PopulationActivityDetection.IndivBurst(j).pSurrogate: (opt_norm = 2) p-value of the peak, against the maxima of the surrogates
%                                                       (1 + number of surrogate maxima >= peak)/(N_SURROGATE + 1)
% PopulationActivityDetection.surrogate:                (opt_norm = 2) thresholds and settings of the surrogates ([] for the SD scan)

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
TAIL_FACTOR = 3;                                % (new, November 2026) the tail of a large event lasts until the signal returns to its median (baseline),
                                                % and at most TAIL_FACTOR times the duration of the event after its end (default: 3)

% (new, November 2026) option 4 of "NetworkBurstDetection" (opt_norm = 2): participation signal with a surrogate threshold
N_SURROGATE = 1000;                             % number of surrogates; each unit is circularly shifted by an independent random offset (default: 1000)
SURROGATE_SEED = 0;                             % seed of the random offsets (the same seed gives the same surrogates)
ALPHA_LARGE = 0.05;                             % large (level-1) events: peak above the (1-ALPHA_LARGE) quantile of the surrogate maxima (max statistic) (default: 0.05)
PERCENTILE_SMALL = 99;                          % small (level-0) events: peak above this percentile of all the surrogate bins pooled (Bocchio et al. 2020) (default: 99)
ALPHA_SECONDARY = 0.05;                         % secondary events (opt_secondary = 1): peak above the (1-ALPHA_SECONDARY) quantile of the surrogate maxima
                                                % within the window (max statistic) (default: 0.05)
apply_parameters('PopulationActivityDetection');                            % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Population Firing rate 

surrogate_info = [];                                                            % (new, November 2026) thresholds of the surrogates (opt_norm = 2)
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
    if opt_norm == 2
        spikes_indiv(:,i) = aux.count > 0;                                      % (new, November 2026) participation: the cluster fires in the bin
    elseif opt_norm
        spikes_indiv(:,i) = aux.normCount;                                      % normalized spike count
    else
        spikes_indiv(:,i) = aux.count;                                          % total spike count per bin
    end
    clear aux
end
spikes = mean(spikes_indiv,2, 'omitnan');                                       % mean Population spike count
PopulationActivityDetection.spikesBin = spikes;                                 % mean Population spike count (normalized or nominal, or participation, depending on "opt_norm"
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
use_surrogate = opt_norm == 2 && threshold ~= 0;                                % (new, November 2026) option 4: threshold from surrogates
                                                                                % (threshold = 0: scan of "Optimal_secondaryPeak", SD threshold)
if use_surrogate
    smoothWindow = SMOOTH_FACTOR(1) + smoothWindowAdd;                          % a single smoothing (each smoothing would need its own surrogates)
    [IndBurst, n_peaks_smooth, n_peaksLevel1_smooth, smoothed_spikesBin, surrogate_info, p_event] = surrogate_events(spikeData, range, bin, ...
        smoothWindow, time_binned, spikes, opt_secondary, N_SURROGATE, SURROGATE_SEED, ALPHA_LARGE, PERCENTILE_SMALL, ALPHA_SECONDARY);
    n_smoothing = 0;                                                            % (the search with the SD threshold below is skipped)
end
if ~opt_secondary && ~use_surrogate
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
        % (revised November 2026; user's choices) the tail of each large event lasts until the signal returns to its baseline (the median
        % of the signal), and at most TAIL_FACTOR times the duration of the event after its end, whichever comes first. The tail is removed
        % to measure the noise (before, the event plus the same duration again, which on some recordings left most of the tail in) and, below,
        % it is closed to the small events (before, they were searched in the whole signal, so the bumps on the tails, where the signal is
        % higher and noisier than the baseline that sets the threshold, became small events)
        tail_end = tail_limits(smoothedSpikes, start_Nburst, end_Nburst, TAIL_FACTOR);
        for j=1:contador
            smoothedSpikes_cropped(start_Nburst(j):tail_end(j)) = nan;
        end
%        threshold0 = max(threshold0, MIN_THRESHOLD);
        threshold0 = max(threshold, MIN_THRESHOLD);
        threshold_peak_low = std(smoothedSpikes_cropped,'omitnan')*threshold0;
        threshold_second = min(threshold_peak_low, threshold_peak);

        %fprintf("Running 2nd round of Population Peak Detection with an absolute threshold reduced by a factor of %f \n", threshold_peak/threshold_second)

        [pks2, locs2, widths2, proms2] = findpeaks(smoothedSpikes, MinPeakProminence=threshold_second);
        % (new, November 2026; user's choice) the small candidates (below the threshold of the large events) whose peak is inside a large
        % event (from its start) or in its tail are part of that event, not small events (from the end only, the bumps between the peak and
        % a late end stayed: the merge of overlapping edges does not join two adjacent peaks, their edges stop at the valley between them)
        in_tail = false(numel(locs2),1);
        for j=1:contador
            in_tail = in_tail | (locs2(:) >= start_Nburst(j) & locs2(:) <= tail_end(j));
        end
        drop = in_tail & proms2(:) <= threshold_peak;
        pks2(drop) = [];    locs2(drop) = [];   widths2(drop) = [];     proms2(drop) = [];

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
%       (revised November 2026, user's choice) the events whose edges (20% level) overlap are merged into one event, which spans all of them;
%       before, one of them was discarded (the lower level one, or the first one when both had the same level). The peak of the merged
%       event is the one of its highest level, and then highest, candidate; each candidate is refined below, and the edges of the merged
%       event are the outermost ones

    contador_final = numel(locs_final);
    start_burst = time_binned(start_Nburst_final)+bin/2;
    end_burst = time_binned(end_Nburst_final)-bin/2;
    group = zeros(contador_final,1);                                        % merged event of each candidate
    n_group = 0;    group_end = -inf;
    for k=1:contador_final
        if start_Nburst_final(k) > group_end                                % no overlap with the event being built: a new event
            n_group = n_group + 1;
            group_end = end_Nburst_final(k);
        else                                                                % overlap: merged
            group_end = max(group_end, end_Nburst_final(k));
        end
        group(k) = n_group;
    end
    peak_final = nan(n_group,1);    level_final = nan(n_group,1);
    start_final = nan(n_group,1);   end_final = nan(n_group,1);
    members = cell(n_group,1);
    for g=1:n_group
        m = find(group == g);
        top = m(peak_level(m) == max(peak_level(m)));                       % candidates of the highest level
        [~, b] = max(smoothedSpikes(locs_final(top)));                      % and then the highest one
        best = top(b);
        peak_final(g) = locs_final(best);
        level_final(g) = peak_level(best);
        start_final(g) = min(start_burst(m));
        end_final(g) = max(end_burst(m));
        members{g} = struct('start', start_burst(m), 'end', end_burst(m), 'best', find(m == best));
    end

    IndBurst(sm).peak = peak_final;
    IndBurst(sm).start = start_final;
    IndBurst(sm).end = end_final;
    IndBurst(sm).level = level_final;
    IndBurst(sm).members = members;
    n_peaks_smooth(sm) = n_group;
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
    start20=nan(n_peaks_final,1); last20=nan(n_peaks_final,1);  start50=nan(n_peaks_final,1); last50=nan(n_peaks_final,1);  % (new, November 2026)
    spike_merged = mergeSpikeTrains(units, spikeTimestamp);                 % (revised October 2026) spike trains of the selected clusters, merged once for all the peaks 
    for k=1:n_peaks_final       
        if use_surrogate                                                    % (new, November 2026) option 4: the edges at 50% of the peak (event window), and at 20%
            peakCharact = struct('peak', IndBurst(sm_optimal).peakTime(k), 'peakAmp', IndBurst(sm_optimal).peakAmp(k), ...
                'start', IndBurst(sm_optimal).start50(k), 'end', IndBurst(sm_optimal).end50(k), ...
                'start20', IndBurst(sm_optimal).start20(k), 'end20', IndBurst(sm_optimal).end20(k), ...
                'start50', IndBurst(sm_optimal).start50(k), 'end50', IndBurst(sm_optimal).end50(k));
        else                                                                % methods 1-3: each candidate of the (merged) event is refined
            cand = IndBurst(sm_optimal).members{k};
            n_cand = numel(cand.start);
            s_c = nan(n_cand,1);    e_c = nan(n_cand,1);    s50_c = nan(n_cand,1);  e50_c = nan(n_cand,1);
            for c=1:n_cand
                refined = refinePeakDetection(units, spikeTimestamp, cand.start(c), cand.end(c), bin, START_LEVEL, END_LEVEL, opt_norm, opt_plot, spike_merged);
                s_c(c) = refined.start;     e_c(c) = refined.end;   s50_c(c) = refined.start50;     e50_c(c) = refined.end50;
                if c == cand.best
                    peakCharact = refined;                                  % peak of the event: its best candidate
                end
            end
            peakCharact.start = min(s_c);   peakCharact.end = max(e_c);     % edges of the event: the outermost ones (20% level)
            peakCharact.start20 = peakCharact.start;    peakCharact.end20 = peakCharact.end;
            peakCharact.start50 = min(s50_c);   peakCharact.end50 = max(e50_c);   % half-width (50% level), reported as a metric
        end
        peak(k) = peakCharact.peak;                                 % temporal location of the peak of event k, in sec 
        peakAmp(k) = peakCharact.peakAmp;                           % Amplitude of the peak of event k 
        peakIdx(k) = IndBurst(sm_optimal).peak(k);                  % location of the peak of the event k
        start(k) = peakCharact.start;                               % beginning of the event k, in sec
        last(k) = peakCharact.end;                                  % end of the event k, in sec 
        level(k) = IndBurst(sm_optimal).level(k);                   % level of the event k (=0, small event; =1, large event)  
        start20(k) = peakCharact.start20;   last20(k) = peakCharact.end20;  % (new, November 2026) edges at 20% and 50% of the height of the peak
        start50(k) = peakCharact.start50;   last50(k) = peakCharact.end50;
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
            PopulationActivityDetection.IndivBurst(idx).start20 = start20(k);       % (new, November 2026) edges at 20% of the height of the peak, in sec
            PopulationActivityDetection.IndivBurst(idx).last20 = last20(k);
            PopulationActivityDetection.IndivBurst(idx).start50 = start50(k);       % (new, November 2026) edges at 50% of the height of the peak (half-width), in sec
            PopulationActivityDetection.IndivBurst(idx).last50 = last50(k);
            if use_surrogate
                PopulationActivityDetection.IndivBurst(idx).pSurrogate = p_event(k);    % (new, November 2026) p-value of the peak (max statistic)
                PopulationActivityDetection.IndivBurst(idx).coreLast = min(last(k), peak(k) + (peak(k) - start(k)));   % (new, November 2026) end of
                                                                                        % the core of the event (start -> peak, plus the same time after
                                                                                        % the peak): the participation of the units is tested on it
            end
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
    PopulationActivityDetection.IndivBurst.start20 = [];
    PopulationActivityDetection.IndivBurst.last20 = [];
    PopulationActivityDetection.IndivBurst.start50 = [];
    PopulationActivityDetection.IndivBurst.last50 = [];
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
    PopulationActivityDetection.IndivBurst.start20 = [];
    PopulationActivityDetection.IndivBurst.last20 = [];
    PopulationActivityDetection.IndivBurst.start50 = [];
    PopulationActivityDetection.IndivBurst.last50 = [];
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
if opt_norm == 2
    PopulationActivityDetection.surrogate = surrogate_info;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IndBurst, n_events, n_large, smoothed, info, p_event] = surrogate_events(spikeData, range, bin, smoothWindow, time_binned, spikes, ...
    opt_secondary, N_SURROGATE, SURROGATE_SEED, ALPHA_LARGE, PERCENTILE_SMALL, ALPHA_SECONDARY)

% function to detect the events of the participation signal with the thresholds from the surrogates (option 4, opt_norm = 2)
% (new, November 2026). Only absolute values of the participation are used (no prominence, user's choice): each stretch of the smoothed
% signal above the lower threshold is a candidate, its highest point is its peak (level and p-value), and its edges are where the signal
% falls to 50% of the height of the peak above the median of the signal (event window), and to 20% (reported). Candidates whose 50% edges
% overlap are merged into one event, which spans all of them, with the higher peak (user's choices). Tail fragments (stretches with no rise
% of their own, right after the signal left the threshold) are dropped

smoothed =smoothdata(spikes, 'gaussian', smoothWindow);
null = participationSurrogate(spikeData, range, bin, smoothWindow, N_SURROGATE, SURROGATE_SEED);
maxima = sort(null.maxima, 'descend');
if opt_secondary
    alpha = ALPHA_SECONDARY;
else
    alpha = ALPHA_LARGE;
end
k_alpha = floor(alpha*(N_SURROGATE + 1));                                   % peak > maxima(k_alpha)  <=>  p-value <= alpha
if k_alpha >= 1
    threshold_large = maxima(k_alpha);
else
    threshold_large = inf;                                                  % (too few surrogates for this alpha)
end
if opt_secondary                                                            % secondary events: max statistic only
    threshold_small = threshold_large;
else                                                                        % small events: percentile of all the surrogate bins pooled
    cdf = cumsum(null.pooledCount)/sum(null.pooledCount);
    threshold_small = min(find(cdf >= PERCENTILE_SMALL/100, 1)*null.histWidth, threshold_large);
end
info = struct('nSurrogate', N_SURROGATE, 'seed', SURROGATE_SEED, 'alpha', alpha, 'percentileSmall', PERCENTILE_SMALL, ...
    'thresholdLarge', threshold_large, 'thresholdSmall', threshold_small);
if opt_secondary
    info.percentileSmall = nan;
end

% stretches above the lower threshold (consecutive bins above it are one candidate), and the peak of each one
y = smoothed(:);
t = time_binned(:);
n_bins = numel(y);
above = y > threshold_small;
d = diff([false; above; false]);
run_start = find(d == 1);
run_end = find(d == -1) - 1;
n_events = numel(run_start);
pks = zeros(n_events,1);    locs = zeros(n_events,1);
for r=1:n_events
    [pks(r), jj] = max(y(run_start(r):run_end(r)));                        % peak: highest participation of the stretch
    locs(r) = run_start(r) + jj - 1;
end
% tail fragments are not events (new, November 2026; user's choice): a stretch with no rise of its own (its highest point is its first bin)
% that starts right where the signal had just left the threshold (a dip of at most MAX_DIP_TAIL bins after the previous stretch, and not
% higher than where the signal left it), or at the first bin of the window (the tail of an event before the window, e.g., the large event
% before a gap in the search for secondary events), is the decaying tail of the event before it. A stretch with its own rise is kept
MAX_DIP_TAIL = 2;                                                           % maximum length of the dip below the threshold, in bins
keep = true(n_events,1);
for r=1:n_events
    if locs(r) == run_start(r)                                              % no rise: the peak is the first bin of the stretch
        if r == 1
            keep(r) = run_start(r) > 1;                                     % tail of an event before the window
        else
            keep(r) = ~((run_start(r) - run_end(r-1) - 1) <= MAX_DIP_TAIL && y(run_start(r)) <= y(run_end(r-1)));
        end
    end
end
run_start = run_start(keep);    run_end = run_end(keep);    pks = pks(keep);    locs = locs(keep);
n_events = numel(pks);
% edges at 50% (event window) and 20% (reported) of the height of the peak above the median of the signal (half-height above the typical
% level of this recording, or of this gap for the secondary events; user's choice), searched up to the neighbouring peaks; when the level
% is not reached, the lowest point in between (as the other methods). Half of the bins of the signal are below its median, so the levels
% are crossed close to the peak even when the participation stays high for a long time (e.g., the persistent firing after large events)
baseline = median(y);
start20 = zeros(n_events,1);    end20 = start20;    start50 = start20;  end50 = start20;
for r=1:n_events
    lb = 1;     rb = n_bins;
    if r > 1, lb = locs(r-1); end
    if r < n_events, rb = locs(r+1); end
    [start20(r), end20(r)] = level_edges(y, t, locs(r), baseline + 0.2*(pks(r) - baseline), lb, rb);
    [start50(r), end50(r)] = level_edges(y, t, locs(r), baseline + 0.5*(pks(r) - baseline), lb, rb);
end
info.baseline = baseline;                                                   % median of the signal (reference of the edges)
% events whose 50% edges overlap are merged into one event, which spans all of them (user's choice), with the higher peak
keep = true(n_events,1);    last = 1;
for r=2:n_events
    if start50(r) <= end50(last)
        start50(last) = min(start50(last), start50(r));     end50(last) = max(end50(last), end50(r));
        start20(last) = min(start20(last), start20(r));     end20(last) = max(end20(last), end20(r));
        if pks(r) > pks(last)
            pks(last) = pks(r);     locs(last) = locs(r);
        end
        keep(r) = false;
    else
        last = r;
    end
end
pks = pks(keep);    locs = locs(keep);
start20 = start20(keep);    end20 = end20(keep);    start50 = start50(keep);    end50 = end50(keep);
n_events = numel(pks);
level = double(pks > threshold_large);
p_event = (1 + sum(null.maxima(:)' >= pks(:), 2))/(N_SURROGATE + 1);      % p-value of each peak (max statistic)

IndBurst.peak = locs;                                                       % bin of the peak
IndBurst.start = start50;                                                   % event window: 50% of the height of the peak
IndBurst.end = end50;
IndBurst.level = level;
IndBurst.peakTime = t(locs);                                                % time of the peak, in sec
IndBurst.peakAmp = pks;                                                     % participation at the peak
IndBurst.start20 = start20;                                                 % edges at 20% of the height of the peak (reported)
IndBurst.end20 = end20;
IndBurst.start50 = start50;
IndBurst.end50 = end50;
n_large = sum(level);
smoothed = {smoothed};

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s_t, e_t] = level_edges(y, t, loc, level, lb, rb)

% function to find where the signal falls below "level" to the left and to the right of the peak (bin "loc"), by linear interpolation
% between bins, searching between the bins "lb" and "rb" (the neighbouring peaks); when the signal does not fall below the level, the 
% lowest point on that side (option 4, new November 2026)

k = find(y(lb:loc) < level, 1, 'last');
if ~isempty(k)
    k = k + lb - 1;
    s_t = t(k) + (level - y(k))/(y(k+1) - y(k))*(t(k+1) - t(k));
else
    [~, m] = min(y(lb:loc));
    s_t = t(lb + m - 1);
end
k = find(y(loc:rb) < level, 1, 'first');
if ~isempty(k)
    k = k + loc - 1;
    e_t = t(k-1) + (y(k-1) - level)/(y(k-1) - y(k))*(t(k) - t(k-1));
else
    [~, m] = min(y(loc:rb));
    e_t = t(loc + m - 1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tail_end = tail_limits(y, start_idx, end_idx, factor)

% (new, November 2026) last bin of the tail of each event: the first bin after its end where the signal is at or below its median (return to
% baseline), and at most factor*(end - start) bins after its end, whichever comes first (user's choice)

n = numel(y);
base = median(y, 'omitnan');
tail_end = zeros(numel(end_idx),1);
for j=1:numel(end_idx)
    cap = min(n, end_idx(j) + round(factor*(end_idx(j) - start_idx(j))));
    back = find(y(end_idx(j):cap) <= base, 1);
    if isempty(back)
        tail_end(j) = cap;
    else
        tail_end(j) = end_idx(j) + back - 1;
    end
end
end
