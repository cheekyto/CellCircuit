function Optimal_secondaryPeak = Optimal_secondaryPeak(units, spikeTimestamp, range, BurstDetected, PopActivityPeak, window_edge, bin1, threshold, opt_norm, opt_distance, opt_binopti)

% function to detect secondary peaks from the population activation in a entire spike train data, 
% Its goal is to run "PopulationActivityDetection" with different settings,
% find the optimal bin size and use it to report the analysis result

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   PopActivityPeak:    output of the function "PopulationActivityDetection", containing the information about primary peaks
%   window_edge:        matrix with the temporal location of each primary peak
%   bin1:               bin size used for the detection of primary peaks, in sec
%   units(j):           =0, if channel j should be ignored;     =1, if channel j should be incoporated for population activity
%   opt_norm:           =0, No normalization; =1: normalization of the firing rate; =2: participation, with surrogate thresholds
%   range:          range of the recording (range(1)=start; range(2)=end), in sec
%   threshold:      in terms of the standard deviation 
%   opt_norm:       =0, No normalization; =1: normalization of the firing rate
%   opt_distance    =0, no cleaning of the primary events based on the distance; =1, remove the primary events that are too close from the previous one (echo events)
%   opt_binopti:    =0, do not run the bin optimization; =1: run the bin optimization
%                   =value, run the detection with that bin value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optimal_secondaryPeak.nBurst 
% Optimal_secondaryPeak.IndivBurst(i).peak       
% Optimal_secondaryPeak.IndivBurst(i).start      
% Optimal_secondaryPeak.IndivBurst(i).last       
% Optimal_secondaryPeak.nBurst_per_gap(j):      total number of secondary peaks in the gap j 
% Optimal_secondaryPeak.IndivBurst(j).parent:   parent primary peak to which the seconday peak j is associated with
% Optimal_secondaryPeak.IndivBurst(j).dropOwn:  (new, November 2026) drop of the signal between the large (parent) event and the secondary
%                                               peak j, at the bin of the secondary search, as a fraction of the height of peak j above the
%                                               median of the gap (1 = the signal falls to the median)
% Optimal_secondaryPeak.IndivBurst(j).dropZ:    (option 4) the same drop, in units of its sampling noise (binomial noise of the participation)
% Optimal_secondaryPeak.IndivBurst(j).searchStart: start of the search in the gap of peak j, in sec
% Optimal_secondaryPeak.searchWindows(k,:):    (new, November 2026) [start end] of the search in the gap k, in sec
% Optimal_secondaryPeak.bin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                            % factor to rescale from msec to sec
%BIN_SECONDARY = 0.50;                           % reduction factor for the bin size (default: 0.30)
%FACTOR_BIN_SECONDARY = 2.5;                     % (default: 3)
SMOOTH_ADD = 0;                                 % increament of the standard smoothing window length, to look for secondary peaks, in number of bins (default: 1; even for symmetric smoothing)
BIN_MIN = 0.010;                                % minimum size of the bin (default: 0.008 sec)                                           

WINDOW_FRACTION = 1/4;                          % to partition the gap between peaks and assess the general level of activity  right after primary peaks (default: 1/4)
THRESHOLD_BROAD = 1.80;                         % threshold to assess the activity right after primary peaks; ratio to the last partition (default: 1.80)
THRESHOLD_WINDOWS = 0.40;                       % minimum number of primary peaks with a signifnicant tail activity (default: 0.50-0.40)

PERCENTILE_2GAP = 0.20;                         % percentile of the distance between potential secondary peaks to determine the "time scale" (default: 0.20)
SEARCH_VALLEY = 1.0;                            % proportion of the "typical_gap" to start the search of secondary peaks (default: 1.0)

% ???                                           % factor to calculate the prominence in the search of secondary peaks, in terms of the std (default: 1.0)
                                                % NOTE: if there is secondary peaks, the std is considerable (it is not the regular baseline std)
THRESHOLD_SECONDARY = 0.60;                     % factor to scale the threshold for the primary peaks and use it for the secondary peaks (default: 0.6-0.5)

MIN_CHAN = 6;                                   % minimum number of participating channels/units to classify an organoid as network bursting (default: 5)
                                                % (revised October 2026: 5 -> 6, the same value as the final quality control of the secondary events in "NetworkBurstDetection")
MIN_CHAN_SECONDARY = 0.25;                      % minimum number of participating channels/units to accept the secondary events, in terms of the number of participating units in the primary events (default: 0.25)
BIN_DEFAULT = 0.100;                             % default size of the temporal bin to calculate the population firing rate, in sec (default: 0.100 sec) 
BIN2_FACTOR_LOW = [0.40; 0.20; 1.00];           % reduction factors tested for the bin size, when the bin of the primary events is < BIN_DEFAULT (opt_binopti = 1)
BIN2_FACTOR_HIGH = [0.30; 0.12; 0.75; 0.05];    % reduction factors tested for the bin size, when the bin of the primary events is >= BIN_DEFAULT (opt_binopti = 1)
BIN2_FACTOR_FIXED = 0.30;                       % reduction factor for the bin size, without optimization (opt_binopti = 0) (default: 0.30)
                                                % (revised October 2026: the three factors above were written in the "switch" below; same values)
% (new, November 2026; user's choices) the secondary events are separated from their large event at the bin of the secondary search, where the
% drop between them is visible (at the bin of the large events, the large event and its reverberations blur into a single mountain)
SMOOTH_FINE = 3;                                % smoothing of the signal at the bin of the secondary search, in number of bins (+SMOOTH_ADD); the same as
                                                % the first smoothing of "PopulationActivityDetection" (SMOOTH_FACTOR(1)) (default: 3)
START_LEVEL_FINE = 0.50;                        % the search starts where the large event falls to this fraction of its height above the median of the
                                                % gap, at the bin of the secondary search (default: 0.50)
SEPARATION_LEVEL = 0.50;                        % a secondary event is separate from its large event if, between the two peaks, the signal falls below this
                                                % fraction of the height of the secondary peak above the median of the gap (the 50% rule used to merge the
                                                % events of option 4) (default: 0.50)
MIN_DROP_Z = 3;                                 % the drop must also be at least MIN_DROP_Z times its sampling noise, which is large with few units and fine
                                                % bins: binomial noise of the participation (option 4), Poisson noise of the spike counts (options 1-3; revised
                                                % November 2026, user's choice: before, option 4 only) (default: 3)
apply_parameters('Optimal_secondaryPeak');                                  % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

switch opt_binopti
    case 1
        if bin1 < BIN_DEFAULT
            bin2 = BIN2_FACTOR_LOW;
        else
            bin2 = BIN2_FACTOR_HIGH;
        end
    case 0 
        bin2 = BIN2_FACTOR_FIXED;                   % reduction factor for the bin size (default: 0.30)
    otherwise
        bin2 = opt_binopti/bin1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Search for secondary peaks between the primary peaks

    
spikeData = spikeTimestamp;
for k=size(spikeData,2):-1:1
    if ~units(k)
        spikeData(k)= [];
    end
end
n_bin = numel(bin2);                                                        % total number of bin sizes to be tested
n_window = PopActivityPeak.nBurst_level1;
n_clusters = size(spikeData,2);                                             % total number of clusters
PopActivityPeak0 = PopActivityPeak;

% (new, November 2026) peak of the large event that opens each gap, and all the spikes of the selected units sorted in time, to build the
% signal at the bin of the secondary search (start of the search and separation of the secondary events from their large event)
use_fine = false;
if n_window > 0 && isfield(PopActivityPeak0, 'IndivBurst') && isfield(PopActivityPeak0.IndivBurst, 'level')
    idx_large = find([PopActivityPeak0.IndivBurst.level] == 1);
    use_fine = numel(idx_large) == n_window;                                % the gaps follow the level-1 primary events (option_level = 1)
end
if use_fine
    peak_large = [PopActivityPeak0.IndivBurst(idx_large).peak];             % peak of the large event of each gap, in sec
    spk_t = cell(n_clusters,1);
    for i=1:n_clusters
        spk_t{i} = spikeData(i).timeStamp(:)*K_MSEC_SEC;                    % spike times, in sec
    end
    n_spk = cellfun(@numel, spk_t);
    all_t = vertcat(spk_t{:});
    all_u = repelem((1:n_clusters)', n_spk);
    all_u = all_u(:);                                                       % (fix, November 2026) a column also with a single cluster ("repelem" of a scalar gives a row)
    [all_t, order] = sort(all_t);   all_u = all_u(order);
    impulse = zeros(4*(SMOOTH_FINE + SMOOTH_ADD) + 1, 1);   impulse(2*(SMOOTH_FINE + SMOOTH_ADD) + 1) = 1;
    noise_factor = sqrt(sum(smoothdata(impulse, 'gaussian', SMOOTH_FINE + SMOOTH_ADD).^2));   % reduction of the noise by the smoothing
    clear spk_t order impulse
end

if n_window > 0
        
    nPeak2 = zeros(n_bin,1);
    for j=1:n_bin

        binSize = max(bin1*bin2(j),BIN_MIN);                                 % bin size j
        
        % Assess if there is significant activity after the main peak
        ratio_activity = nan(n_window,3);
        largest_gap = 0;
        for k=1:n_window                                                    % across all the gaps between primary peaks
            gap = window_edge(k,2) - window_edge(k,1);
            largest_gap = max([largest_gap, gap]);                          % largest gap (between primary peaks)
            aux1 = bin_spike_time(spikeData(1).timeStamp, binSize/K_MSEC_SEC, [window_edge(k,1), window_edge(k,1)+gap*WINDOW_FRACTION]);    
            aux2 = bin_spike_time(spikeData(1).timeStamp, binSize/K_MSEC_SEC, [window_edge(k,2)-gap*WINDOW_FRACTION*2, window_edge(k,2)]);    
            spikes_indiv1 = nan(numel(aux1.binMidPoint),n_clusters);         % matrix to save the binned data
            spikes_indiv2 = nan(numel(aux2.binMidPoint),n_clusters);         % matrix to save the binned data
            
            for i=1:n_clusters                                              % across all clusters
                aux1 = bin_spike_time(spikeData(i).timeStamp, binSize/K_MSEC_SEC, [window_edge(k,1), window_edge(k,1)+gap*WINDOW_FRACTION]);    
                aux2 = bin_spike_time(spikeData(i).timeStamp, binSize/K_MSEC_SEC, [window_edge(k,2)-gap*WINDOW_FRACTION*2, window_edge(k,2)]);    
%                aux2 = bin_spike_time(spikeData(i).timeStamp, binSize/K_MSEC_SEC, [window_edge(k,2)-gap*WINDOW_FRACTION*2, window_edge(k,2)-gap*WINDOW_FRACTION]);    
                if opt_norm == 2                                            % (new, November 2026) participation: the cluster fires in the bin
                    spikes_indiv1(:,i) = aux1.count > 0;
                    spikes_indiv2(:,i) = aux2.count > 0;
                elseif opt_norm
                    spikes_indiv1(:,i) = aux1.normCount;                    % normalized spike count
                    spikes_indiv2(:,i) = aux2.normCount;                    % normalized spike count
                else
                    spikes_indiv1(:,i) = aux1.count;                        % total spike count per bin
                    spikes_indiv2(:,i) = aux2.count;                        % total spike count per bin
                end
                clear aux1 aux2
            end
            spikes1 = mean(spikes_indiv1,2, 'omitnan');                     % mean Population spike count
            spikes2 = mean(spikes_indiv2,2, 'omitnan');                     % mean Population spike count
            ratio_activity(k,1) = mean(spikes1);                            % mean activity at the beginning of the gap                      
            ratio_activity(k,2) = mean(spikes2);                            % mean activity at the end of the gap
            ratio_activity(k,3) = mean(spikes1)/mean(spikes2);
        end

        aux_high = find(ratio_activity(:,3)>THRESHOLD_BROAD);
        proportion_high = numel(aux_high)/n_window;                         % proportion of gaps with significant initial activity

        % determine the analysis window after the main peaks
        if proportion_high > THRESHOLD_WINDOWS 
    
            aux = bin_spike_time(spikeData(1).timeStamp, binSize/K_MSEC_SEC, [window_edge(1,1), window_edge(1,1)+largest_gap]);    
            time_binned = aux.binMidPoint*K_MSEC_SEC - window_edge(1,1);    % time array with the midpoint for each bin, in sec
            activityMatrix = nan(n_window, numel(aux.binMidPoint));
            for k=1:n_window
                endtime = min(window_edge(k,2), range(2));
                aux = bin_spike_time(spikeData(1).timeStamp, binSize/K_MSEC_SEC, [window_edge(k,1), endtime]);    
                spikes_indiv = nan(numel(aux.binMidPoint),n_clusters);      % matrix to save the binned data
                for i=1:n_clusters                                          % across all clusters
                    aux = bin_spike_time(spikeData(i).timeStamp, binSize/K_MSEC_SEC, [window_edge(k,1), endtime]);    
                    if opt_norm == 2
                        spikes_indiv(:,i) = aux.count > 0;                  % (new, November 2026) participation: the cluster fires in the bin
                    elseif opt_norm
                        spikes_indiv(:,i) = aux.normCount;                  % normalized spike count
                    else
                        spikes_indiv(:,i) = aux.count;                      % total spike count per bin
                    end
                    clear aux
                end            
                spikes = mean(spikes_indiv,2, 'omitnan');                   % mean Population spike count
                activityMatrix(k,1:numel(spikes)) = spikes(:);
                clear spikes
            end

            popActivity = mean(activityMatrix,1,'omitnan');                 % averaging population firing across all gaps between primary peaks
            umbral = mean(ratio_activity(:,2), 'omitnan');
            ruido = std(ratio_activity(:,2), 'omitnan');
        
            % -- Use the distance between pseudo peaks to determine the "temporal scale" --
            idx = 0;
            second_gap = [];
            for k=1:n_window 
                if (window_edge(k,2)-window_edge(k,1)) > binSize*5
                    % NOTE: the threshold 0 is intended; "PopulationActivityDetection" sets it internally to its minimum value (MIN_THRESHOLD, in SD)
                    secondary_peak = PopulationActivityDetection(units, spikeTimestamp, window_edge(k,:), binSize, 0, SMOOTH_ADD, opt_norm, 1, 0);
                    if secondary_peak.nBurst >= 2                                   % (revised October 2026: ">= 2"; the windows with exactly 2 pseudo-peaks, 
                                                                                    %  which give one distance, were not used with "> 2")
                        for i = 2:secondary_peak.nBurst
                            idx = idx+1;
                            second_gap(idx) = secondary_peak.IndivBurst(i).peak - secondary_peak.IndivBurst(i-1).peak;
                        end
                    end
                end
                clear secondary_peak
            end
            if idx>0
                aux_gap = sort(second_gap, 'ascend');
                typical_gap = aux_gap(ceil(PERCENTILE_2GAP*numel(second_gap)));
                search_start = SEARCH_VALLEY * typical_gap;
            else
                search_start = largest_gap  / 10;
            end
            ActDrop = searchActDrop(popActivity, time_binned, activityMatrix, umbral, ruido, search_start);
                    
            % using the window found above, run the search for secondary peaks
            nbursts = 0;        idx = 0;
            PopActivityPeak = PopActivityPeak0;
            nburst_per_gap = nan(n_window,1);
            search_windows = nan(n_window,2);                                % (new, November 2026) window of the search in each gap [start end], in sec
            for k=1:n_window 
                gap = window_edge(k,2) - window_edge(k,1);                  % new line
                if k<n_window
                    endtime = min([window_edge(k,1)+ActDrop.time, window_edge(k,2)-gap*WINDOW_FRACTION]);
                else
                    endtime = min([window_edge(k,1)+ActDrop.time, range(2)]);
                end
                if ratio_activity(k,3)<1
                    threshold2 = threshold;
                else
                    threshold2 = threshold*THRESHOLD_SECONDARY;
                end
                % (revised November 2026; user's choice) the search starts where the large event falls to START_LEVEL_FINE of its height, at the
                % bin of the secondary search; before, it started at a point set by the bin of the large events (after its end, at its end, or
                % halfway down its decay), which missed the first reverberations or included the decay. If that point is not found, as before
                search_start = window_edge(k,1);
                fine = [];
                if use_fine
                    fine = fine_signal(all_t, all_u, n_spk, n_clusters, peak_large(k) - bin1, min(window_edge(k,2), range(2)), binSize, opt_norm, ...
                        SMOOTH_FINE + SMOOTH_ADD, peak_large(k), bin1, START_LEVEL_FINE, noise_factor);
                    if ~isempty(fine) && ~isnan(fine.endTime) && fine.endTime < endtime - 5*binSize
                        search_start = fine.endTime;
                    end
                end
                secondary_peaks(k,1) = PopulationActivityDetection(units, spikeTimestamp, [search_start; endtime], binSize, threshold2, SMOOTH_ADD, opt_norm, 1, 0);
                search_windows(k,:) = [search_start, endtime];
                % (new, November 2026; user's choice) each secondary event must be separate from its large event at this bin: between the two
                % peaks, the signal falls below SEPARATION_LEVEL of the height of the secondary peak (and, for option 4, by more than MIN_DROP_Z
                % times its sampling noise); otherwise it is part of the large event. This replaces the comparison with the end of the large
                % event (measured at the bin of the large events) in "NetworkBurstDetection"
                drop_own = nan(secondary_peaks(k,1).nBurst,1);   drop_z = drop_own;
                if ~isempty(fine) && secondary_peaks(k,1).nBurst > 0
                    separate = false(secondary_peaks(k,1).nBurst,1);
                    for i=1:secondary_peaks(k,1).nBurst
                        [separate(i), drop_own(i), drop_z(i)] = separation(fine, secondary_peaks(k,1).IndivBurst(i).peak, SEPARATION_LEVEL, ...
                            MIN_DROP_Z, n_clusters, opt_norm, noise_factor);
                    end
                    secondary_peaks(k,1).IndivBurst = secondary_peaks(k,1).IndivBurst(separate);
                    secondary_peaks(k,1).nBurst = sum(separate);
                    drop_own = drop_own(separate);   drop_z = drop_z(separate);
                end
                nbursts = nbursts + secondary_peaks(k,1).nBurst;
              
               nburst_per_gap(k) = secondary_peaks(k,1).nBurst;
               if secondary_peaks(k,1).nBurst>0
                    for i=1:secondary_peaks(k,1).nBurst
                        idx = idx+1;
                        PopActivityPeak.Secondary.IndivBurst(idx).start = secondary_peaks(k,1).IndivBurst(i).start;
                        PopActivityPeak.Secondary.IndivBurst(idx).last = secondary_peaks(k,1).IndivBurst(i).last;

                        Current_PopActivityPeak.IndivBurst(idx).peak = secondary_peaks(k,1).IndivBurst(i).peak;
                        Current_PopActivityPeak.IndivBurst(idx).peakAmp = secondary_peaks(k,1).IndivBurst(i).peakAmp;
                        Current_PopActivityPeak.IndivBurst(idx).start = secondary_peaks(k,1).IndivBurst(i).start;
                        Current_PopActivityPeak.IndivBurst(idx).last = secondary_peaks(k,1).IndivBurst(i).last;
                        Current_PopActivityPeak.IndivBurst(idx).tail = secondary_peaks(k,1).IndivBurst(i).tail;
                        Current_PopActivityPeak.IndivBurst(idx).parent = k;
                        Current_PopActivityPeak.IndivBurst(idx).level = 0;
                        Current_PopActivityPeak.IndivBurst(idx).start20 = secondary_peaks(k,1).IndivBurst(i).start20;   % (new, November 2026) edges at 20%
                        Current_PopActivityPeak.IndivBurst(idx).last20 = secondary_peaks(k,1).IndivBurst(i).last20;     % and 50% of the height of the peak
                        Current_PopActivityPeak.IndivBurst(idx).start50 = secondary_peaks(k,1).IndivBurst(i).start50;
                        Current_PopActivityPeak.IndivBurst(idx).last50 = secondary_peaks(k,1).IndivBurst(i).last50;
                        if isfield(secondary_peaks(k,1).IndivBurst, 'pSurrogate')                                       % (option 4) p-value of the peak,
                            Current_PopActivityPeak.IndivBurst(idx).pSurrogate = secondary_peaks(k,1).IndivBurst(i).pSurrogate; % and end of the core of the event
                            Current_PopActivityPeak.IndivBurst(idx).coreLast = secondary_peaks(k,1).IndivBurst(i).coreLast;
                        end
                        Current_PopActivityPeak.IndivBurst(idx).dropOwn = drop_own(i);      % (new, November 2026) separation from the large event
                        Current_PopActivityPeak.IndivBurst(idx).dropZ = drop_z(i);
                        Current_PopActivityPeak.IndivBurst(idx).searchStart = search_start;
                    end
                end
            end
            Current_PopActivityPeak.nBurst = nbursts; 
            Current_PopActivityPeak.nBurst_per_gap = nburst_per_gap;
            Current_PopActivityPeak.searchWindows = search_windows;          % (new, November 2026) used by the recruitment check (local null of the secondary events)
            Current_PopActivityPeak.bin = binSize; 
               
            % Quality Control on secondary events
            if nbursts > 0  && nbursts >= max(nPeak2)                  
                PopActivityPeak.Secondary.nBurst = nbursts;
                NetworkBurstPlayer = NetworkBurstPlayers(spikeTimestamp, range, PopActivityPeak, BurstDetected, 3, bin1, threshold, opt_distance, 0); % statistical assessment                 
                n_participating = numel(NetworkBurstPlayer.Channels);      % total number of participating units for primary events
                unit_min = max(MIN_CHAN, n_participating*MIN_CHAN_SECONDARY);   % minimum number of participating units for secondary events
                % it will check the secondary events, one by one     
                for k=nbursts:-1:1                                              % across all secondary events  
                    contador = 0;
                    for i=1:n_participating                                     % across all participating units for primary events
                        ch = NetworkBurstPlayer.Channels(i,1);
                        if numel(NetworkBurstPlayer.Secondary.timeStamp{k,ch})>0  
                            contador = contador + 1;                            % total number of clusters with activity during event k
                        end
                    end
                    if contador < unit_min                                     % if the number of participating units is enough                 
                        nbursts  = nbursts  - 1;
                        parent = Current_PopActivityPeak.IndivBurst(k).parent;
                        Current_PopActivityPeak.nBurst_per_gap(parent) = Current_PopActivityPeak.nBurst_per_gap(parent) - 1;
                        Current_PopActivityPeak.IndivBurst(k) = [];             % delete the secondary event k
                    end
                    clear contador
                end
                Current_PopActivityPeak.nBurst = nbursts;
            end
            nPeak2(j) = nbursts;
            all_PopActivityPeak{j} = Current_PopActivityPeak;       
            clear secondary_peaks n_bursts Current_PopActivityPeak          
        end
    end
   [max_peak2, best_bin] = max(nPeak2);                                    % identify the binsize that produces the larger number of secondary peaks
    
   if max_peak2>0
       Optimal_secondaryPeak = all_PopActivityPeak{best_bin};
   else
       Optimal_secondaryPeak.nBurst = nan; 
       Optimal_secondaryPeak.IndivBurst.peak = [];       
       Optimal_secondaryPeak.IndivBurst.start = [];      
       Optimal_secondaryPeak.IndivBurst.last = [];
       Optimal_secondaryPeak.IndivBurst.tail = [];
       Optimal_secondaryPeak.IndivBurst.level = [];
       Optimal_secondaryPeak.nBurst_per_gap = [];
   end
else
    Optimal_secondaryPeak.nBurst = nan; 
    Optimal_secondaryPeak.IndivBurst.peak = [];       
    Optimal_secondaryPeak.IndivBurst.start = [];      
    Optimal_secondaryPeak.IndivBurst.last = [];
    Optimal_secondaryPeak.IndivBurst.tail = [];
    Optimal_secondaryPeak.IndivBurst.level = [];
    Optimal_secondaryPeak.nBurst_per_gap = [];
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fine = fine_signal(all_t, all_u, n_spk, n_units, t0, t1, binSize, opt_norm, smoothWindow, t_large, half_window, start_level, noise_factor)

% (new, November 2026) population signal from the large event to the end of its gap, at the bin of the secondary search, with the same type of
% signal (opt_norm) and smoothing as the search, and where the large event falls to "start_level" of its height above the median of the gap
%   all_t, all_u:   spike times (in sec, sorted) and their unit;    n_spk(i): number of spikes of unit i (normalization, as "bin_spike_time")
%   t_large:        peak of the large event (at the bin of the large events), searched within +-half_window at this bin

fine = [];
edges = t0:binSize:t1;
if numel(edges) < 4
    return
end
n_bins = numel(edges) - 1;
first = find(all_t >= t0, 1);   last = find(all_t < edges(end), 1, 'last');
if isempty(first) || isempty(last) || last < first
    return
end
t = all_t(first:last);   u = all_u(first:last);
b = min(floor((t - t0)/binSize) + 1, n_bins);                               % bin of each spike
switch opt_norm
    case 2                                                                  % participation: fraction of units with at least one spike
        pairs = unique((u - 1)*n_bins + b);
        y = accumarray(mod(pairs - 1, n_bins) + 1, 1, [n_bins, 1])/max(n_units, 1);
        v = zeros(n_bins,1);                                                % (binomial noise: computed from y in "separation")
    case 1                                                                  % normalized rate: spikes of each unit / its total spikes, mean across units
        y = accumarray(b, 1./n_spk(u), [n_bins, 1])/max(sum(n_spk > 0), 1);
        v = accumarray(b, 1./n_spk(u).^2, [n_bins, 1])/max(sum(n_spk > 0), 1)^2;   % Poisson variance of y (each spike count used as its variance)
    otherwise                                                               % mean spike count
        y = accumarray(b, 1, [n_bins, 1])/max(n_units, 1);
        v = accumarray(b, 1, [n_bins, 1])/max(n_units, 1)^2;               % Poisson variance of y
end
y = smoothdata(y, 'gaussian', smoothWindow);
v = smoothdata(v, 'gaussian', smoothWindow)*noise_factor^2;                % variance after the smoothing (approximately)
tc = edges(1:end-1)' + binSize/2;
near = find(abs(tc - t_large) <= half_window);
if isempty(near)
    return
end
[peak_amp, j] = max(y(near));   i_large = near(j);                         % peak of the large event at this bin
base = median(y(tc >= t_large));                                           % median of the gap
end_time = nan;
below = find(y(i_large:end) < base + start_level*(peak_amp - base), 1);
if ~isempty(below) && peak_amp > base
    end_time = tc(i_large + below - 1);
end
fine = struct('t', tc, 'y', y, 'v', v, 'iLarge', i_large, 'base', base, 'endTime', end_time);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [separate, drop_own, drop_z] = separation(fine, t_peak2, level, min_z, n_units, opt_norm, noise_factor)

% (new, November 2026) is the secondary peak at t_peak2 separate from the large event of its gap? Between the two peaks the signal must fall
% below "level" of the height of the secondary peak above the median of the gap (drop_own > 1 - level), and the drop must be at least min_z
% times its sampling noise: binomial for the participation (option 4; a fraction of n_units), Poisson for the spike counts (options 1-3;
% revised November 2026), reduced by the smoothing

separate = false;   drop_own = nan;   drop_z = nan;
i2 = find(fine.t <= t_peak2, 1, 'last');
if isempty(i2) || i2 <= fine.iLarge
    return
end
h2 = max(fine.y(max(i2 - 1, fine.iLarge + 1):min(i2 + 1, numel(fine.y))));   % height of the secondary peak at this bin
[valley, iv] = min(fine.y(fine.iLarge:i2));                                    % lowest point between the two peaks
iv = fine.iLarge + iv - 1;
nb = max(i2 - 1, fine.iLarge + 1):min(i2 + 1, numel(fine.y));   [~, ip] = max(fine.y(nb));   ip = nb(ip);   % bin of the secondary peak
height = h2 - fine.base;
if height <= 0
    drop_own = 0;
    return
end
drop_own = (h2 - valley)/height;
separate = valley < fine.base + level*height;
if opt_norm == 2
    se = @(p) noise_factor*sqrt(max(p*(1 - p), 0)/max(n_units, 1));
    drop_z = (h2 - valley)/max(sqrt(se(min(h2, 1))^2 + se(min(valley, 1))^2), eps);
else                                                                        % (revised November 2026, user's choice) options 1-3: Poisson noise
    drop_z = (h2 - valley)/max(sqrt(fine.v(ip) + fine.v(iv)), eps);
end
separate = separate && drop_z >= min_z;
end
