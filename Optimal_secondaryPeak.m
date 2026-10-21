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
%   opt_norm:           =0, No normalization; =1: normalization of the firing rate
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
                if opt_norm
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
                    if opt_norm
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
                secondary_peaks(k,1) = PopulationActivityDetection(units, spikeTimestamp, [window_edge(k,1); endtime], binSize, threshold2, SMOOTH_ADD, opt_norm, 1, 0);
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
                    end
                end
            end
            Current_PopActivityPeak.nBurst = nbursts; 
            Current_PopActivityPeak.nBurst_per_gap = nburst_per_gap;
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