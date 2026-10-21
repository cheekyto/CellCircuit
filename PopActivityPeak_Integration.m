function PopActivityPeak_Integration = PopActivityPeak_Integration(units, spikeTimestamp, BurstDetected, range, bin, threshold, opt_norm, opt_distance, opt_binopti, option_level, opt_plot)

% function to detect population activation in a entire spike train data, using the average population firing rate and calling the function "PopulationActivityDetection" 
% Its goal is to run "PopulationActivityDetection" with different settings to perform the analysis of an entire recording

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   structure array spikeTimestamp (its length is the number of recording channels), with the following fields
% 	.recLength:	length of the recording, in seconds.
% 	.nUnits:	total number of identified units in the channel
% 	.timeStamp:	time stamp for all units (nUnits) in the channel, in msec

%   units(j):       =0, if channel j should be ignored;     =1, if channel j should be incoporated for population activity
%   range:          range of the recording (range(1)=start; range(2)=end), in sec
%   bin:            size of the temporal bin to calculate the population firing rate, in sec        
%   threshold:      in terms of the standard deviation 
%   opt_norm:       =0, No normalization; =1: normalization of the firing rate
%   opt_distance    =0, no cleaning of the primary events based on the distance; =1, remove the primary events that are too close from the previous one (echo events)
%   opt_binopti     =0, do not run the bin optimization for the secondary peaks; =1, run the bin optimization
%                   =value, run the detection with that bin value
%   option_level    =1, just analyze the level1 peaks to determine the optimal bin size and the secondary peaks
%                   =0, analyze all peaks (level0 and level1) to determine the optimal bin size
%   opt_plot:       =0, no plot; =1: plot the data and the identified network bursts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PopActivityPeak_Integration.spikesBin(i,1):           total number of spikes in bin i
% PopActivityPeak_Integration.time(i):                  midpoint of the bin i, in sec 
% PopActivityPeak_Integration.nBurst:                   total number of network bursts 
% PopActivityPeak_Integration.IndivBurst(j).peak:       peak point of network burst j, in sec (after refinement)       
% PopActivityPeak_Integration.IndivBurst(j).start:      start point of network burst j, in sec (after refinement)
% PopActivityPeak_Integration.IndivBurst(j).last:       last point of network burst j, in sec (after refinement)
% PopActivityPeak_Integration.IndivBurst(j).peakIdx:    bin index of the peak point of network burst j

% PopulationActivityDetection.Secondary.nBurst 
% PopulationActivityDetection.Secondary.IndivBurst(i).peak       
% PopulationActivityDetection.Secondary.IndivBurst(i).start      
% PopulationActivityDetection.Secondary.IndivBurst(i).last       
% PopActivityPeak_Integration.Secondary.nBurst_per_gap(j):      total number of secondary peaks in the gap j 
% PopActivityPeak_Integration.Secondary.IndivBurst(j).parent:   parent primary peak to which the seconday peak j is associated with

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                            % factor to rescale from msec to sec
BIN_DEFAULT1 = 0.100;                             % default size of the temporal bin to calculate the population firing rate, in sec (default: 0.10 sec) 
BIN_DEFAULT2 = 0.150;                             % default size of the temporal bin to calculate the population firing rate, in sec (default: 0.15 sec) 
apply_parameters('PopActivityPeak_Integration');                            % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Search for the primary peaks

PopActivityPeak_Integration = PopulationActivityDetection(units, spikeTimestamp, range, bin, threshold, 0, opt_norm, 0, 0);
if opt_distance
    PopActivityPeak_Integration = evaluatePrimaryBursts(PopActivityPeak_Integration);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Search for secondary peaks between the primary peaks

switch option_level
    case 0      % all peaks (large and small)
        n_Burst = PopActivityPeak_Integration.nBurst;
        burst_start = nan(n_Burst,1);    burst_last = nan(n_Burst,1);    burst_peak = nan(n_Burst,1);
        for i=1:n_Burst
            burst_start(i) = PopActivityPeak_Integration.IndivBurst(i).start;    
            burst_last(i) = PopActivityPeak_Integration.IndivBurst(i).last;    
            burst_peak(i) = PopActivityPeak_Integration.IndivBurst(i).peak ;
        end

    case 1      % only consider level1 peaks (large peaks)
        n_Burst = PopActivityPeak_Integration.nBurst_level1;
        burst_start = nan(n_Burst,1);    burst_last = nan(n_Burst,1);    burst_peak = nan(n_Burst,1);
        idx = 0;
        for i=1:PopActivityPeak_Integration.nBurst
            if PopActivityPeak_Integration.IndivBurst(i).level == 1
                idx = idx + 1;
                burst_start(idx) = PopActivityPeak_Integration.IndivBurst(i).start;    
                burst_last(idx) = PopActivityPeak_Integration.IndivBurst(i).last;    
                burst_peak(idx) = PopActivityPeak_Integration.IndivBurst(i).peak ;
            end
        end
end

if n_Burst > 0
    
    spikeData = spikeTimestamp;
    for k=size(spikeData,2):-1:1
        if ~units(k)
            spikeData(k)= [];
        end
    end
    n_clusters = size(spikeData,2);                                            % total number of clusters
    n_window = n_Burst;
    window_edge = nan(n_window, 2);
    
    for k=1:n_window
        if bin < BIN_DEFAULT1 
%            window_edge(k,1) = PopActivityPeak_Integration.IndivBurst(k).last + (PopActivityPeak_Integration.IndivBurst(k).last-PopActivityPeak_Integration.IndivBurst(k).start)/2;        
            window_edge(k,1) = burst_last(k) + (burst_last(k)-burst_start(k))/2;        
        else
            if bin < BIN_DEFAULT2 
%                window_edge(k,1) = PopActivityPeak_Integration.IndivBurst(k).last;        
                window_edge(k,1) = burst_last(k);        
            else
%                window_edge(k,1) = PopActivityPeak_Integration.IndivBurst(k).peak + (PopActivityPeak_Integration.IndivBurst(k).last-PopActivityPeak_Integration.IndivBurst(k).peak)/2;
                window_edge(k,1) = burst_peak(k) + (burst_last(k)-burst_peak(k))/2;
            end
        end
        if k == n_window 
            window_edge(k,2) = range(2);                                        % end of the recording
        else
%            window_edge(k,2) = PopActivityPeak_Integration.IndivBurst(k+1).start;            
            window_edge(k,2) = burst_start(k+1);            
        end
    end
    
    PopActivityPeak_Integration.Secondary = Optimal_secondaryPeak(units, spikeTimestamp, range, BurstDetected, PopActivityPeak_Integration, window_edge, bin, threshold, opt_norm, opt_distance, opt_binopti);

  if 0
    figure()
    bar(time_binned, popActivity); hold on
    plot(time_binned, ones(numel(time_binned))*umbral,'r--', 'LineWidth',2); hold on
    plot([search_start; search_start], [0 max(popActivity)], 'c--','LineWidth',2);    hold on
    plot([ActDrop.time; ActDrop.time], [0 max(popActivity)], 'r--','LineWidth',2);    hold on

    figure()
    for k=1:n_clusters
        data = spikeData(k);
        plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+k-1,'|',"MarkerEdgeColor",[0 0 0]);     hold on;     clear data
    end 
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     options to plot the process

if opt_plot
    figure()
    spikes = PopActivityPeak_Integration.spikesBin;
    bar(PopActivityPeak_Integration.time, spikes); hold on;
    smoothedSpikes = smoothdata(spikes, 'gaussian', 3);                 
    plot(PopActivityPeak_Integration.time, smoothedSpikes,'m--' ); hold on;
    for i=1:PopActivityPeak_Integration.nBurst
        plot(PopActivityPeak_Integration.IndivBurst(i).peak, spikes(PopActivityPeak_Integration.IndivBurst(i).peakIdx)*1.1,'ro', 'MarkerFaceColor','r'); hold on
        plot(PopActivityPeak_Integration.IndivBurst(i).start, spikes(PopActivityPeak_Integration.IndivBurst(i).peakIdx)/2,'k*'); hold on
        plot(PopActivityPeak_Integration.IndivBurst(i).last, spikes(PopActivityPeak_Integration.IndivBurst(i).peakIdx)/2,'g*'); hold on
    end

    if PopActivityPeak_Integration.Secondary.nBurst>0
        for j=1:PopActivityPeak_Integration.Secondary.nBurst
            aux = find(PopActivityPeak_Integration.time > PopActivityPeak_Integration.Secondary.IndivBurst(j).peak);
            if numel(aux>0)
                plot(PopActivityPeak_Integration.Secondary.IndivBurst(j).peak, spikes(aux(1))*1.1,'co', 'MarkerFaceColor','c'); hold on
                plot(PopActivityPeak_Integration.Secondary.IndivBurst(j).start, spikes(aux(1))/2,'k*'); hold on
                plot(PopActivityPeak_Integration.Secondary.IndivBurst(j).last, spikes(aux(1))/2,'g*'); hold on
            end 
        end
    end
end

else
        PopActivityPeak_Integration.Secondary.nBurst = nan; 
        PopActivityPeak_Integration.Secondary.IndivBurst.peak = [];       
        PopActivityPeak_Integration.Secondary.IndivBurst.start = [];      
        PopActivityPeak_Integration.Secondary.IndivBurst.last = [];
        PopActivityPeak_Integration.Secondary.nBurst_per_gap = [];
        PopActivityPeak_Integration.Secondary.bin = []; 
end
end