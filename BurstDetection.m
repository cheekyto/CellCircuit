function BurstDetection = BurstDetection(signalSpike, option_method)

% function to detect bursts in a spike train data, using the method specified
% (revised October 2026: same results, but faster -> no "clear" inside the loops, no loop across the spikes of each burst)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   signalSpike(i,1):   time stamp of the spike i, in ms
% 	option_method:	    option input to specified the algorithm to be used for burst detection
%               =1      Maximal Interval

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BurstDetection.spikes(i,1)                    % time stamp of spike i, in ms
%               .spikes(i,2)                    % isi corresponding to the spike i, in ms     [isi(j)=T(j+1)-T(j)]     
%               .spikes(i,3)                    % # of burst to which the  spike i belongs (=0 if spike i is NOT part of a burst)
% BurstDetection.IndivBurst(j).start            % time stamp of the first spike of burst j, in ms
% BurstDetection.IndivBurst(j).last             % time stamp of the last spike of burst j, in ms
% BurstDetection.IndivBurst(j).peak             % time stamp of the peak of burst j, in ms
% BurstDetection.IndivBurst(j).median           % time stamp of the median/middle spike of the burst j, in ms
% BurstDetection.IndivBurst(j).isi(i)           % isi corresponding to the spike i of the burst j
% BurstDetection.parameter                      % structure with information about the burst detection method
% BurstDetection.nBurst                         % total number of bursts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RATIO_PEAK = 0.5;                               % proportion of spikes in the burst that are grouped to look for the peak of the burst (default: 0.5)
apply_parameters('BurstDetection');                                         % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

switch option_method

    case 1          % Maximal Interval
        MIN_SPIKE_BURST =  3;                   % minimum number of spikes within a burst (default: 3)
        ISI_START_BURST = 200;                  % threshold (maximum) ISI to identify the beginning of a burst, in ms (default: 200 -> 5Hz)
        ISI_PEAK_BURST = 150;                   % minimum peak ISI to accept a burst, in ms (default: 150 -> 6.7Hz [100ms -> 10Hz])
        ISI_END_BURST_INI = 300;                % initial threshold (maximum) ISI to identify the end of a burst, in ms (default: 300ms -> 3.3Hz)
        FACTOR_MATURATION = 1.05;               % factor to increase the ISI_END_BURST as the burst has more units (default: 1.05 -> 5% increase per unit, starting from MIN_SPIKE_BURST)
        ISI_END_BURST_MAX = 500;                % maximum threshold (maximum) ISI to identify the end of a burst, in ms (500 -> 2Hz) 
        MIN_IBI = 500;                          % minimum inter-burst interval to identify two consecutive bursts, in ms (default: 500); from finish to start of the next one
        MIN_BURST_DURATION = 25;                % minimum burst duration, in ms (default: 25 -> with 3 spikes, it means 80 Hz)
        apply_parameters('BurstDetection', 'method1');                      % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")
        
        parameters = struct('method','Maximal Interval','RATIO_PEAK',RATIO_PEAK,'MIN_SPIKE_BURST',MIN_SPIKE_BURST,'ISI_START_BURST',ISI_START_BURST,'ISI_END_BURST_INI',ISI_END_BURST_INI,'MIN_IBI',MIN_IBI ...
            ,'FACTOR_MATURATION',FACTOR_MATURATION,'ISI_END_BURST_MAX',ISI_END_BURST_MAX,'MIN_BURST_DURATION',MIN_BURST_DURATION,'ISI_PEAK_BURST',ISI_PEAK_BURST);
        BurstDetection.parameter = parameters;

    case 2
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the inter-spike interval

signal_ISI = diff(signalSpike);                                             % inter-spike interval, in msec
n_spikes = size(signalSpike,1);                                             % total number of events in the spike train                        
if size(signalSpike,2) > 1                                                  % if the variable signalSpike is not a column array
    error("The input variable for *BurstDetection* has the wrong dimensions (must be a column array")
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Burst Detection using the method "option_method"

switch option_method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1          % Maximal Interval

    index_iniPairs = find(signal_ISI  < ISI_START_BURST);                   % index of all pairs of spikes close enough to be the beginning of a burst
    n_potential = numel(index_iniPairs);                                    % number of potential bursts (the criteria for the beginning of the burst was met)
    potentialBurst = nan(n_potential,2);                                    % variable to save the initial and last point of each potential burst event
    n_burst = 0;                                                            % total number of burst events (initialized at zero)

    for i=1:n_potential                                                     % across all potential starting points of a burst (if there is no potential burst, nothing is done: n_burst = 0)
        start_idx = index_iniPairs(i);                                      % index of the potential starting point of the burst i
        if n_burst>0 && start_idx <= potentialBurst(n_burst,2)              % the starting point is contained in the previous burst -> next potential starting point
            continue
        end
        end_idx = start_idx + MIN_SPIKE_BURST - 1;                          % index of the potential end point of the burst (the burst can not be shorter than this)
        if end_idx > n_spikes                                               % the potential burst goes beyond the recording length (and so do all the next ones, since they start later)
            break
        end
        if any(signal_ISI(start_idx:(end_idx-1)) > ISI_END_BURST_INI)       % there is an isi larger than ISI_END_BURST_INI among the first MIN_SPIKE_BURST spikes: the burst would be cut
            continue                                                        % with less than MIN_SPIKE_BURST spikes, so it is always rejected (same result as before, where the burst was cut and then rejected)
        end

        threshold_isi = ISI_END_BURST_INI*FACTOR_MATURATION;                % updated threshold isi to identify the end of the burst
        while end_idx<n_spikes && signal_ISI(end_idx)<threshold_isi                 % while the last spike of the burst has a subthreshold isi
            threshold_isi = min(threshold_isi*FACTOR_MATURATION, ISI_END_BURST_MAX);% update threshold isi to identify the end of the burst
            end_idx = end_idx+1;                                                    % update the end of the burst
        end                                                                 % (the "clear threshold_isi" that was here was the slowest line of the function)

        % required features for a burst
        duration = signalSpike(end_idx) - signalSpike(start_idx);           % duration of the burst, in msec
        isi_burst = signal_ISI(start_idx:(end_idx-1));                      % vector with the isi values for the burst, in msec
        if numel(isi_burst)+1 >= MIN_SPIKE_BURST && ...                     % if the spike train meets the minimum number of spikes
           duration > MIN_BURST_DURATION && ...                             % and the minimum duration criteria
           min(isi_burst) < ISI_PEAK_BURST                                  % and the minimum peak spiking rate criteria (at least a clear peak in firing rate)
            n_burst = n_burst+1;                                            % update the total number of burst detected
            potentialBurst(n_burst,1) = start_idx;                          % index of the first spike of the burst
            potentialBurst(n_burst,2) = end_idx;                            % index of the last spike of the burst
        end
    end

    % Merge the bursts that are too close 
    if n_burst>1                                                            % at least two potential bursts                                                              
        for i=n_burst-1:-1:1                                                % across all potential burst but the last one, in reverse order
            t_ini = signalSpike(potentialBurst(i+1,1));                     % first spike of the next potential burst (i+1), in msec
            t_fin = signalSpike(potentialBurst(i,2));                       % last spike of the current potential burst (i), in msec
            ibi = t_ini - t_fin;                                            % gap between the bursts i and i+1, in msec

            duration_post = signalSpike(potentialBurst(i+1,2))-signalSpike(potentialBurst(i+1,1));  % duration of the next potential burst (i+1)  
            duration_pre = signalSpike(potentialBurst(i,2))-signalSpike(potentialBurst(i,1));       % duration of the current potential burst (i)
            ibi_threshold = max(min(duration_pre, duration_post),MIN_IBI);  % in case the two potential bursts are quite long (compared to MIN_IBI)

            if ibi < ibi_threshold                                          % if the bursts are too close
               potentialBurst(i,2) = potentialBurst(i+1,2);                 % merge the burst i and i+1
               potentialBurst(i+1,:) = [];                                  % eliminate the burst i+1 (as independent burst)
               n_burst = n_burst - 1;                                       % update the total number of bursts
            end
        end
    end
    aux = zeros(n_spikes,3);
    aux(:,1) = signalSpike;                                                 % time stamp for the spike train
    aux(:,2) = [signal_ISI; nan];                                           % isi values

    % Identify the peak of the burst and the median/middle spike-time
    if n_burst > 0                                                          % pre-allocate the structure array (same fields and order as before)
        BurstDetection.IndivBurst = struct('peak', cell(1,n_burst), 'start', [], 'last', [], 'isi', [], 'median', []);
    end
    for i=1:n_burst                                                         % across all identified bursts
        ini = potentialBurst(i,1);                                          % index of the first spike of the burst i
        fin = potentialBurst(i,2);                                          % index of the last spike of the burst i
        n_spike_burst = fin-ini+1;                                          % total number of spikes in the burst i
        n_peak = max(MIN_SPIKE_BURST, ceil(RATIO_PEAK*n_spike_burst));      % number of spikes that are grouped to look for the peak of the burst (n_peak <= n_spike_burst)
        shift = floor(n_peak/2);

        aux(ini:fin,3) = i;                                                 % load the burst number (i) for all the spikes of the burst i
        bin = signalSpike(ini+n_peak-1:fin) - signalSpike(ini:fin-n_peak+1);% temporal distance between "n_peak" consecutive spikes (one value per group of spikes within the burst)
        [~, idx_min] = min(bin);                                            % search for the peak frequency within the burst i (shortest distance; "nanmin" is not needed, there are no NaN)
        BurstDetection.IndivBurst(i).peak = signalSpike(idx_min+shift+ini-1);   % time stamp of the peak of burst i
        BurstDetection.IndivBurst(i).start = signalSpike(ini);                  % time stamp of the first spike of burst i, in ms
        BurstDetection.IndivBurst(i).last = signalSpike(fin);                   % time stamp of the last spike of burst i, in ms
        BurstDetection.IndivBurst(i).isi = signal_ISI(ini:fin-1);               % isi corresponding to all spikes of the burst i
        BurstDetection.IndivBurst(i).median = median(signalSpike(ini:fin));     % median/middle spike-time, in ms
    end
    
    BurstDetection.spikes = aux;                                            % raw data
    BurstDetection.nBurst = n_burst;                                        % total number of burst events

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2          %
end

end