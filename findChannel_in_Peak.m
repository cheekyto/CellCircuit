function [Ch_in_Peak, aux_burst, burst_id, meaningful_activity] = findChannel_in_Peak(BurstDetected, n_channels, t_step, Npeak, PopActivity, range, opt_test, option)

% function to identify the channels firing during a population peak
% (revised October 2026: the spikes within all the events are counted at once -see "windowCounts"-, same results; no "clear" inside the loops)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   structure array BurstDetected 
%   n_channels:                         total number of channels
%   NPeak:                              total number of population peaks
%   PopActivity.IndivBurst(j).start:    start time of the population peak j, in sec
%   PopActivity.IndivBurst(j).end:      end time of the population peak j, in sec
%   range:                              range of the recording (range(1)=start; range(2)=end), in sec
%   option = 1                          Genereal screening 
%          = 2                          Statistical assessement

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Ch_in_Peak:                                         structure "PopActivity" updated with new fields
%   aux_burst(i):                                       total number of bursting neurons within burst i
%   burst_id(i,ch):                                     =0, if the neuron ch does not fire within the event i 
%                                                       =1, if the neuron ch fires during the event i
%                                                       =2, if the neuron ch has a burst during the event i
%                                                       =3, if the neuron ch has multiple spikes during the event i (but they are not part of a burst)
%   meaningful_activity(ch):                            =1, if the activity during the peaks is higher than the rest of the recording

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SAFE_K = 0.40;                                  % when saving the spikes for each network burst event, the extra time is added on each side (default: 0.4)                          
ACTIVITY_THRESHOLD = 3.0;                       % threshold for the ratio between the firing rate within potential network bursts and overall firing rate, for individual neurons (default: 3.0)
K_ONSET = 0.20;                                 % when assessing the statistical significance of network burst participation, the extra time is added at the begining (default: 0.2)

MIN_MULTIPLE = 2;                               % minimum number of spikes a unit has to display during a network burst to be considered a multiple unit (default: 2) 

% for satistical assesment 
N_SHUFFLING = 500;                              % number of shuffling repetitions for statistical assessment (default: 500)
RNG_SEED = 0;                                   % base seed for the shuffling; channel ch uses RNG_SEED + ch (change it to get a different, but still reproducible, set of shuffles)
STAT_THRESHOLD = 0.95;                          % threshold of statistical significance (default: 0.95)
MIN_PARTICIPATION = 0.20;                       % minimum proportion of network bursts in which a channel has to participate (default: 0.20)
MIN_PARTICIPATION_ABS = 2.0;                    % minimum number of network bursts in which a channel has to participate (default: 2)
ACTIVITY_THRESHOLD_FINAL = 1.80;                % (used during statistical assessment) threshold for the ratio between the firing rate within potential network bursts and overall firing rate (default: 2.0, NOT higher)
                                                % lower it (1.60) when there is high level of population activity (it is not a few individual units with high level of activity)
apply_parameters('findChannel_in_Peak');                                    % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

K_MSEC_SEC = 1/1000;                            % factor to rescale from msec to sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch option

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    General screening 

    case 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Identify the channels involved in each population peak

burst_id = zeros(Npeak, n_channels);                                        % =0, if the neuron does not fire within the event 
                                                                            % =1, if the neuron fires during the event
                                                                            % =2, if the neuron has a burst during the event
burst_labels = cell(Npeak, n_channels);                                     % labels of the bursts for each cluster that participate from the different peaks
timeStamp = cell(Npeak, n_channels);                                        % time stamps of the spikes from each neurons that happened during the different peaks

ini_Nburst = zeros(Npeak,1);    fin_Nburst = zeros(Npeak,1);                % (revised October 2026) windows of all events, calculated once (before, once per cluster)
for k=1:Npeak                                                               % across all events 
    ini_Nburst(k) = PopActivity.IndivBurst(k).start/K_MSEC_SEC;             % starting point of the event k, in msec
    fin_Nburst(k) = PopActivity.IndivBurst(k).last/K_MSEC_SEC;              % last point of the event k, in msec
    time_added = (fin_Nburst(k)-ini_Nburst(k))*SAFE_K;                      % extra time added on both ends
    ini_Nburst(k) = ini_Nburst(k) - time_added;                             % just in case, add more data points to the left
    fin_Nburst(k) = fin_Nburst(k) + time_added;                             % just in case, add more data points to the right
    if n_channels > 0
        PopActivity.IndivBurst(k).window = fin_Nburst(k) - ini_Nburst(k);   % time window for the spike train extraction (based on the core duration of the burst)
    end
end

for ch=1:n_channels                                                         % across all clusters
            
    switch opt_test
        case 1  % old method

            for k=1:Npeak                                                           % across all events 
                aux1 = find(BurstDetected{ch}.spikes(:,1)>=ini_Nburst(k));          % index for all spikes of neuron ch, after starting point of event k
                aux2 = find(BurstDetected{ch}.spikes(:,1)<fin_Nburst(k));           % index for all spikes of neuron ch, before last point of event k
                aux = intersect(aux1, aux2);                                        % index for all spikes of neuron ch, within the event k

                if numel(aux)>0
                    burst_id(k,ch)=1;                                               % label the channels that have at least one spike within the event k    
                    timeStamp{k,ch} = BurstDetected{ch}.spikes(aux,1);              % time stamp of all spikes from neuron ch within the event k
                    burst_label = BurstDetected{ch}.spikes(aux,3);                  % for all spikes of neuron ch that are within event k, # of burst to which each spike belongs (=0 if the spike is NOT part of a burst)
                    if sum(burst_label)>0                                           % if at least one spike within event k belong to an individual burst of neuron ch
                        burst_id(k,ch) = 2;                                         % label the channels that have at least one spike withing the event k that also belongs to an individual burst of the neuron
                    else
                        if numel(aux) >= MIN_MULTIPLE                               % the unit ch displayed multiple spikes during burst k
                            burst_id(k,ch) = 3;
                        end
                    end
                    burst_labels{k, ch} = unique(burst_label);                      % label of individual bursts of cluster ch that happens within network burst k
                else
                    timeStamp{k,ch} = [];                                           % do not have spikes within the network burst
                end
            end
            
        case 2  % binnary search (revised October 2026: all the events at once with "windowCounts" -same result as "windowMembership" event 
                %  by event-, and no "clear" inside the loops)
            
            [n_hits, first] = windowCounts(BurstDetected{ch}.spikes(:,1), ini_Nburst, fin_Nburst);    % number of spikes of neuron ch within each event, and index of the first one
            for k=find(n_hits(:)' > 0)                                              % across the events with at least one spike (the others keep timeStamp{k,ch} = [])
                index = first(k):1:(first(k)+n_hits(k)-1);                          % index of the spikes of neuron ch within the event k
                burst_id(k,ch) = 1;                                                 % label the channels that have at least one spike within the event k 
                timeStamp{k,ch} = BurstDetected{ch}.spikes(index,1);                % time stamp of all spikes from neuron ch within the event k
                burst_label = BurstDetected{ch}.spikes(index,3);                    % for all spikes of neuron ch that are within event k, # of burst to which each spike belongs (=0 if the spike is NOT part of a burst)
                if any(burst_label > 0)                                             % if at least one spike within event k belong to an individual burst of neuron ch
                    burst_id(k,ch) = 2;                                             % label the channels that have at least one spike withing the event k that also belongs to an individual burst of the neuron
                else
                    if n_hits(k) >= MIN_MULTIPLE                                    % the unit ch displayed multiple spikes during burst k
                        burst_id(k,ch) = 3;
                    end
                end
                burst_labels{k,ch} = unique(burst_label);                           % label of individual bursts of cluster ch that happens within network burst k
            end
    end     % end for the swtich based on opt_test
end     % end of the for loop across channels
        
PopActivity.timeStamp = timeStamp;                                          % time stamp of the spikes within the network bursts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify the neurons that are too noisy (activity during network burst events is similar to the rest)
        
meaningful_activity = zeros(n_channels,1);
for ch=1:n_channels                                                         % across all clusters
            
    spikes_in = 0;                                                          % total number of spikes during the potential network bursts, for channel ch
    window_in = 0;                                                          % total duration of all potential network burst, in sec
            
    for k=1:Npeak
        spikes_in = spikes_in + numel(PopActivity.timeStamp{k,ch});         
        window_in = window_in + PopActivity.IndivBurst(k).window*K_MSEC_SEC;  % extended duration of the event k, in sec: the same window used to collect its spikes
                                                                            % (revised October 2026: it was the core duration, last-start, while the spikes were counted
                                                                            %  in the window extended by SAFE_K on each side -> rate overestimated by up to 1.8x)
    end

    switch opt_test
            
        case 1  % old method

            aux1 = find(BurstDetected{ch}.spikes(:,1)>=range(1)/K_MSEC_SEC);         % index for all spikes of neuron ch, after starting point of analysis window
            aux2 = find(BurstDetected{ch}.spikes(:,1)<range(2)/K_MSEC_SEC);         % index for all spikes of neuron ch, before last point of analysis window
            aux = intersect(aux1, aux2);                                            % index for all spikes of neuron ch, within the event k
            spikes_total = numel(aux);                                              % total number of spikes of neuron ch during the analysis window 
    
        case 2  % binnary search
                
            [spikes_total, ~, ~] = windowMembership(BurstDetected{ch}.spikes(:,1), range(1)/K_MSEC_SEC, range(2)/K_MSEC_SEC, 1, 0);
    end

    freq_in = spikes_in/window_in;
    freq_total = spikes_total/(range(2)-range(1));
    if freq_in > freq_total * ACTIVITY_THRESHOLD
        meaningful_activity(ch) = 1;                                        % the activity during the peaks is higher than the rest of the recording
    else
        burst_id(:,ch) = 0;                                                 % remove channel ch from the participating clusters
    end
end

if Npeak>0
    aux_burst = nan(Npeak,1);
    aux_multiple = nan(Npeak,1);
else
    aux_burst = nan;
    aux_multiple = nan;
end

for k=1:Npeak                                                               % across all events 
    aux1 = find(burst_id(k,:)>0);                                           % index of neurons that fire within burst k                  
    aux2 = find(burst_id(k,:)==2);                                          % index of neurons with individual burst within network burst k    
    aux3 = find(burst_id(k,:)>1);                                           % index of neurons with multiple spikes within burst k
    aux_burst(k) = numel(aux2);                                             % total number of bursting neurons within burst k (to use later)
    aux_multiple(k) = numel(aux3);                                          % total number of neuron with multiple spikes within burst k

    PopActivity.IndivBurst(k).nChannel.total = numel(aux1);
    PopActivity.IndivBurst(k).nChannel.bursting = numel(aux2);
    PopActivity.IndivBurst(k).nChannel.multiple = numel(aux3);
    PopActivity.IndivBurst(k).channels.total(:,1) = aux1;                   % index of all neurons involved in event k
    PopActivity.IndivBurst(k).channels.bursting(:,1) = aux2;                % index of bursting neurons involved in event k
    PopActivity.IndivBurst(k).channels.multiple(:,1) = aux3;                % index of neurons with multiple spikes involved in event k
            
    clear aux1 aux2
end

Ch_in_Peak = PopActivity;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Statistical assessement
    
    case 2

        particip_ind = nan(n_channels,1);
        event_ind = nan(n_channels,1); 
        freq_ratio_ind = nan(n_channels,1);

score = nan(n_channels,2);                                                  % matrix to collect the information about the statistical assessment
if opt_test == 1                                                            % (revised October 2026) with opt_test = 2, the bursts are characterized inside "shuffleCache",
    for i=1:n_channels                                                      % only when the shuffles of the cluster are generated
        BurstFeature{i} = BurstFeatures(BurstDetected{i});                  % characterize the bursts for cluster i
    end
end

window = nan(Npeak,2);                                                      % temporal window for each burst event
window2 = nan(Npeak,2);                                                     % auxiliary temporal window for each burst event        
for k=1:Npeak                                                               % across all network bursts
    ini0 = PopActivity.IndivBurst(k).start;
    fin0 = PopActivity.IndivBurst(k).last;
    window(k,1) = (ini0 - (fin0-ini0)*K_ONSET)/K_MSEC_SEC;                  % starting point of the burst k, in msec
    window(k,2) = (fin0 + (fin0-ini0)*K_ONSET)/K_MSEC_SEC;                  % last point of the burst k, in msec
    window2(k,1) = (ini0 - (fin0-ini0)*SAFE_K)/K_MSEC_SEC;                  % starting point of the burst k, in msec
    window2(k,2) = (fin0 + (fin0-ini0)*SAFE_K)/K_MSEC_SEC;                  % last point of the burst k, in msec
end
        
% Identify the participating clusters
participation = zeros(n_channels,1);
for ch=1:n_channels                                                         % across all clusters
    event = 0;
    spikes_in = 0;                                                          % total number of spikes during the potential network bursts
    window_in = 0;                                                          % total duration of all potential network burst, in sec

    data = BurstDetected{ch}.spikes(:,1);                                   % spike train data for channel ch
    if numel(data)>0                                                        % there is at least one unit
        hits = zeros(N_SHUFFLING,1);
                
        switch opt_test
            case 1  % old method
                ShuffleSpike = ShuffleSpikeTrain(data, range/K_MSEC_SEC, t_step, 3, N_SHUFFLING, 0, BurstDetected{ch}, BurstFeature{ch}, RNG_SEED + ch);
                for i=1:N_SHUFFLING                                         % across all shuffling iterations
                    signal = ShuffleSpike{i}(:,1);                          % shuffled spike train, in msec
                    for k=1:Npeak                                           % across all burst events
                        aux1 = find(signal>=window(k,1));
                        aux2 = find(signal<window(k,2));
                        aux = intersect(aux1, aux2);                        % index of spikes within the burst k
                        if numel(aux)>0                                     % there is at least one unit within the burst k
                            hits(i) = hits(i) + 1;                          % number of bursts in which the shuffled signal has a participating unit
                        end
                    end
                end

            case 2  % binnary search (revised October 2026: all the burst events at once with "windowCounts" -same result as "windowMembership"
                    %  event by event-; before, a "clear" inside the innermost loop took ~half of the total time of "Analyze_spikeTrain")
                    % (revised October 2026) the shuffles of each cluster are the same every time this function is called (same seed per cluster), so they
                    % are generated only once per run of "NetworkBurstDetection" (see "shuffleCache"), and all of them are counted at once ("shuffleHits");
                    % same result as before: for each shuffle i, hits(i) = sum(windowCounts(ShuffleSpike{i}(:,1), window(:,1), window(:,2)) > 0)
                shuffles = shuffleCache('get', ch, data, range/K_MSEC_SEC, t_step, N_SHUFFLING, BurstDetected{ch}, RNG_SEED + ch);
                hits = shuffleHits(shuffles, window(:,1), window(:,2));     % number of bursts in which each shuffled signal has a participating unit
        end
                
        switch opt_test
            case 1  % old method
                for k=1:Npeak                                               % across all burst events
                    aux1 = find(data>=window(k,1));
                    aux2 = find(data<window(k,2));
                    aux = intersect(aux1, aux2); 
                    if numel(aux)>0
                        event = event + 1;                                  % number of bursts in which the original signal has a participating unit 
                    end
                    aux3 = find(data>window2(k,1));
                    aux4 = find(data<window2(k,2));
                    auxbis = intersect(aux3, aux4); 
                    spikes_in = spikes_in + numel(auxbis);
                end
                                        
            case 2  % binnary search (revised October 2026: all the events at once, as above)
                n_events1 = windowCounts(data, window(:,1), window(:,2));  % number of spikes of the original signal within each burst event
                event = sum(n_events1 > 0);                                 % number of bursts in which the original signal has a participating unit 
                n_events2 = windowCounts(data, window2(:,1), window2(:,2));% number of spikes within each (extended) burst event
                spikes_in = sum(n_events2);                                 % total number of spikes during the potential network bursts
        end
        for k=1:Npeak                                                       % across all burst events
            window_in = window_in + (window2(k,2)-window2(k,1))*K_MSEC_SEC; % in sec
        end
        
        score(ch,2) = event;                                                % number of bursts in which the original signal from cluster ch has a participating unit 
        score(ch,1) = numel(find(hits<event))/N_SHUFFLING;                  % normalized number of times in which the shuffled signal has fewer participating units than the original spike train
    else                                                                    % There are no spikes in cluster ch
        score(ch,2) = 0;
        score(ch,1) = 0;
    end      

        
    switch opt_test
        case 1  % old method
    
            aux1 = find(data>=range(1)/K_MSEC_SEC);                                  % index for all spikes of neuron ch, after starting point of analysis window
            aux2 = find(data<range(2)/K_MSEC_SEC);                                  % index for all spikes of neuron ch, before last point of analysis window
            aux = intersect(aux1, aux2);                                            % index for all spikes of neuron ch, within the event k
            spikes_total = numel(aux);                                              % total number of spikes of neuron ch during the analysis window 
    
        case 2  % binnary search
                
            [spikes_total, ~, ~] = windowMembership(data, range(1)/K_MSEC_SEC, range(2)/K_MSEC_SEC, 1, 0);
    end

    freq_in = spikes_in/window_in;                                          % to compare with the general spiking rate
    freq_total = spikes_total/(range(2)-range(1));

    particip_ind(ch) = event/Npeak;
    event_ind(ch) = event; 
    freq_ratio_ind(ch) = freq_in/freq_total;

    if score(ch,1) >= STAT_THRESHOLD
        particip = event/Npeak;
        freq_ratio = freq_in/freq_total;
        if particip > MIN_PARTICIPATION && event >= MIN_PARTICIPATION_ABS && freq_ratio > ACTIVITY_THRESHOLD_FINAL
            participation(ch) = 1;                                          % the participation of cluster ch in network bursting is statistically signicant 
        end
    end
end

Ch_in_Peak = PopActivity;
Ch_in_Peak.participation = participation;
aux_burst = [];     burst_id = [];      meaningful_activity = [];
end
end