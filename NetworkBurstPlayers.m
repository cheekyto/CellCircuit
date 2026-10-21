function NetworkBurstPlayers = NetworkBurstPlayers(spikeTimestamp, range, PopulationActivityDetected, BurstDetected, option, bin, threshold, opt_distance, option_plot)

% function to identify and classify the units involved in network bursting
% (spikeTimestamp, PopulationActivityDetected, and BurstDetected should be restricted to the same temporal window)

% the option = 1 runs a first screening, classifying the potential events into WEAK and SOLID. It also identify the robustly pariticipating clusters

option_test = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   structure array spikeTimestamp (its length is the number of recording channels), with the following fields
% 	.recLength:	length of the recording, in seconds.
% 	.nUnits:	total number of identified units in the channel
% 	.timeStamp:	time stamp for all units (nUnits) in the channel, in msec
% 	.position:	position of the electrode

%   PopulationActivityDetected.spikesBin:            binned population spike
%   PopulationActivityDetected.time:                 midpoint of the bin i, in sec
%   PopulationActivityDetected.nBurst:               total number of network bursts 
%   PopulationActivityDetected.IndivBurst(j).peak:   peak point of network burst j, in sec       
%   PopulationActivityDetected.IndivBurst(j).start:  start point of network burst j, in sec
%   PopulationActivityDetected.IndivBurst(j).last:   last point of network burst j, in sec
%   PopulationActivityDetected.IndivBurst(j).peakIdx:peak point of network burst j
%   PopulationActivityDetected.IndivBurst(j).level:  level of network burst j

%   range:      range of the recording (range(1)=start; range(2)=end), in sec
%   bin:        size of the temporal bin to calculate the population firing rate, in sec (0.2-0.4)        
%   threshold:  in terms of the standard deviation (ventral: 3.0 - 5.0 ; dorsal: 5.0 - 8.0)
%   option = 1          First screening: robust vs weak network bursts 
%          = 2          Statistical assessement
%          = 3          Only participating channels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   NetworkBurstPlayers.IndivBurst(j).peak:                     peak point of network burst j, in sec
%   NetworkBurstPlayers.IndivBurst(j).peakAmp:                  amplitude of the peak of network burst j
%   NetworkBurstPlayers.IndivBurst(j).start:                    start point of network burst j, in sec
%   NetworkBurstPlayers.IndivBurst(j).last:                     last point of network burst j, in sec
%   NetworkBurstPlayers.IndivBurst(j).peakIdx:                  peak point of network burst j
%   NetworkBurstPlayers.IndivBurst(j).window:                   temporal window to extract spike train from network burst j (revised October 2026: in ms after
%                                                               the screening -option 1, SAFE_K on each side-; in sec after the statistical assessment -option 2,
%                                                               SAFE_K*duration before and 2*SAFE_K*duration after the event, intended-)
%   NetworkBurstPlayers.IndivBurst(j).nChannel.total:           number of neurons involved in network burst j 
%   NetworkBurstPlayers.IndivBurst(j).nChannel.bursting:        number of bursting neurons involved in network burst j 
%   NetworkBurstPlayers.IndivBurst(j).channels.total(i,1):      list of neurons involved in network burst j
%   NetworkBurstPlayers.IndivBurst(j).channels.bursting(i,1):   list of bursting neurons involved in network burst j
%   NetworkBurstPlayers.IndivBurst(j).solid:                    =1, if it is a solid burst; =0, weak burst (only across a few number of channels)
%   NetworkBurstDetection.nBurst(1):                            total number of network bursts (weak + solid)
%                               (2):                            total number of solid network bursts
%   NetworkBurstPlayers.Channels(i,1):                          list of channels that robustly participate of network bursting     
%   NetworkBurstPlayers.timeStamp{i,j}:                         time stamp of the spikes for the burst i coming from channel j (all channels), in msec

%%%%  same field for NetworkBurstPlayers.Secondary
%   NetworkBurstPlayers.Secondary.ch_parent{j}                  list of all primary peaks the neuron j is associated to

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SAFE_K = 0.40;                                  % when saving the spikes for each network burst event, the extra time is added on each side (default: 0.4)                          
                                                % (revised October 2026: in the statistical assessment, SAFE_K*duration before the event and 2*SAFE_K*duration after it, intended)

% to classify the events SOLID vs WEAK
MIN_CHAN_BURST =  5;                            % minimum number of channels involved in a single network burst (default: 5)
MIN_CHAN_BURST_RATIO = 0.05;                    % minimum proportion of bursting channels involved in a single network burst, relative to the median of bursting participating units (default: 0.4)
MIN_MULTIPLE = 3;                               % minimum number of channels with multiple spikes involved in a single network burst (default: 3)
MIN_MULTIPLE_RATIO = 0.03;                      % minimum proportion of channels with multiple spikes involved in a single network burst, relative to the total number of units (default: 0.03)

% to classify the units
MIN_PARTICIPATION_TONIC = 0.65;                 % proportion of solid network burst in which a non-bursting channel has to participated to be counted (default: 0.95)
MIN_PARTICIPATION_PHASIC = 0.25;                % proportion of solid network burst in which a bursting channel has to participated to be counted (default: 0.6)
MIN_SOLID_BURST = 3;                            % minimum number of solid burst to classify the channels as solid vs weak (default: 3)

% to assess the secondary peaks
MIN_PEAK2_GAP = 1;                              % minimum number of secondary peaks per gap/window (default: 1)
GAP_WITH_PEAK = 0.75;                           % minimum proportion of primary peaks followed by secondary peaks (default: 0.75)                    

MIN_PARTICIPATION_SECONDARY = 0.40;             % minimum proportion of secondary network bursts in which a channel has to participate (default: 0.40)        

K_MSEC_SEC = 1/1000;                            % factor to rescale from msec to sec
N_PLOT = 20;
SHIFT = 1;

FAST_MODE = 1;                                  % =1, don't search for participating units for secondary peaks until the last round
                                                % =0, seach for participating units for secondary peaks
apply_parameters('NetworkBurstPlayers');                                    % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Prepare the data and the settings       %%%

n_channels = size(spikeTimestamp,2);                                        % total number of clusters
t_step = spikeTimestamp(1).t_step;
peak1 = zeros(PopulationActivityDetected.nBurst,1);
for i=1:PopulationActivityDetected.nBurst
    peak1(i) = PopulationActivityDetected.IndivBurst(i).peakIdx;            % index of the primary population peaks
end
Npeak1 = PopulationActivityDetected.nBurst;                                 % total number of primary peaks
Npeak1_level1 = PopulationActivityDetected.nBurst_level1;                   % total number of level1 primary peaks 

if isnan(PopulationActivityDetected.Secondary.nBurst) || PopulationActivityDetected.Secondary.nBurst == 0
    Npeak2 = 0;
else
    Npeak2 = PopulationActivityDetected.Secondary.nBurst;                   % total number of secondary peaks
end
if Npeak1 < 1 || isnan(Npeak1)                                              % (revised October 2026: after the block above, which overwrote it; with no primary
    Npeak2 = 0;                                                             %  events but some secondary events, "participation = [] + participation2" crashed)
end

if option == 2
    FAST_MODE = 0;
end
opt_shorcut = 0;
if option == 3
    option = 2;
    opt_shorcut = 1;
    FAST_MODE = 1;
end

switch option

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Initial screening: robust vs weak network bursts 
    
    case 1

        NetworkBurstPlayers = PopulationActivityDetected;                   % load the information for the population peaks
        meaningful_activity = [];                                           % (revised October 2026) defined also without primary events (used by the plot)
        if Npeak1 > 0    
            [NetworkBurstPlayers, n_bursting1, burst_id1, meaningful_activity] = findChannel_in_Peak(BurstDetected, n_channels, t_step, Npeak1, NetworkBurstPlayers, range, option_test, option);
        end
        if Npeak2 > 0 && ~FAST_MODE
            [NetworkBurstPlayers.Secondary, n_bursting2, burst_id2, ~] = findChannel_in_Peak(BurstDetected, n_channels, t_step, Npeak2, NetworkBurstPlayers.Secondary, range, option_test, option);
        end

        %% Classify the primary network bursts based on the number of bursting neurons participating on it
        
        if Npeak1 > 0
            % Assess overall number of participating units to classify the peaks
            threshold = max(MIN_CHAN_BURST, median(n_bursting1)*MIN_CHAN_BURST_RATIO);  % minimum number of neurons to classify as a solid network burst
            aux_weak = find(n_bursting1<threshold);                                     % index of events that do not meet the threshold; weak events (potential false positive)
            threshold_multiple = max(MIN_MULTIPLE, n_channels*MIN_MULTIPLE_RATIO);      % minimum number of neurons with multiple spike to classify as a solid network burst
            for k=1:Npeak1
                if ismember(k,aux_weak)
                    NetworkBurstPlayers.IndivBurst(k).solid = 0;                        % =0, label for weak events    
                else
                    NetworkBurstPlayers.IndivBurst(k).solid = 1;                        % =1, label for solid/strong events
                end
                if NetworkBurstPlayers.IndivBurst(k).nChannel.multiple < threshold_multiple
                    NetworkBurstPlayers.IndivBurst(k).solid = 0;                        % =0, label for weak events
                    fprintf("\n primary event discarded for multiple spikes")
                end
            end
            NetworkBurstPlayers.nBurst(2) = sum([NetworkBurstPlayers.IndivBurst(1:Npeak1).solid]);  % number of solid events (revised October 2026: counted from the labels,
                                                                                    % so both rules are included; before, the events downgraded for too few units with
                                                                                    % multiple spikes were still counted as solid)

        %% Identify the neurons that robustly participate in solid primary network bursts

            participating_channel1 = zeros(n_channels,1);
            n_solid = max(NetworkBurstPlayers.nBurst(2), MIN_SOLID_BURST);              % total number of solid primary network burst events
            solid_event = [NetworkBurstPlayers.IndivBurst(1:Npeak1).solid]' == 1;       % (revised October 2026) participation is counted only in the solid events, for both rules
                                                                                        % (the thresholds are proportions of the solid events; before, weak events also counted)
            for ch=1:n_channels                                                         % across all clusters
                aux1 = find((burst_id1(:,ch)==1 | burst_id1(:,ch)==3) & solid_event);   % index of solid burst events where the neuron ch fires without bursting: one spike (=1) or several spikes (=3)
                                                                                        % (revised October 2026: the events with several spikes but no burst, =3, were not counted)
                aux2 = find(burst_id1(:,ch)==2 & solid_event);                          % index of solid burst events where the neuron ch displays individual bursts
                if numel(aux2) > numel(aux1)                                            % bursting neuron
                    if numel(aux2)>n_solid*MIN_PARTICIPATION_PHASIC                     % enough bursting during network events
                        participating_channel1(ch) = 1;                                 % label the neuron as robust participant
                    end
                else                                                                    % tonic neuron
                    if numel(aux1)+numel(aux2)>n_solid*MIN_PARTICIPATION_TONIC           % enough activity during network events (revised October 2026: per-unit rule restored;
                                                                                        %  the line used "IndivBurst(k).nChannel.total", where k was left over from a previous loop,
                                                                                        %  so every tonic unit got the same decision, based on the last event)
                        participating_channel1(ch) = 1;                                 % label the neuron as robust participant
                    end
                end
            end
            NetworkBurstPlayers.Channels(:,1) = find(participating_channel1);           % index of all robustly participating neurons
        else
            NetworkBurstPlayers.nBurst(2) = 0;
            NetworkBurstPlayers.Channels = []; 
        end

        %% Identify the neurons that robustly participate in secondary network bursts
        % NOTE (October 2026): this block does not run (FAST_MODE = 1 for option 1). The rule MIN_PEAK2_GAP / GAP_WITH_PEAK is applied
        % in "NetworkBurstDetection", after the final quality control of the secondary events
        if Npeak2>0 && ~FAST_MODE

            aux  = find(PopulationActivityDetected.Secondary.nBurst_per_gap >= MIN_PEAK2_GAP);
%            if numel(aux) > Npeak1*GAP_WITH_PEAK                            % enough number of gaps with secondary peaks
            if numel(aux) > Npeak1_level1*GAP_WITH_PEAK                     % enough number of gaps with secondary peaks
                participating_channel2 = zeros(n_channels,1);
                for ch=1:n_channels                                         % across all clusters
                    aux2 = find(burst_id2(:,ch)>0);                         % index of burst events where the neuron ch fires
                    parent = [];
                    if numel(aux2)>0
                        for j=1:numel(aux2)                                 % across all secondary peaks where neuron ch fires
                            parent = vertcat(parent, PopulationActivityDetected.Secondary.IndivBurst((aux2(j))).parent);
                        end
                        parent = unique(parent);
%                        total_parent = numel(parent)/Npeak1;                % number of different gaps (defined by primary peaks) where the neuron ch fires    
                        total_parent = numel(parent)/Npeak1_level1;         % number of different gaps (defined by primary peaks) where the neuron ch fires    
                        if total_parent > MIN_PARTICIPATION_SECONDARY
                            participating_channel2(ch) = 1;                 % label the neuron as participant
                        end
                    end
                    ch_parent{ch,1} = parent;
                    clear parent aux2
                end                
            else
                ch_parent = [];
                participating_channel2 = [];
            end
            
            if numel(participating_channel2) > 0
                NetworkBurstPlayers.Secondary.Channels(:,1) = find(participating_channel2);    % index of all participating neurons in secondary peaks
                NetworkBurstPlayers.Secondary.ch_parent = ch_parent;                           % list of gaps (defined by primary peaks) where the neuron ch fires
            else
                NetworkBurstPlayers.Secondary.Channels = [];    % index of all participating neurons in secondary peaks
                NetworkBurstPlayers.Secondary.ch_parent = [];                           % list of gaps (defined by primary peaks) where the neuron ch fires
            end
        else
            NetworkBurstPlayers.Secondary.Channels = [];
        end
                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Statistical assessement
    
    case 2

        NetworkBurstPlayers = PopulationActivityDetected;                   % load the information for the population peaks
        if Npeak1 > 0    
            [NetworkBurstPlayers, a, b, c] = findChannel_in_Peak(BurstDetected, n_channels, t_step, Npeak1, NetworkBurstPlayers, range, option_test, option);
            participation1 = NetworkBurstPlayers.participation;
            participation = participation1;
        else
            participation = [];
        end
        if Npeak2 > 0 && ~FAST_MODE
            [NetworkBurstPlayers.Secondary, a, b, c] = findChannel_in_Peak(BurstDetected, n_channels, t_step, Npeak2, NetworkBurstPlayers.Secondary, range, option_test, option);
            participation2 = NetworkBurstPlayers.Secondary.participation;
            participation = participation + participation2;
        end
        % run the final search using the participating clusters (it might change the original "NetworkBurstPlayers" because of the updated "participation")
        if sum(participation)>0
            
            if opt_shorcut == 0     % Normal statistical assessment (not option #3)
                if Npeak2 == 0
                    NetworkBurstPlayers.Secondary.bin = 0;
                end
                NetworkBurstPlayers = PopActivityPeak_Integration(participation, spikeTimestamp, BurstDetected, range, bin, threshold, 1, opt_distance, 1, 1, option_plot);  % information about the new population peaks                        
            end

            NetworkBurstPlayers.ChannelsIntegrated(:,1) = find(participation>0);            % index of all participating neurons
            NetworkBurstPlayers.Channels(:,1) = find(participation1>0);                     % index of participating neurons in primary events
            
            Selected_spikeTimestamp = spikeTimestamp(NetworkBurstPlayers.ChannelsIntegrated(:,1));     % for plotting
            NonSelected_spikeTimestamp = spikeTimestamp(find(~participation));                         % for plotting    

            n_peak1 = NetworkBurstPlayers.nBurst;                               % total number of primary network burst events
            n_peak1_level1 = NetworkBurstPlayers.nBurst_level1;                 % total number of level1 primary peaks 

            % load the timeStamp for units in primary events
            timeStamp = cell(n_peak1, n_channels);
            if ~FAST_MODE
                % (revised October 2026) the window of each event is calculated once (before, once per cluster), and there is no "clear" inside
                % the loops (it took most of the time of these loops); same results
                ini_Nburst = nan(n_peak1,1);     fin_Nburst = nan(n_peak1,1);
                for k=1:n_peak1                                                 % across all network bursts
                    ini_Nburst(k) = NetworkBurstPlayers.IndivBurst(k).start;    % starting point of the event k, in sec
                    fin_Nburst(k) = NetworkBurstPlayers.IndivBurst(k).last;     % last point of the event k, in sec
                    time_added = (fin_Nburst(k)-ini_Nburst(k))*SAFE_K;          % extra time added: SAFE_K*duration before the event, and 2*SAFE_K*duration after it
                    ini_Nburst(k) = ini_Nburst(k) - time_added;                 % just in case, add more data points to the left
                    fin_Nburst(k) = fin_Nburst(k) + time_added*2;               % just in case, add more data points to the right (twice as many: intended)
                    if n_channels > 0
                        NetworkBurstPlayers.IndivBurst(k).window = fin_Nburst(k) - ini_Nburst(k); % time window for the spike train extraction (based on the core duration of the burst), in sec
                    end
                end
                for ch=1:n_channels                                             % across all clusters (participating and non-participating)
                    for k=1:n_peak1                                             % across all network bursts
                        switch option_test
                            case 1  % old method
                                aux1 = find(BurstDetected{ch}.spikes(:,1)>=ini_Nburst(k)/K_MSEC_SEC);   % index for all spikes of neuron ch, after starting point of analysis window
                                aux2 = find(BurstDetected{ch}.spikes(:,1)<fin_Nburst(k)/K_MSEC_SEC);    % index for all spikes of neuron ch, before last point of analysis window
                                aux = intersect(aux1, aux2);                                            % index for all spikes of neuron ch, within the event k
    
                            case 2  % binnary search
                                [~, aux, ~] = windowMembership(BurstDetected{ch}.spikes(:,1), ini_Nburst(k)/K_MSEC_SEC, fin_Nburst(k)/K_MSEC_SEC, 1, 1);
                        end

                        if numel(aux)>0
                            timeStamp{k,ch} = BurstDetected{ch}.spikes(aux,1);      % time stamp of all spikes from neuron ch within the event k
                        else
                            timeStamp{k,ch} = [];                                   % do not have spikes within the network burst
                        end                
                    end            
                end
            end            
            NetworkBurstPlayers.timeStamp = timeStamp;                          % time stamp of the spikes within the network bursts

            % load the timeStamp for units in secondary events
            if n_peak1_level1 > 0 && ~isnan(NetworkBurstPlayers.Secondary.nBurst) 
                    
                n_peak2 = NetworkBurstPlayers.Secondary.nBurst;             % total number of primary network burst events
                timeStamp2 = cell(n_peak2, n_channels);
                    
                if n_peak2 > 0
                    [NetworkBurstPlayers.Secondary, a, b, c] = findChannel_in_Peak(BurstDetected, n_channels, t_step, n_peak2, NetworkBurstPlayers.Secondary, range, option_test, option);
                    NetworkBurstPlayers.Secondary.Channels(:,1) = find(NetworkBurstPlayers.Secondary.participation>0);   % index of participating neurons in secondary events
                    % (revised October 2026) windows calculated once and no "clear" inside the loops, as for the primary events; same results
                    ini_Nburst = nan(n_peak2,1);     fin_Nburst = nan(n_peak2,1);
                    for k=1:n_peak2                                                             % across all secondary network bursts
                        ini_Nburst(k) = NetworkBurstPlayers.Secondary.IndivBurst(k).start;      % starting point of the event k, in sec
                        fin_Nburst(k) = NetworkBurstPlayers.Secondary.IndivBurst(k).last;       % last point of the event k, in sec
                        time_added = (fin_Nburst(k)-ini_Nburst(k))*SAFE_K;                      % extra time added: SAFE_K*duration before the event, and 2*SAFE_K*duration after it
                        ini_Nburst(k) = ini_Nburst(k) - time_added;                             % just in case, add more data points to the left
                        fin_Nburst(k) = fin_Nburst(k) + time_added*2;                           % just in case, add more data points to the right (twice as many: intended)
                        if n_channels > 0
                            NetworkBurstPlayers.Secondary.IndivBurst(k).window = fin_Nburst(k) - ini_Nburst(k); % time window for the spike train extraction (based on the core duration of the burst), in sec
                        end
                    end
                    for ch=1:n_channels     
                        for k=1:n_peak2                                                         % across all secondary network bursts
                            switch option_test
                                case 1  % old method
                                    aux1 = find(BurstDetected{ch}.spikes(:,1)>=ini_Nburst(k)/K_MSEC_SEC);   % index for all spikes of neuron ch, after starting point of analysis window
                                    aux2 = find(BurstDetected{ch}.spikes(:,1)<fin_Nburst(k)/K_MSEC_SEC);    % index for all spikes of neuron ch, before last point of analysis window
                                    aux = intersect(aux1, aux2);                                            % index for all spikes of neuron ch, within the event k
    
                                case 2  % binnary search
                                    [~, aux, ~] = windowMembership(BurstDetected{ch}.spikes(:,1), ini_Nburst(k)/K_MSEC_SEC, fin_Nburst(k)/K_MSEC_SEC, 1, 1);
                            end

                            if numel(aux)>0
                                timeStamp2{k,ch} = BurstDetected{ch}.spikes(aux,1);             % time stamp of all spikes from neuron ch within the event k
                            else
                                timeStamp2{k,ch} = [];                                          % do not have spikes within the network burst
                            end
                        end
                    end
                end
                NetworkBurstPlayers.Secondary.timeStamp = timeStamp2;                            % time stamp of the spikes within the network bursts
            else
                NetworkBurstPlayers.Secondary.nBurst = 0;
            end
        else       % no participating units 
            NetworkBurstPlayers.ChannelsIntegrated = [];
            NetworkBurstPlayers.Channels = [];
            NetworkBurstPlayers.nBurst = 0;                                     
            NetworkBurstPlayers.timeStamp = [];
            NetworkBurstPlayers.Secondary.nBurst = 0;
            NonSelected_spikeTimestamp = spikeTimestamp(find(~participation));
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Only participating channels in primary peaks
    
    case 3

        NetworkBurstPlayers = PopulationActivityDetected;                   % load the information for the population peaks
        if Npeak1 > 0    
            [NetworkBurstPlayers, a, b, c] = findChannel_in_Peak(BurstDetected, n_channels, t_step, Npeak1, NetworkBurstPlayers, range, option_test, option);
            participation = NetworkBurstPlayers.participation;
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if  option_plot

   switch option
       case 1       % Initial screening: robust vs weak network bursts

           if Npeak1 > 0
                figure;       % population firing
                spikes = PopulationActivityDetected.spikesBin;
                time_bin = PopulationActivityDetected.time;
                bar(time_bin(:),spikes(:,1));    hold on
                if NetworkBurstPlayers.nBurst(1)>0
                    scatter(time_bin(peak1(:,1)),spikes(peak1(:,1),1)*1.2,'or');
                end    
                if numel(aux_weak)>0
                    for i=1:numel(aux_weak)
                        hold on;        scatter(time_bin(peak1(aux_weak(i))),spikes(peak1(aux_weak(i),1),1)*1.2,'*k');
                    end
                    text(time_bin(peak1(aux_weak(i))),spikes(peak1(aux_weak(i),1),1),"WEAK burst",'Color','k')
                end
                title(sprintf('Population firing - Network Burst Detection (%g potential / %g solid)',NetworkBurstPlayers.nBurst(1),NetworkBurstPlayers.nBurst(2)))
                xlabel('Time (sec)');       ylabel('Number of spikes')
           end
               
           if Npeak2>0 && ~FAST_MODE && numel(participating_channel2)>0
                for j=1:Npeak2
                    hold on; plot(NetworkBurstPlayers.Secondary.IndivBurst(j).peak, max(spikes)/3,'co', 'MarkerFaceColor','c'); hold on
                end
           end

            figure;       % raster plot of few example neurons
            if Npeak1>0
                n_examples = min(sum(participating_channel1), N_PLOT);
                for i=1:n_examples
                    idx = NetworkBurstPlayers.Channels(i,1);
                    n_unit = size(BurstDetected{idx}.spikes,1);
                    plot(BurstDetected{idx}.spikes(:,1), ones(n_unit,1)+SHIFT*i,'|',"MarkerEdgeColor",[0 0 0]);     hold on
                end
                title('Raster Plot - Robustly participating neurons and solid network bursts')
            else
                n_examples = min(N_PLOT,size(BurstDetected,1));
                for i=1:n_examples
                    n_unit = size(BurstDetected{i}.spikes,1);
                    plot(BurstDetected{i}.spikes(:,1), ones(n_unit,1)+SHIFT*i,'|',"MarkerEdgeColor",[0 0 0]);       hold on
                end
                title('Raster Plot - NO network bursting!!')
            end
            if n_examples==0
                i=1;
            end
            pico = 1+SHIFT*i;
            for k=1:numel(peak1)
                if  NetworkBurstPlayers.IndivBurst(k).solid               % only for solid/strong events
                    plot([NetworkBurstPlayers.IndivBurst(k).start/K_MSEC_SEC NetworkBurstPlayers.IndivBurst(k).start/K_MSEC_SEC], [0 pico], ':b');    hold on
                    plot([NetworkBurstPlayers.IndivBurst(k).last/K_MSEC_SEC NetworkBurstPlayers.IndivBurst(k).last/K_MSEC_SEC], [0 pico], ':r');      hold on
                end
            end
            
            aux = find(meaningful_activity==0);
            j=0;
            if numel(aux)>0
                for j=1:numel(aux)
                    idx = aux(j);
                    n_unit = size(BurstDetected{idx}.spikes,1);
                    plot(BurstDetected{idx}.spikes(:,1), ones(n_unit,1)+SHIFT*i+j,'|',"MarkerEdgeColor",'red');     hold on
                end
            end
            ylim([0.5 1.5+SHIFT*i+j])

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       case 2       % Statistical assessement
        
           figure;
           if sum(participation)>0
                subplot(3,1,[1,2]);
                for idx=1:size(Selected_spikeTimestamp,2)
                    n_unit = size(Selected_spikeTimestamp(idx).timeStamp,1);
                    plot(Selected_spikeTimestamp(idx).timeStamp(:,1)*K_MSEC_SEC, ones(n_unit,1)+idx,'|',"MarkerEdgeColor",[0 0 0]); hold on
                end
                ylabel('Neurons');  title('Selected units - Participating in network burst');  xlim([range(1) range(2)]); ylim([1 idx+2])
                
                subplot(3,1,3);
                if numel(NetworkBurstPlayers.timeStamp)>0
                    bar(NetworkBurstPlayers.time(:),NetworkBurstPlayers.spikesBin(:,1));
                    xlabel('Time (sec)');   ylabel('Number of spikes');  xlim([range(1) range(2)])
                end
           end

            figure;
            subplot(3,1,[1,2]);
            if size(NonSelected_spikeTimestamp,2)>0
                for idx=1:size(NonSelected_spikeTimestamp,2)
                    n_unit = size(NonSelected_spikeTimestamp(idx).timeStamp,1);
                    if n_unit>0
                        plot(NonSelected_spikeTimestamp(idx).timeStamp(:,1)*K_MSEC_SEC, ones(n_unit,1)+idx,'|',"MarkerEdgeColor",[0 0 0]);  hold on
                    end
                end
                ylabel('Neurons');  title('Non-participating units');  xlim([range(1) range(2)]);   ylim([1 idx+2]) 
            end
            
           if sum(participation)>0
                subplot(3,1,3);                
                if numel(NetworkBurstPlayers.timeStamp)>0
                    bar(NetworkBurstPlayers.time(:),NetworkBurstPlayers.spikesBin(:,1));
                    xlabel('Time (sec)');   ylabel('Number of spikes');  xlim([range(1) range(2)])
                end
           end
   end
end
end