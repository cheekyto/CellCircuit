function PopulationActivitySize = PopulationActivitySize(units, spikeTimestamp, t_start, t_fin, opt_factor)

% function to check if each population event has enough active units (and enough units with multiple spikes)
% (revised October 2026: the spikes within all the events are counted at once with "windowCounts"; same result)

MIN_CHAN =  6;                              % minimum number of channels involved in a single network burst (default: 5)
MIN_MULTIPLE = 5;                           % minimum number of channels with multiple spikes involved in a single network burst (default: 3)
MULTIPLE_SPIKE = 2;                         % minimum number of spikes to label as multiple (default: 2)

SAFE_K = 0.40;                              % when saving the spikes for each network burst event, the extra time is added on each side (default: 0.4)                          
apply_parameters('PopulationActivitySize');                                 % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spikeData = spikeTimestamp;
for k=size(spikeData,2):-1:1
    if ~units(k)
        spikeData(k)= [];
    end
end                                                                         % (revised October 2026) this "end" was missing ("else" + "end" only closed the "if"), so the whole 
                                                                            % function ran once per cluster (same final result, since the last pass had all clusters filtered,
                                                                            % but N times slower); the last "end" of the file closes the function now
n_events = numel(t_start);
PopulationActivitySize = true(n_events, 1); 
n_channels = size(spikeData,2);
score = zeros(n_events, n_channels);

% (revised October 2026) the windows are calculated once, and the spikes of each cluster within all the events are counted at once with
% "windowCounts" (same result as "windowMembership" event by event, which took most of the time of this function); no loop across events
time_added = (t_fin(:)-t_start(:))*SAFE_K;                                  % extra time added on both ends
ini_Nburst = t_start(:) - time_added;                                       % just in case, add more data points to the left
fin_Nburst = t_fin(:) + time_added;                                         % just in case, add more data points to the right
for ch=1:n_channels                                                         % across all clusters         
    score(:, ch) = windowCounts(spikeData(ch).timeStamp, ini_Nburst, fin_Nburst);  % number of spikes of cluster ch within each event
end     % end of the for loop across channels

n_active = sum(score > 0, 2);                                               % number of active units in each event
n_multiple = sum(score > MULTIPLE_SPIKE, 2);                                % number of multiple-spike units in each event
PopulationActivitySize(n_active < MIN_CHAN * opt_factor) = false;           % not enough number of active units
PopulationActivitySize(n_multiple < MIN_MULTIPLE * opt_factor) = false;     % not enough number of multiple-spike units
            
end
