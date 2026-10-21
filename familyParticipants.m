function familyParticipants = familyParticipants(NetworkBurstDetected, BurstDetected, netburst_member, units_overall)

% function to identify the units that participate in a family of primary network bursts (revised October 2026)
% A unit participates in the family if it is a participating unit of the primary network bursts (it passed the statistical assessment
% of "findChannel_in_Peak" across all the primary events) AND it fires in the events of the family with the same participation rule
% used in that assessment: in more than MIN_PARTICIPATION of the family's events and in at least MIN_PARTICIPATION_ABS of them (or in
% its only event, if the family has a single event). The same fixed rule is used for every family and every recording (no statistical
% test per family, whose power would depend on the number of events of the family, and no fall back)
% It is used only to calculate the features of each family of primary events ("NetworkBurstFeatures", in "Analyze_spikeTrain")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   NetworkBurstDetected:   output of "NetworkBurstDetection" (.IndivBurst(j).start / .last, in sec)
%   BurstDetected{ch}:      output of "BurstDetection" for the cluster ch (.spikes(:,1): time stamp of the spikes, in msec)
%   netburst_member(j):     =1 (true), if the primary event j belongs to the family
%   units_overall:          index of the participating units of all the primary events (NetworkBurstDetected.Channels)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   familyParticipants(i,1):    index of the participating units of the family (a subset of units_overall, same order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (the same values as in "findChannel_in_Peak": if they are changed there, change them here too)

K_ONSET = 0.20;                                 % extra time added on each side of the event, in terms of its duration (as in "findChannel_in_Peak")
MIN_PARTICIPATION = 0.20;                       % minimum proportion of the family's events in which a unit has to participate (as in "findChannel_in_Peak")
MIN_PARTICIPATION_ABS = 2.0;                    % minimum number of the family's events in which a unit has to participate (as in "findChannel_in_Peak")
                                                % (if the family has fewer events, all of them)
apply_parameters('familyParticipants');                                     % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")
K_MSEC_SEC = 1/1000;                            % factor to rescale from msec to sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

units_overall = units_overall(:);
events = find(netburst_member);                                             % primary events of the family
n_events = numel(events);
if n_events == 0 || isempty(units_overall)
    familyParticipants = zeros(0,1);
    return
end

% windows of the events, as in the statistical assessment of "findChannel_in_Peak"
window = nan(n_events, 2);
for j=1:n_events
    ini0 = NetworkBurstDetected.IndivBurst(events(j)).start;
    fin0 = NetworkBurstDetected.IndivBurst(events(j)).last;
    window(j,1) = (ini0 - (fin0-ini0)*K_ONSET)/K_MSEC_SEC;                  % starting point of the event j, in msec
    window(j,2) = (fin0 + (fin0-ini0)*K_ONSET)/K_MSEC_SEC;                  % last point of the event j, in msec
end

n_active = zeros(numel(units_overall), 1);                                  % number of the family's events in which each unit fires
for i=1:numel(units_overall)
    spikes = BurstDetected{units_overall(i)}.spikes(:,1);
    n_active(i) = sum(windowCounts(spikes, window(:,1), window(:,2)) > 0);  % (same counting as "findChannel_in_Peak")
end

participate = n_active/n_events > MIN_PARTICIPATION & n_active >= min(MIN_PARTICIPATION_ABS, n_events);
familyParticipants = units_overall(participate);

end
