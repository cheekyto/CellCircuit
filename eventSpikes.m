function [timeStamp, window] = eventSpikes(BurstDetected, n_channels, starts, lasts, safe_k, option_test)

% (new, November 2026) spikes of each cluster within each event, in a window from SAFE_K x duration before its start to 2 x SAFE_K x duration
% after its last time. It is the code of "NetworkBurstPlayers" (primary and secondary events), moved here so that the events added by
% option 5 of "NetworkBurstDetection" ("add_large_option4") get their spikes from exactly the same code: a change here applies to all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   BurstDetected{ch}.spikes(:,1):  spike times of the cluster ch, in msec
%   n_channels:         number of clusters
%   starts(k), lasts(k):    start and last time of the event k, in sec
%   safe_k:             extra time added on each side, as a proportion of the duration of the event (SAFE_K of "NetworkBurstPlayers")
%   option_test:        =1, search with "find"; =2, binary search ("windowMembership")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   timeStamp{k,ch}:    spike times of the cluster ch within the window of the event k, in msec ([] without spikes)
%   window(k):          length of the window of the event k, in sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;
n_events = numel(starts);
ini_Nburst = nan(n_events,1);     fin_Nburst = nan(n_events,1);
for k=1:n_events                                                            % across all events
    ini_Nburst(k) = starts(k);                                              % starting point of the event k, in sec
    fin_Nburst(k) = lasts(k);                                               % last point of the event k, in sec
    time_added = (fin_Nburst(k)-ini_Nburst(k))*safe_k;                      % extra time added: SAFE_K*duration before the event, and 2*SAFE_K*duration after it
    ini_Nburst(k) = ini_Nburst(k) - time_added;                             % just in case, add more data points to the left
    fin_Nburst(k) = fin_Nburst(k) + time_added*2;                           % just in case, add more data points to the right (twice as many: intended)
end
window = fin_Nburst - ini_Nburst;                                           % time window for the spike train extraction, in sec

timeStamp = cell(n_events, n_channels);
for ch=1:n_channels                                                         % across all clusters (participating and non-participating)
    for k=1:n_events                                                        % across all events
        switch option_test
            case 1  % old method
                aux1 = find(BurstDetected{ch}.spikes(:,1)>=ini_Nburst(k)/K_MSEC_SEC);   % index for all spikes of neuron ch, after starting point of analysis window
                aux2 = find(BurstDetected{ch}.spikes(:,1)<fin_Nburst(k)/K_MSEC_SEC);    % index for all spikes of neuron ch, before last point of analysis window
                aux = intersect(aux1, aux2);                                            % index for all spikes of neuron ch, within the event k

            case 2  % binary search
                [~, aux, ~] = windowMembership(BurstDetected{ch}.spikes(:,1), ini_Nburst(k)/K_MSEC_SEC, fin_Nburst(k)/K_MSEC_SEC, 1, 1);
        end

        if numel(aux)>0
            timeStamp{k,ch} = BurstDetected{ch}.spikes(aux,1);              % time stamp of all spikes from neuron ch within the event k
        else
            timeStamp{k,ch} = [];                                           % do not have spikes within the event
        end
    end
end

end
