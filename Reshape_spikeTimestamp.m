 function Reshape_spikeTimestamp = Reshape_spikeTimestamp(spikeTimestamp, window)

% function to reshape the structure "spikeTimestamp" for the specified temporal window (in sec)
% (revised October 2026: logical mask instead of "find" + "intersect"; same result)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   structure array spikeTimestamp (its length is the number of recording channels), with the following fields
% 	.recLength:	length of the recording, in seconds.
% 	.nUnits:	total number of identified units in the channel
% 	.timeStamp:	time stamp for all units (nUnits) in the channel, in msec
%   .window:    [start end] time window for the unit extraction, in sec
% 	.position:	position of the electrode
% 	.cluster:	cluster id
%   .t_step:    temporal resolution, in msec
%   .positionXY: coordinates (X,Y) of the peak amplitude signal for the cluster 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                        % factor to rescale from msec to sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_cluster = size(spikeTimestamp,2);                                         % total number of clusters in the dataset
Reshape_spikeTimestamp = spikeTimestamp;                                    % initialize the output variable

for k=1:n_cluster                                                           % across all clusters
    Reshape_spikeTimestamp(k).window = window;                              % redefine the window
    Reshape_spikeTimestamp(k).recLength = window(2) - window(1);            % redefine the lenght of the window
    Reshape_spikeTimestamp(k).timeStamp = [];                               % initialize the spike time field    

    train = spikeTimestamp(k).timeStamp;                                    % original spike time data
    aux = train>window(1)/K_MSEC_SEC & train<window(2)/K_MSEC_SEC;          % all data points within the window (logical mask; same result as "intersect" of two "find", but faster)
    Reshape_spikeTimestamp(k).nUnits = sum(aux);                            % total number of events within the window
    if any(aux)
        Reshape_spikeTimestamp(k).timeStamp = train(aux);                   % spike time of all events within the window
    end
end

end