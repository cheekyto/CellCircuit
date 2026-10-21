function ShapeSpikeTrainData(inputs)

% function to shape the spike train data coming out from spike sorting (performed by kilosort)
% (revised October 2026: the spikes are grouped by cluster only once -faster-; the spikes exactly on the boundary between two conditions are no longer lost)
% (optional) inputs:    structure to overwrite any variable of the sections INPUTS and PARAMETERS (e.g., inputs.folders, inputs.files_output); see "override_inputs"
%                       when it is not provided, the values defined in the code are used

if ~parameter_set('active'), clc; end; close all;      % (revised October 2026) no "clc" while "run_pipeline" runs (it would erase its progress messages)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system = 3;                 % =1, Ubuntu; =2, Yang's PC; =3, Yang's MacBook; =0, do not modify the MATLAB path (e.g., test harness)
if nargin < 1, inputs = struct(); end                                       % no inputs provided -> use the values defined in the code
if isfield(inputs, 'system'), system = inputs.system; end                   % the system can also be provided as an input
switch system
    case 0      % the MATLAB path is set up by the caller
        delimiter = filesep;
    case 1      % Anlysis PCP - Ubuntu
        addpath(genpath('/home/cortex/Documents/Yang'));
        delimiter = '/';
    case 2      % Yang's PC - GitHub
        addpath(genpath('C:\Users\Yang\OneDrive - Harvard University\Documents\GitHub\MEA-analysis\'));
        delimiter = '\\';
    case 3      % Yang's MacBook
        addpath(genpath('/Users/sungminyang/Documents/GitHub/CellCircuit'));
        delimiter = '/';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folders = [                                                                 % Folder where to find the original data (spike_times.npy and spike_clusters.npy)
    "/path/to/recording_folder"];

files_output = [                                                            % Name of the output file where the variable "spikeTimestamp" will be saved in                               
    "recording_name"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                        % factor to rescale from msec to sec

unitType_good = [1 2 3];                    % classification of the unit type, performed by "bombcell", that will be considered for downstream analysis
                                            % =0, noise; =1, good; =2, MUA; =3, non-somatic

apply_parameters('ShapeSpikeTrainData');                                    % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")
override_inputs(inputs);                    % overwrite the INPUTS / PARAMETERS above with the ones provided as input (if any)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   structure array spikeTimestamp (its length is equal to the number of sorted units), with the following fields
% 	.recLength:	length of the recording, in seconds.
% 	.nUnits:	total number of identified events or spikes for the unit
% 	.timeStamp:	time stamp of all events (nUnits) for the unit, in msec
%   .window:    [start end] time window for the unit extraction, in sec
% 	.position:	position of the electrode assigned to the unit
% 	.cluster:	cluster id (unit #; can be used to trace back the waveform data)
%   .t_step:    temporal resolution of the recording, in msec

%   rec_conditions:     list of the filename for the different conditions that were concatenated for the spike sorting
%   Indiv_Cond:         cell with individual "spikeTimestamp" structure arrays, one per recording condition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check the inputs (folder)

n_recordings = numel(files_output);                                        % # recording files
for k=1:n_recordings
    if ~isfolder(folders(k))
        error("The folder %s (index: %g) can not be found!! ", folders(k), k)
    end
end

for k=1:n_recordings                                                        % across all recordings/folders
    folder = folders(k)                                                     % folder with the raw data
    file_output = files_output(k);                                          % name for the output file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read the original data

file_time = sprintf('%s%sspike_times.npy', folder, delimiter);
time_spike = double(readNPY(file_time));                                    % single array with the time stamp for all units from all channels        

file_cluster = sprintf('%s%sspike_clusters.npy', folder, delimiter);
cluster = readNPY(file_cluster);                                            % codes for the cluster id

metadata = readtable(sprintf('%s%smetadata.csv', folder, delimiter));
unitType = metadata.unitType;                                               % type of unit; classification made by "bombcell"
cluster_id2 = metadata.new_cluster_id;
peak_channel = metadata.peak_channel;                                       % channel with the peak amplitude signal for each cluster
peak_chan_x = metadata.peak_chan_x;                                         % position x of the channel with the peak amplitude signal for each cluster
peak_chan_y = metadata.peak_chan_y;                                         % position y of the channel with the peak amplitude signal for each cluster

recording_data = readtable(sprintf('%s%srecording_data.csv', folder, delimiter));
sampling_f = recording_data.sampling_freq;                                  % sampling frequency, in Hz
rec_conditions = recording_data.recording_name;                             % list of the filename for the different conditions that were concatenated for the spike sorting
sample_start = recording_data.sample_start;                                 % start sample for each recording condition
sample_end = recording_data.sample_end;                                     % end sample for each recording condition
rec_samples = [sample_start sample_end];
n_conditions = size(rec_conditions,1);                                      % total number of recording conditions

fprintf('The number of recording conditions for %s is %g (at %g Hz) \n', file_output, n_conditions, median(sampling_f));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify the right type of unit

n_goodType = numel(unitType_good);                                          % total number of "good" unit types 
cluster_id = [];                                                            % array that will contain the id of all "good" clusters
for k=1:n_goodType                                                          % across all "good" unit types
    good = unitType==unitType_good(k);                                      % row index for clusters of type k
    cluster_id = unique(vertcat(cluster_id, cluster_id2(good)));            % id of all good clusters 
    clear good
end
n_cluster = numel(cluster_id);                                              % total number of good clusters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Group the spikes by cluster (revised October 2026)
% (before, all the spikes were searched once per cluster and per condition, which was the slowest part of the function)

[~, cluster_pos] = ismember(double(cluster), cluster_id);                   % position of the cluster of each spike within "cluster_id" (=0 for the clusters that are not "good")
is_good = find(cluster_pos > 0);                                            % index of all the spikes from "good" clusters
spike_index = accumarray(cluster_pos(is_good), is_good, [n_cluster 1], @(x) {sort(x)});    % spike_index{i}: index of all events/spikes in cluster i (in ascending order, as "find" returned them; empty if no spikes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Deconstruct the data

Indiv_Cond = cell(n_conditions,1);
for cond=1:n_conditions                                                     % across all recording conditions
    s_start = rec_samples(cond,1);                                          % start sample for condition cond          
    s_end = rec_samples(cond,2);                                            % end sample for condition cond
    sampling_rate = sampling_f(cond);                                       % sampling frequency, in Hz
    t_step = 1/sampling_rate/K_MSEC_SEC;                                    % temporal resolution, in msec
    
    spikeTimestamp = struct('timeStamp',cell(1,n_cluster),'nUnits',cell(1,n_cluster),'recLength',cell(1,n_cluster),'window',cell(1,n_cluster),'cluster',cell(1,n_cluster),'position',cell(1,n_cluster),'t_step',cell(1,n_cluster),'positionXY',cell(1,n_cluster));
    
    for i=1:n_cluster                                                       % across all different clusters
        cluster_idx = cluster_id(i);                                        % label of the cluster i
        aux_idx = spike_index{i};                                           % index of all events/spikes in cluster i (to read the raw data file)
        n_original = numel(aux_idx);                                        % number of total events/spikes in cluster i

        if n_original>0                                                     % there is at least one event/spike from cluser i
            spikes = time_spike(aux_idx);                                   % sample stamp of all events/spikes in cluster i
            idx_t = find(spikes>=s_start & spikes<s_end);                   % the units in the range [s_start, s_end)
                                                                            % (revised October 2026: before, the range was (s_start, s_end) and a spike exactly at the
                                                                            %  boundary between two conditions was lost from both of them)
            if numel(idx_t)>0                                               % there is at least one event/spike from cluster i within the condition cond
                spikeTimestamp(i).timeStamp = spikes(idx_t)*t_step;         % only the units in the window of interest, in msec
                spikeTimestamp(i).nUnits = numel(idx_t);                    % number of units in the window of interest
            else                                                            % if there is NO units the window of interest
                spikeTimestamp(i).timeStamp = [];
                spikeTimestamp(i).nUnits = 0;
            end
        else                    
            spikeTimestamp(i).timeStamp = [];                               % no events/spikes
            spikeTimestamp(i).nUnits = 0;
        end
    
        spikeTimestamp(i).recLength = (s_end-s_start)*t_step*K_MSEC_SEC;        % length of the window, in sec 
        spikeTimestamp(i).window = [s_start s_end]*t_step*K_MSEC_SEC;           % initial and last point of the window, in sec
        spikeTimestamp(i).cluster = cluster_idx;                                % label of the cluster i
        spikeTimestamp(i).position = peak_channel(cluster_id2==cluster_idx);    % position of the peak amplitude signal for the cluster i
        spikeTimestamp(i).t_step = t_step;                                      % temporal resolution, in msec
        spikeTimestamp(i).positionXY = [peak_chan_x(cluster_id2==cluster_idx) peak_chan_y(cluster_id2==cluster_idx)];    % coordinates (X,Y) of the peak amplitude signal for the cluster i
    end
    Indiv_Cond{cond} = spikeTimestamp;                                      % load the processed data for the condition cond
    clear spikeTimestamp
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Concatenate the spike train from different conditions (to finally generate a single spikeTimestamp structure array)

spikeTimestamp = Indiv_Cond{1};                                             % initialize with the first condition
if n_conditions>1
    for cond=2:n_conditions                                                 % across all recording conditions, except the first one
        for i=1:n_cluster                                                   % across all different clusters
            spikeTimestamp(i).window(2) = Indiv_Cond{cond}(i).window(2);                                        % extend the end of the window    
            spikeTimestamp(i).recLength = spikeTimestamp(i).recLength + Indiv_Cond{cond}(i).recLength;          % add the duration of the condition cond
            spikeTimestamp(i).nUnits = spikeTimestamp(i).nUnits + Indiv_Cond{cond}(i).nUnits;                   % add the number of units of the condition cond
            spikeTimestamp(i).timeStamp = vertcat(spikeTimestamp(i).timeStamp, Indiv_Cond{cond}(i).timeStamp);  % append the spike time recorded during condition cond
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the output (original data)

file_save = sprintf('%s%s%s_raw.mat', folder,delimiter,file_output);
save(file_save, 'spikeTimestamp','rec_conditions','Indiv_Cond');

fprintf('Total number of clusters identified by the spike sorting algorithm: %g\n\n', n_cluster);
    
end
end
