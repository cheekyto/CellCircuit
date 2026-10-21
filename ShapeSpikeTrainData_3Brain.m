function ShapeSpikeTrainData_3Brain(inputs)

% function to generate the spike train data coming out directly from the 3Brain system
% (optional) inputs:    structure to overwrite any variable of the sections INPUTS and PARAMETERS (e.g., inputs.folderPath, inputs.condition_order); see "override_inputs"
%                       when it is not provided, the values defined in the code are used
% (revised October 2026: the waveforms of each cluster are extracted at once -before, a loop across all the spikes of the recording for every 
%  cluster, which took hours-, and the spikes are grouped by cluster in one pass; same results)

if ~parameter_set('active'), clc; end; close all;      % (revised October 2026) no "clc" while "run_pipeline" runs (it would erase its progress messages)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system = 3;                 % =1, Ubuntu; =2, Yang's PC; =3, Yang's MacBook; =0, do not modify the MATLAB path (e.g., test harness)
if nargin < 1, inputs = struct(); end                                       % no inputs provided -> use the values defined in the code
if isfield(inputs, 'system'), system = inputs.system; end                   % the system can also be provided as an input
switch system
    case 1      % Anlysis PCP - Ubuntu
        addpath(genpath('/home/cortex/Documents/Yang'));
    case 2      % Yang's PC - GitHub
        addpath(genpath('C:\Users\Yang\OneDrive - Harvard University\Documents\GitHub\MEA-analysis\'));
    case 3      % Yang's MacBook
        addpath(genpath('/Users/sungminyang/Documents/GitHub/CellCircuit'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folderPath = '/path/to/recording_folder';             % Folder where to find the original data (xxx.bxr files)
% mainFiles = dir(fullfile(folderPath, '*Baseline.bxr'));                 % to identify ONLY the main/baseline file of each experiment/dataset (revised October 2026: below, after "override_inputs")

suffix = '_Baseline.bxr';

option_explore_conditions = 0;                                              % =1, to check the original order of the conditions (these files might have been saved in no chronological order)
condition_order = [2;1;3;4];                                                  % set the temporal order of the conditions

option_version = 2;                                                         % =1, BrainWave v5.x or older
                                                                            % =2, BrainWave v6.x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                        % factor to rescale from msec to sec
MATRIX_SIZE = 64;

apply_parameters('ShapeSpikeTrainData_3Brain');                             % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")
override_inputs(inputs);                    % overwrite the INPUTS / PARAMETERS above with the ones provided as input (if any)
mainFiles = dir(fullfile(folderPath, '*Baseline.bxr'));                     % to identify ONLY the main/baseline file of each experiment/dataset (after "override_inputs", since it depends on folderPath)

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
%% Check the experiments in the folder

mainFiles = mainFiles(~[mainFiles.isdir]);                                  % If you want to exclude directories (just in case):
mainNames = cell(numel(mainFiles), 1);
L = length(suffix);
if numel(mainFiles)>0
    for k = 1:numel(mainFiles)
        fname = mainFiles(k).name;                                          % Identify the main filename of the different experiments in the folder
        if length(fname) >= L && strcmp(fname(end-L+1:end), suffix)
            mainNames{k} = fname(1:end-L);                                  % everything before "filter.bxr"
        else
            mainNames{k} = '';                                              % just in case, shouldn't happen
        end
    end
else
    error('There is NO Files with that termination!!!')
end
mainNames = mainNames(~cellfun(@isempty, mainNames));
mainNames = unique(mainNames);

groups = struct('main_filename', {}, 'files', {});
for i = 1:numel(mainNames)
    main = mainNames{i};
    pattern = [main '*.bxr'];
    fileList = dir(fullfile(folderPath, pattern));                          % list of files that belong to the group i
    fileList = fileList(~[fileList.isdir]);                                 % just files
    if isempty(fileList)
        continue;
    end
    [~, sortIdx] = sort([fileList.datenum]);                                % ascending: oldest -> newest
    fileList = fileList(sortIdx);
    allNames = {fileList.name}.';              
    fileNamesSorted = cell(size(allNames));
    for k = 1:numel(allNames)
        [~, baseName, ~] = fileparts(allNames{k});                          % baseName = name without extension
        fileNamesSorted{k} = baseName;
    end

    n_conditions = numel(allNames);
    fprintf('\n Processing the experiment: %s (%g conditions)\n', main, n_conditions);

    Indiv_Cond = cell(n_conditions,3);
    rec_conditions_pilot = fileNamesSorted;
    if size(rec_conditions_pilot,1)==numel(condition_order)
       for j=1:numel(condition_order)
           fileNamesSorted{j} = rec_conditions_pilot{condition_order(j)};
       end
    else
        fprintf("\n The number of conditions selected for re-organization does not match the original list")
    end
    
    if option_explore_conditions
        rec_conditions_pilot
        fileNamesSorted
    else                                                                    % re-organize the order of the conditions
        rec_conditions = fileNamesSorted;
        waveformLengths = nan(n_conditions,1);  samplingRates = nan(n_conditions,1); n_clusters = nan(n_conditions,1); dacFactors = nan(n_conditions,1); offsetValues = nan(n_conditions,1);
        shift = 0;

        for k=1:n_conditions                                                    % across all recordings/folders
            fprintf('       %s \n', fileNamesSorted{k});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the original data

file_input = fullfile(folderPath, sprintf('%s.bxr', fileNamesSorted{k}));
wellID = 'Well_A1';
%exploreHDF5(file_input);  % replace with your actual filenam

switch option_version
    case 1      % BrainWave v5.x or older

        samplingRate = h5readatt(file_input, '/','SamplingRate');
        chIdxs = h5read(file_input, '/Well_A1/StoredChIdxs');
        spikeTOC = h5read(file_input, '/Well_A1/SpikeTOC');
        toc = h5read(file_input, '/TOC');

        info = h5info(file_input);
        spikeSortingPerformed = 0;
        for i=1:size(info.Groups.Datasets,1)
            if strcmp(info.Groups.Datasets(i).Name,'SpikeUnits')
                spikeSortingPerformed = 1;
            end
        end
 
        minDigitalValue = h5readatt(file_input, '/','MinDigitalValue');
        maxDigitalValue = h5readatt(file_input, '/','MaxDigitalValue');
        minAnalogValue = h5readatt(file_input, '/','MinAnalogValue');
        maxAnalogValue = h5readatt(file_input, '/','MaxAnalogValue');

        spikeData = h5read(file_input, '/Well_A1/SpikeTimes');
        spikeDataCh = h5read(file_input, '/Well_A1/SpikeChIdxs');
        waveformLength = int64(h5readatt(file_input,'/Well_A1/SpikeForms','Wavelength'));
        wave = h5read(file_input, '/Well_A1/SpikeForms');

    case 2      % BrainWave v6.x    

        rawJson = h5read(file_input, '/ExperimentInfo');
        experimentInfo = jsondecode(rawJson{1});
        if isfield(experimentInfo.SignalConverter.SampleToTimeConverter, 'FrameRateHertz')     % (revised October 2026: by name; before, it relied on the order of the fields)
            samplingRate = experimentInfo.SignalConverter.SampleToTimeConverter.FrameRateHertz;
        else                                                                % (BrainWave v6 writes it as "FrameRateHerz")
            samplingRate = experimentInfo.SignalConverter.SampleToTimeConverter.FrameRateHerz;
        end
        chIdxs = experimentInfo.CorePlateData.StoredPlateElIdxs;  % 0-based electrode indexes
        toc = h5read(file_input, '/TOC');                              
        spikeTOC = h5read(file_input, ['/' wellID '/SpikeTOC']);            % [N] int64
        % toc = toc';   % now [N x 2]: 

        fileInfo = h5info(file_input, ['/' wellID]);
        datasetNames = {fileInfo.Datasets.Name};
        spikeSortingPerformed = 0;
        spikeSortingPerformed  = any(strcmp(datasetNames, 'SpikeUnits'));
        
        dac = experimentInfo.SignalConverter.DigitalToAnalogConverter;
        % (revised October 2026) each value in its own variable; before, the analog limits (in microvolts) were saved as the digital ones and vice versa,
        % so "dacFactor" was the inverse (0.256 digital levels per microvolt instead of 3.907 microvolts per level) and "offsetValue" was 2047.5 instead
        % of -8000 microvolts (for the reference recording); the waveforms plotted by "Classify3BrainWaveform" and "Compare_Kilosort_vs_3Brain" were not in microvolts
        minDigitalValue = double(dac.MinDigitalValue);
        maxDigitalValue = double(dac.MaxDigitalValue);
        minAnalogValue = double(dac.MinAnalogValueMicroVolt);
        maxAnalogValue = double(dac.MaxAnalogValueMicroVolt);

        spikeData = h5read(file_input, ['/' wellID '/SpikeTimes']);
        spikeDataCh = h5read(file_input, ['/' wellID '/SpikePlateElIdxs']);
        waveformInfo   = h5info(file_input, ['/' wellID '/SpikeForms']);
        waveformLength = 0;
        for i = 1:length(waveformInfo.Attributes)
            if strcmp(waveformInfo.Attributes(i).Name, 'Wavelength')
                waveformLength = int64(waveformInfo.Attributes(i).Value);
            end
        end
        wave = h5read(file_input, ['/' wellID '/SpikeForms']);
end

if spikeSortingPerformed                                                    % (revised October 2026) this code is meant for raw data (spikes per electrode)
    warning('ShapeSpikeTrainData_3Brain:spikeSorted', ['%s contains units sorted by BrainWave ("SpikeUnits"); this code is meant for raw data without spike sorting. ' ...
        'The spike times are organized per unit, but "Classify3BrainWaveform" and "Compare_Kilosort_vs_3Brain" look up the waveforms by electrode ' ...
        '(cluster+1), so they would show the waveforms of other units'], file_input);
end

dacFactor = (maxAnalogValue - minAnalogValue) / (maxDigitalValue - minDigitalValue);
offsetValue = minAnalogValue - dacFactor * minDigitalValue;
samplingRates(k) = samplingRate;
dacFactors(k) = dacFactor; 
offsetValues(k) = offsetValue;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% From the given start position and duration (in frames), find the corresponding range of spike positions using the TOC

tocStartIdx = 1;
FS = double(toc(1,:));                                                      % Frame Start (included)
tocEndIdx = size(toc,2);
FE = double(toc(2,:));                                                      % Frame End (included)

spikeStartPosition = spikeTOC(tocStartIdx);
spikeStartPosition = max(spikeStartPosition,1);
% spikeEndPosition = spikeTOC(tocEndIdx);
spikeEndPosition = int64(numel(spikeData));                                 % (revised October 2026) the last spike of the file; before, "spikeTOC(tocEndIdx)", the position 
                                                                            % where the LAST block of the TOC starts, so the spikes of the last block (last ~1 sec) were 
                                                                            % not read (329 of 248,376 spikes in the reference recording 82_83)
expDuration = (max(FE) - min(FS)) / samplingRate;                           % total duration of the recording, in sec

if spikeSortingPerformed
   spikeChUnit = h5read(file_input, ['/' wellID '/SpikeUnits']);
   spikeDataChUnits = spikeChUnit(spikeStartPosition:spikeEndPosition);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Collect the required spike data and waveform data

spikeDataTimestamps = double(spikeData(spikeStartPosition:spikeEndPosition));
spikeDataChIdxs = spikeDataCh(spikeStartPosition:spikeEndPosition);
dataLength = spikeEndPosition - spikeStartPosition;                         % (number of spikes - 1)

waveformLengths(k) = waveformLength;
spikeDataWaveforms = wave(((spikeStartPosition-1)*waveformLength+1):spikeEndPosition*waveformLength);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the raster dataset with spike times and convert spike times from frame position to seconds for each channel
% Create the structure "spikeTimestamp"

if spikeSortingPerformed
    n_cluster = double(max(spikeDataChUnits));
else
    n_cluster = double(numel(chIdxs));
end
n_clusters(k) = n_cluster; 
spikeTimestamp = struct('idx_toc',cell(1,n_cluster),'timeStamp',cell(1,n_cluster),'nUnits',cell(1,n_cluster),'recLength',cell(1,n_cluster),'window',cell(1,n_cluster),'cluster',cell(1,n_cluster),'position',cell(1,n_cluster),'t_step',cell(1,n_cluster),'positionXY',cell(1,n_cluster));

% (revised October 2026) the spikes are grouped by cluster in one pass (before, one "find" across all the spikes for every cluster); same indexes,
% in ascending order: spike_group{i} = find(spikeDataChUnits==i) (spike-sorted data) or find(spikeDataChIdxs==chIdxs(i)) (electrodes)
if spikeSortingPerformed
    spike_label = double(spikeDataChUnits(:));                              % cluster of each spike (1..n_cluster; other values are ignored, as before)
else
    [~, spike_label] = ismember(spikeDataChIdxs(:), chIdxs(:));             % position of the electrode of each spike in "chIdxs" (0 if it is not there)
    if numel(unique(chIdxs)) < numel(chIdxs)                                % (an electrode listed twice would get its spikes twice; "ismember" gives the first one)
        error('ShapeSpikeTrainData_3Brain: repeated electrodes in the list of stored electrodes');
    end
end
[sorted_label, spike_order] = sort(spike_label);                            % stable sort: the spikes of each cluster keep their original (ascending) order
valid = sorted_label >= 1 & sorted_label <= n_cluster & sorted_label == fix(sorted_label);
spike_group = accumarray(sorted_label(valid), spike_order(valid), [n_cluster 1], @(x) {sort(x)}, {zeros(0,1)});

for i=1:n_cluster                                                           % across all different clusters
    if spikeSortingPerformed
        aux_idx = spike_group{i};                                           % index of all units in cluster i
        aux = spikeDataChIdxs(aux_idx);
        if numel(aux)>0
            cluster_id = double(mode(aux));
            pos_x = fix(cluster_id/MATRIX_SIZE)+1;
            pos_y = rem(cluster_id, MATRIX_SIZE)+1;
         else
            cluster_id = nan;
            pos_x = nan;
            pos_y = nan;
        end
    else
        aux_idx = spike_group{i};                                           % index of all units in channel i
        cluster_id = chIdxs(i);
        pos_x = fix((i-1)/MATRIX_SIZE)+1;
        pos_y = rem(i-1, MATRIX_SIZE)+1;
    end
    n_units = numel(aux_idx);                                               % number of total units in cluster i
 
    if n_units>0                                                            % there is at least one unit from cluser i 
       spikes = spikeDataTimestamps(aux_idx); 
       spikes = double(spikes)/samplingRate;         
       spikeTimestamp(i).timeStamp = spikes/K_MSEC_SEC;                     % spike time stamp for cluster i, in msec
       spikeTimestamp(i).nUnits = n_units;                                  % number of units in the window of interest
       spikeTimestamp(i).idx_toc = aux_idx;                                 % index of all units in channel/cluster i (to be used as index for "spikeDataChIdxs")
    else                    
        spikeTimestamp(i).timeStamp = [];                                   % no units
        spikeTimestamp(i).nUnits = 0;
        spikeTimestamp(i).idx_toc = [];
    end
   
    spikeTimestamp(i).recLength = expDuration;                              % length of the window, in sec 
    spikeTimestamp(i).window = [shift shift+expDuration];                   % initial and last point of the window, in sec
    spikeTimestamp(i).cluster = cluster_id;                                 % label of the cluster i (0-4095)
    spikeTimestamp(i).position = nan;                                       % position of the peak amplitude signal for the cluster i
    spikeTimestamp(i).t_step = 1/samplingRate/K_MSEC_SEC;                   % temporal resolution, in msec
    spikeTimestamp(i).positionXY = [pos_x-1; pos_y-1];                      % coordinates (X,Y) of the peak amplitude signal for the cluster i ([0, 63], to make it compatible with the Kilosort data)

    clear spikes n_units aux_idx aux cluster_id
end

Indiv_Cond{k,1} = spikeTimestamp;                                           % load the processed data for the condition k
shift = shift + expDuration;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the waveform dataset with individual traces for each electrode/unit

% (revised October 2026) all the waveforms of a cluster are extracted at once (before, a loop across all the spikes of the recording with "ismember",
% for every cluster); same rows, in the same order. "waveform2" (the waveform of each spike) includes all the spikes of the cluster (before, the
% last spike of the file, dataLength+1, was not included); "waveform1" keeps the spikes i = 1..dataLength (its samples extend into the next spike)
wl = double(waveformLength);                                                % length of each waveform, in samples
for idx=1:n_cluster
    aux_idx = spikeTimestamp(idx).idx_toc;                                  % index of all units in channel/cluster i (to be used as index for "spikeDataChIdxs")
    sel2 = double(aux_idx(:));                                              % all the spikes of the cluster
    sel1 = double(aux_idx(aux_idx <= dataLength));   sel1 = sel1(:);        % spikes with a complete "waveform1" window (as in the loop "for i=1:dataLength")
    waveform1 = zeros(0, wl);
    waveform2 = zeros(0, wl);
    if ~isempty(sel1)
        waveform1 = double(spikeDataWaveforms(sel1*wl + (0:wl-1)));         % samples i*wl : (i+1)*wl-1 of each spike i (one row per spike)
        waveform1 = reshape(waveform1, numel(sel1), wl);                    % (a single spike: make sure it is a row)
    end
    if ~isempty(sel2)
        waveform2 = double(spikeDataWaveforms((sel2-1)*wl + (1:wl)));       % samples (i-1)*wl+1 : i*wl of each spike i (one row per spike)
        waveform2 = reshape(waveform2, numel(sel2), wl);
    end
    waveformData1{idx} = waveform1;
    waveformData2{idx} = waveform2;
    clear waveform1 waveform2 sel1 sel2
end

Indiv_Cond{k,2} = waveformData2;                                           % load the waveform2 for the condition k
Indiv_Cond{k,3} = waveformData1;                                           % load the waveform1 data for the condition k
clear spikeTimestamp

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Concatenate the spike train from different conditions (to finally generate a single spikeTimestamp structure array)

spikeTimestamp = Indiv_Cond{1};                                             % initialize with the first condition
if n_conditions>1
    for cond=2:n_conditions                                                 % across all recording conditions, except the first one
        for i=1:n_cluster                                                   % across all different clusters
            spikeTimestamp(i).window(2) = Indiv_Cond{cond,1}(i).window(2);                                        % extend the end of the window    
            spikeTimestamp(i).recLength = spikeTimestamp(i).recLength + Indiv_Cond{cond,1}(i).recLength;          % add the duration of the condition cond
            spikeTimestamp(i).nUnits = spikeTimestamp(i).nUnits + Indiv_Cond{cond,1}(i).nUnits;                   % add the number of units of the condition cond
            spikeTimestamp(i).timeStamp = vertcat(spikeTimestamp(i).timeStamp, Indiv_Cond{cond,1}(i).timeStamp+Indiv_Cond{cond,1}(i).window(1)/K_MSEC_SEC);  % append the spike time recorded during condition cond
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the output (original data)

file_save = fullfile(folderPath, sprintf('%s_3Brain_raw.mat',rec_conditions{1}));
save(file_save, 'spikeTimestamp','rec_conditions','Indiv_Cond','waveformLengths','samplingRates','n_clusters','dacFactors','offsetValues');

    end             % end of the if condition to explore the conditions
end                 % end of the for loops across experiments

end

