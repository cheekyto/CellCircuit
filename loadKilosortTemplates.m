function [templates, unitType, samplingRate] = loadKilosortTemplates(folder, spikeTimestamp)

% function to read the Kilosort template of each unit of "spikeTimestamp" on its peak electrode, to measure its waveform with
% "WaveformFeatures" (new, November 2026)
% The templates are the ones written next to the spike times after the spike sorting: "waveforms.mat" (cell array "waveforms", one
% template per row of "metadata.csv"), or, if that file is missing, the column "waveform_template" of "metadata.csv" (same values)
% Kilosort units only: the units of a 3Brain dataset (or of an integrated Kilosort + 3Brain dataset) have no Kilosort template

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   folder:             folder of the recording, with "metadata.csv" and "waveforms.mat"
%   spikeTimestamp:     structure array of the units (e.g., from the "_QCmanual.mat" file); the field "cluster" is the Kilosort cluster
%                       id (column "new_cluster_id" of "metadata.csv", see "ShapeSpikeTrainData"), and "t_step" the temporal resolution, in msec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   templates:          matrix with the template of each unit, one row per unit of "spikeTimestamp" (same order)
%   unitType:           type of each unit given by bombcell (=0, noise; =1, good; =2, MUA; =3, non-somatic)
%   samplingRate:       sampling rate of the recording, in Hz (from spikeTimestamp.t_step; NaN without units)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                                                        % factor to rescale from msec to sec

file_meta = fullfile(folder, 'metadata.csv');
file_wave = fullfile(folder, 'waveforms.mat');
if ~isfile(file_meta)
    error('loadKilosortTemplates: the file %s does not exist', file_meta);
end
opts = detectImportOptions(file_meta, 'FileType', 'text', 'Delimiter', ',');   % (the spaces inside the text of the templates must not be
                                                                            % taken as delimiters)
if ~all(ismember({'unitType', 'new_cluster_id'}, opts.VariableNames))
    error('loadKilosortTemplates: %s has no columns "unitType" and "new_cluster_id"', file_meta);
end
use_mat = isfile(file_wave);
if use_mat
    opts.SelectedVariableNames = {'unitType', 'new_cluster_id'};            % the text of the templates is not needed
end
metadata = readtable(file_meta, opts);

if use_mat
    S = load(file_wave, 'waveforms');
    if ~isfield(S, 'waveforms')
        error('loadKilosortTemplates: %s has no variable "waveforms"', file_wave);
    end
    rows = S.waveforms(:);
    if numel(rows) ~= height(metadata)
        error('loadKilosortTemplates: %s has %g templates, but metadata.csv has %g rows', file_wave, numel(rows), height(metadata));
    end
else                                                                        % the templates written as text in metadata.csv ("[ 0. 0. ... ]")
    text = string(metadata.waveform_template);
    rows = cell(numel(text),1);
    for k=1:numel(text)
        rows{k} = sscanf(char(regexprep(text(k), '[\[\]]', ' ')), '%f')';
    end
end

clusters = [];
if ~isempty(spikeTimestamp)
    clusters = double([spikeTimestamp.cluster]);
end
[found, row] = ismember(clusters, double(metadata.new_cluster_id));
if ~all(found)
    error('loadKilosortTemplates: %g units of spikeTimestamp are not in metadata.csv (e.g., cluster %g)', sum(~found), clusters(find(~found,1)));
end

n_samples = numel(rows{1});
templates = zeros(numel(clusters), n_samples);
for i=1:numel(clusters)
    w = double(rows{row(i)}(:))';
    if numel(w) ~= n_samples
        error('loadKilosortTemplates: the templates do not have the same number of samples');
    end
    templates(i,:) = w;
end
unitType = double(metadata.unitType(row));
unitType = unitType(:);

if isempty(spikeTimestamp)
    samplingRate = NaN;
else
    samplingRate = 1 / (spikeTimestamp(1).t_step * K_MSEC_SEC);
end

end
