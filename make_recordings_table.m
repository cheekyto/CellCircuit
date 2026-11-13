function T = make_recordings_table(main_folder, output_file, varargin)

% function to build the recordings table of "run_pipeline" from a folder (new, November 2026): one row per recording found in the folder and
% in all its subfolders (at any depth), with the columns
%   folder:     path of the folder of the recording
%   file_name:  name of the recording, taken from the name of its file "<name>_QCmanual.mat" (the file name without "_QCmanual.mat")
%   picture:    (November 2026) photo of the chip found in the folder of the recording ("recording_picture"; empty when there is none, or
%               several: the folders with several photos are listed, to choose one by hand)
%   picture_organoid:   empty, to be filled by hand with the number of the organoid in the photo (needed to run "SelectElectrodes")
% The datasets written by "SelectElectrodes" ("<name>_<organoid_label>_QCmanual.mat", and the later files of the cropped dataset) are not
% recordings of their own and are left out ("cropped_datasets"): the cropping is run, or run again, from the row of the recording

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   main_folder:        folder to search (it and all its subfolders)
%   output_file:        (optional) Excel file to write (default: "recordings.xlsx" in main_folder); it is replaced if it exists
%   options (name-value pairs):
%       'suffix'            end of the name of the files that identify a recording (default: "_QCmanual.mat"); the rest of the file name is
%                           the name of the recording
%       'pipeline_names'    =true (default): the files that the pipeline writes later for the same recording, "<name>_curated_QCmanual.mat"
%                           ("IntegrateDataSets") and "<name>_integrated_QCmanual.mat", also give "<name>", the common name of the files of the
%                           recording that "run_pipeline" expects in the column file_name (one row per recording, even with several of
%                           these files in its folder); =false: the name is the file name without the suffix, as it is
%       'group'             (optional) value of a column "group" for all the rows ("run_pipeline" needs a group for each recording)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   T:      the table, sorted by folder and name (also written to output_file)

% Example:  make_recordings_table("/path/to/data")
%           make_recordings_table("/path/to/data", "/path/to/data/recordings.xlsx", 'group', "wt")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'suffix', "_QCmanual.mat");
addParameter(p, 'pipeline_names', true);
addParameter(p, 'group', "");
parse(p, varargin{:});
opt = p.Results;

main_folder = char(main_folder);
if ~isfolder(main_folder)
    error('The folder %s does not exist', main_folder);
end
if nargin < 2 || isempty(output_file) || strlength(string(output_file)) == 0
    output_file = fullfile(main_folder, 'recordings.xlsx');
end
suffix = char(opt.suffix);

% files "<name><suffix>" in the folder and all its subfolders (without the hidden folders, and the "._" files that macOS writes on
% external drives)
files = dir(fullfile(main_folder, '**', ['*' suffix]));
files = files(~[files.isdir]);
folder = strings(0,1);  file_name = strings(0,1);
for k = 1:numel(files)
    below = extractAfter(string(files(k).folder), strlength(string(main_folder)));    % part of the path below main_folder
    if startsWith(files(k).name, '._') || contains(below, [filesep '.'])
        continue
    end
    name = string(files(k).name);
    name = extractBefore(name, strlength(name) - strlength(suffix) + 1);   % file name without the suffix
    if opt.pipeline_names
        name = regexprep(name, '_(curated|integrated)$', '');              % files of the same recording written later by the pipeline
    end
    if strlength(name) == 0
        continue
    end
    folder(end+1,1) = string(files(k).folder); %#ok<AGROW>
    file_name(end+1,1) = name; %#ok<AGROW>
end
T = table(folder, file_name);
T = unique(T, 'rows');                                                      % one row per recording, sorted by folder and name
cropped = cropped_datasets(T.folder, T.file_name);                          % (November 2026) datasets written by "SelectElectrodes": the
n_cropped = nnz(cropped);                                                   % cropping is run (again) from the row of the recording
T(cropped, :) = [];
if strlength(string(opt.group)) > 0
    T.group = repmat(string(opt.group), height(T), 1);
end
% photo of the chip of each recording (November 2026)
T.picture = strings(height(T), 1);
T.picture_organoid = strings(height(T), 1);
several = strings(0, 1);
for k = 1:height(T)
    [T.picture(k), candidates] = recording_picture(T.folder(k), T.file_name(k));
    if strlength(T.picture(k)) == 0 && numel(candidates) > 1
        several(end+1,1) = T.folder(k) + ": " + strjoin(candidates, ", "); %#ok<AGROW>
    end
end

% folders right inside main_folder without any recording (e.g., a recording not curated yet)
sub = dir(main_folder);
sub = sub([sub.isdir] & ~startsWith({sub.name}, '.'));
empty = strings(0,1);
for k = 1:numel(sub)
    path_k = string(fullfile(main_folder, sub(k).name));
    if ~any(T.folder == path_k | startsWith(T.folder, path_k + filesep))
        empty(end+1,1) = path_k; %#ok<AGROW>
    end
end

writetable(T, output_file, 'WriteMode', 'replacefile');
fprintf('%d recordings found in %s; table written to %s\n', height(T), main_folder, output_file);
if n_cropped > 0
    fprintf(['%d cropped datasets (written by SelectElectrodes, "<recording>_<organoid_label>") left out: the cropping is run, or run ' ...
        'again, from the row of their recording (one row per organoid, with its organoid_label)\n'], n_cropped);
end
if ~isempty(empty)
    fprintf('%d folders without a file "*%s":\n', numel(empty), suffix);
    fprintf('   %s\n', empty);
end
fprintf('%d of %d recordings with a photo of the chip (column picture; write the number of the organoid in picture_organoid)\n', ...
    nnz(strlength(T.picture) > 0), height(T));
if ~isempty(several)
    fprintf('%d folders with several photos (write the one to use in the column picture):\n', numel(several));
    fprintf('   %s\n', several);
end

end
