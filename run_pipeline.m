function run_pipeline(recordings, varargin)

% MASTER FUNCTION of the pipeline (revised October 2026)
% It runs the stages of the pipeline for the recordings listed in a table, with all the parameters taken from a single parameter file
% ("pipeline_parameters.m"), and it saves the parameters used for each recording
%
% End to end (Kilosort data):   ShapeSpikeTrainData -> QualityControl_General -> QualityControl_HyperSynch -> QualityControl_Manual ->
%                               SelectElectrodes (only the recordings with cropping) -> IntegrateDataSets -> Analyze_spikeTrain,
%                               for each recording; then Compile_MEAproject, for each group
% Step by step:                 only the stages listed in 'stages' (or from the stage 'from'), including the 3Brain stage
%                               "ShapeSpikeTrainData_3Brain"
%
% Examples:
%   run_pipeline("recordings.csv")                                                      % end to end, with "pipeline_parameters.m"
%   run_pipeline("recordings.csv", 'stages', "Analyze_spikeTrain")                      % a single stage (its messages on screen)
%   run_pipeline("recordings.csv", 'from', "IntegrateDataSets")                         % from a stage to the end
%   run_pipeline("recordings.csv", 'parameters', "parameters_threshold5.m")             % a different parameter set
%   run_pipeline("recordings.csv", 'parameters', "saved")                               % the same set as the last run of each recording
%   run_pipeline("recordings.csv", 'changes', {'NetworkBurstDetection.method3.THRESHOLD', 5.0})   % a few values changed
%   run_pipeline("recordings.csv", 'rows', [2 3])                                       % only some recordings of the table

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   recordings:     path of the recordings table (.csv or .xlsx; see "recordings_template.csv"), one row per recording, with the columns
%       folder              folder of the recording (absolute, or relative to the folder of the table)
%       file_name           common name of the files of the recording (e.g., "recording_01"); for a 3Brain recording,
%                           the name of its files without the suffix (e.g., "XXX_3Brain")
%       group               group of the recording: the analyzed files of a group are collected in the folder <project>/<group>,
%                           and compiled together by "Compile_MEAproject"
%       (optional columns)
%       to_discard          clusters discarded by "QualityControl_Manual" (e.g., "3 17 25"; the indexes are defined after eliminating
%                           the silent clusters); empty: none
%       organoid_label      label of the organoid, with crop_ellipse ("x0 y0 radiusx radiusy", in pixels) or crop_excel (file of the
%                           3Brain electrode map): "SelectElectrodes" runs only for the recordings with one of them
%       file_name_3Brain    name of a 3Brain dataset to integrate with the Kilosort dataset ("IntegrateDataSets")
%       source              "Kilosort" (default) or "3Brain" (step by step only)
%       include             =0, skip the recording (default: 1)
%
%   options (name-value pairs):
%       'parameters'        parameter file (default: "pipeline_parameters.m" of this folder), or "saved" (the set saved for each
%                           recording by its last run)
%       'changes'           values changed for this run, on top of the parameter set: {'Function.PARAMETER', value, ...}
%       'stages'            stages to run (they always run in the order of the pipeline) (default: all, end to end)
%       'from'              first stage to run; the following ones also run
%       'rows'              recordings to run: rows of the table (1 = first recording) or file names (default: all, with include ~= 0)
%       'project'           folder for the group folders and the compiled files (default: the folder of the table)
%       'verbose'           =0, only the progress is printed (the messages of the stages go to the log file of each recording)
%                           =1, the messages of the stages are also printed (default: 1 for a single stage, 0 otherwise)
%       'stop_on_error'     =1, stop at the first error; =0, continue with the next recording (default: 0)
%       'extra_path'        folders to add to the MATLAB path (e.g., npy-matlab, needed by "ShapeSpikeTrainData")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   the files of each stage (as when the stage is run by hand), and
%   <folder>/<file_name>_parameters.mat:    parameters (set used by the last run), parameters_stage.(stage) (set used by the last run of
%                                           each stage), history (one entry per run: date, stages, parameter file, changes, version of
%                                           the code, status)
%   <folder>/<file_name>_pipeline_log.txt:  messages of the stages (appended at each run)
%   <project>/<group>/:                     copy of the "_Analyzed.mat" file of each recording of the group
%   <project>/<group>_Compiled.mat/.xlsx, <group>_Grouped&Compiled.xlsx, <group>_parameters.mat, <group>_pipeline_log.txt

% NOTE: "QualityControl_General" processes every "*_raw.mat" file of the folder, so each recording must be in its own folder
% NOTE: while run_pipeline runs, the parameters written in the code of each function are replaced by the parameter set (see
%       "apply_parameters"); a stage run by hand still uses the values in its code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STAGES = ["ShapeSpikeTrainData", "ShapeSpikeTrainData_3Brain", "QualityControl_General", "QualityControl_HyperSynch", "QualityControl_Manual", ...
          "SelectElectrodes", "IntegrateDataSets", "Analyze_spikeTrain", "Compile_MEAproject"];      % order of the pipeline
END_TO_END = STAGES(STAGES ~= "ShapeSpikeTrainData_3Brain");                % end to end (Kilosort)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options

opts = struct('parameters', "", 'changes', {{}}, 'stages', [], 'from', "", 'rows', [], 'project', "", 'verbose', [], ...
    'stop_on_error', 0, 'extra_path', {{}});
if mod(numel(varargin), 2) ~= 0
    error("The options must be name-value pairs (e.g., run_pipeline(""recordings.csv"", 'stages', ""Analyze_spikeTrain""))");
end
for k=1:2:numel(varargin)
    name = char(varargin{k});
    if ~isfield(opts, name)
        error("Unknown option '%s' (options: %s)", name, strjoin(fieldnames(opts), ', '));
    end
    opts.(name) = varargin{k+1};
end

code_dir = fileparts(mfilename('fullpath'));
extra_path = cellstr(string(opts.extra_path));
for k=1:numel(extra_path)
    if strlength(extra_path{k}) > 0
        addpath(extra_path{k});
    end
end

% stages to run
if ~isempty(opts.stages)
    stages = string(opts.stages);
    unknown = stages(~ismember(stages, STAGES));
    if ~isempty(unknown)
        error("Unknown stage '%s' (stages: %s)", unknown(1), strjoin(STAGES, ', '));
    end
    stages = STAGES(ismember(STAGES, stages));                              % in the order of the pipeline
elseif strlength(string(opts.from)) > 0
    idx = find(END_TO_END == string(opts.from));
    if isempty(idx)
        error("Unknown stage '%s' for 'from' (stages: %s)", string(opts.from), strjoin(END_TO_END, ', '));
    end
    stages = END_TO_END(idx:end);
else
    stages = END_TO_END;
end
end_to_end = isempty(opts.stages) && strlength(string(opts.from)) == 0;
verbose = opts.verbose;
if isempty(verbose)
    verbose = isscalar(stages);
end

% recordings
table_file = char(recordings);
if ~is_absolute(table_file)
    table_file = fullfile(pwd, table_file);
end
if ~isfile(table_file)
    error("The recordings table '%s' was not found", table_file);
end
table_dir = fileparts(table_file);
recs = read_recordings(table_file, table_dir);
recs = select_rows(recs, opts.rows);
if isempty(recs)
    error("No recording to run (check the column 'include' and the option 'rows')");
end
n_rec = numel(recs);
project = string(opts.project);
if strlength(project) == 0
    project = string(table_dir);
end

% parameter set
use_saved = strcmpi(string(opts.parameters), "saved");
if use_saved
    source_text = "the set saved for each recording (last run)";
    P_file = [];
else
    param_file = string(opts.parameters);
    if strlength(param_file) == 0
        param_file = string(fullfile(code_dir, 'pipeline_parameters.m'));
    end
    P_file = load_parameters(param_file, opts.changes);                     % (it also checks the changes)
    source_text = param_file;
end
changes_text = changes_description(opts.changes);
if any(stages == "ShapeSpikeTrainData") && any([recs.source] ~= "3Brain") && exist('readNPY', 'file') ~= 2
    error("The function readNPY (npy-matlab), needed by ShapeSpikeTrainData, is not on the MATLAB path: add it, or use the option 'extra_path'");
end
code_ver = code_version(code_dir);
cleanup = onCleanup(@() parameter_set('clear'));                            % the set is deactivated at the end (also after an error)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the stages, recording by recording

t_start = tic;
fprintf("\n=== run_pipeline (%s): %d recording(s); stages: %s\n", now_text('yyyy-MM-dd HH:mm'), n_rec, strjoin(stages, ' > '));
fprintf("    parameters: %s%s\n", source_text, changes_text);
if ~verbose
    fprintf("    the messages of the stages are saved in the log file of each recording (<file_name>_pipeline_log.txt)\n");
end

rec_stages = stages(stages ~= "Compile_MEAproject");
errors = strings(0,1);
ok_rec = true(n_rec,1);
for r=1:n_rec
    if isempty(rec_stages)
        break
    end
    rec = recs(r);
    fprintf("\n[%d/%d] %s (group: %s) - %s\n", r, n_rec, rec.file_name, rec.group, now_text('HH:mm'));
    if end_to_end && rec.source == "3Brain"
        fprintf("      skipped: 3Brain recordings run step by step (option 'stages')\n");
        ok_rec(r) = false;
        continue
    end
    file_params = fullfile(rec.folder, rec.file_name + "_parameters.mat");
    file_log = fullfile(rec.folder, rec.file_name + "_pipeline_log.txt");
    if ~isfolder(rec.folder)
        msg = sprintf("the folder %s does not exist", rec.folder);
        fprintf("      ERROR: %s\n", msg);
        errors(end+1) = sprintf("[%d/%d] %s: %s", r, n_rec, rec.file_name, msg);   %#ok<AGROW>
        ok_rec(r) = false;
        continue
    end

    % parameter set of the recording
    try
        if use_saved
            S = load_saved(file_params);
            P = load_parameters(S.parameters, opts.changes);
        else
            P = P_file;
        end
    catch ME
        fprintf("      ERROR: %s\n", ME.message);
        errors(end+1) = sprintf("[%d/%d] %s: %s", r, n_rec, rec.file_name, ME.message);   %#ok<AGROW>
        ok_rec(r) = false;
        continue
    end
    parameter_set('set', P);
    record = start_record(file_params, P, rec_stages, source_text, changes_text, opts.changes, code_ver);

    for st = rec_stages
        [inputs, skip] = stage_inputs(st, rec, end_to_end, verbose);
        if strlength(skip) > 0
            fprintf("      %-28s skipped (%s)\n", st, skip);
            continue
        end
        [ME, elapsed] = run_stage(st, inputs, verbose, file_log);
        if isempty(ME)
            fprintf("%s\n", done_text(elapsed, verbose, st));
            record.parameters_stage.(st) = P;
            record.history(end).stages_done(end+1) = st;
            if st == "Analyze_spikeTrain"                                   % collect the analyzed file in the folder of its group
                collect_analyzed(rec, project);
            end
        else
            fprintf("%s\n", error_text(ME, elapsed, verbose, st, file_log));
            record.history(end).status = "error in " + st + ": " + first_line(ME.message);
            errors(end+1) = sprintf("[%d/%d] %s, %s: %s", r, n_rec, rec.file_name, st, first_line(ME.message));   %#ok<AGROW>
            ok_rec(r) = false;
        end
        save_record(file_params, record);
        if ~isempty(ME)
            if opts.stop_on_error
                rethrow(ME);
            end
            break                                                           % the next stages depend on this one
        end
    end
    if ok_rec(r)
        record.history(end).status = "done";
        save_record(file_params, record);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compile the results, group by group

if any(stages == "Compile_MEAproject")
    groups = unique([recs.group], 'stable');
    for g = groups
        fprintf("\n[group %s] %d recording(s) in %s - %s\n", g, sum([recs.group] == g), fullfile(project, g), now_text('HH:mm'));
        if any(~ok_rec([recs.group] == g))
            fprintf("      WARNING: some recordings of this group failed; their previous results (if any) are still in the group folder\n");
        end
        file_params = fullfile(project, g + "_parameters.mat");
        file_log = fullfile(project, g + "_pipeline_log.txt");
        try
            if ~isfolder(fullfile(project, g)) || isempty(dir(fullfile(project, g, '*_Analyzed.mat')))
                error("there are no ""_Analyzed.mat"" files in %s", fullfile(project, g));
            end
            if use_saved
                S = load_saved(file_params);
                P = load_parameters(S.parameters, opts.changes);
            else
                P = P_file;
            end
        catch ME
            fprintf("      ERROR: %s\n", ME.message);
            errors(end+1) = sprintf("[group %s]: %s", g, ME.message);      %#ok<AGROW>
            continue
        end
        parameter_set('set', P);
        record = start_record(file_params, P, "Compile_MEAproject", source_text, changes_text, opts.changes, code_ver);
        inputs = struct('folderPath', char(project), 'group', char(g), 'system', 0);
        [ME, elapsed] = run_stage("Compile_MEAproject", inputs, verbose, file_log);
        if isempty(ME)
            fprintf("%s\n", done_text(elapsed, verbose, "Compile_MEAproject"));
            record.parameters_stage.Compile_MEAproject = P;
            record.history(end).stages_done(end+1) = "Compile_MEAproject";
            record.history(end).status = "done";
        else
            fprintf("%s\n", error_text(ME, elapsed, verbose, "Compile_MEAproject", file_log));
            record.history(end).status = "error in Compile_MEAproject: " + first_line(ME.message);
            errors(end+1) = sprintf("[group %s] Compile_MEAproject: %s", g, first_line(ME.message));   %#ok<AGROW>
        end
        save_record(file_params, record);
        if ~isempty(ME) && opts.stop_on_error
            rethrow(ME);
        end
    end
end

fprintf("\n=== run_pipeline finished in %s: %d recording(s), %d error(s)\n", time_text(toc(t_start)), n_rec, numel(errors));
for k=1:numel(errors)
    fprintf("    %s\n", errors(k));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%   LOCAL FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [inputs, skip] = stage_inputs(stage, rec, end_to_end, verbose)

% function to define the inputs of a stage for one recording (folders, file names, and the options that must not stop the run)

skip = "";
base = dataset_name(rec);
switch stage
    case "ShapeSpikeTrainData"
        if rec.source == "3Brain"
            skip = "3Brain recording: use ShapeSpikeTrainData_3Brain";
        end
        inputs = struct('folders', string(rec.folder), 'files_output', string(rec.file_name));
    case "ShapeSpikeTrainData_3Brain"
        if rec.source ~= "3Brain"
            skip = "Kilosort recording";
        end
        inputs = struct('folderPath', char(rec.folder));
    case "QualityControl_General"
        inputs = struct('folder', char(rec.folder), 'condition_3Brain', double(rec.source == "3Brain"), 'option_picture', 0, 'ellipse_cropped', nan);
        if ~verbose
            inputs.option_pause = 0;                                        % (it would wait for a key press, without any message on screen)
        end
    case "QualityControl_HyperSynch"
        inputs = struct('folders', string(rec.folder), 'files_data', string(rec.file_name), 'testing', 0);
    case "QualityControl_Manual"
        inputs = struct('folder', char(rec.folder), 'file_data', char(rec.file_name), 'to_discard', rec.to_discard);
        if end_to_end
            inputs.option_main = 3;                                         % no manual curation: the clusters in "to_discard" are discarded
            inputs.option_saveComplete = 1;
        end
    case "SelectElectrodes"
        if ~has_crop(rec)
            skip = "no cropping for this recording";
        end
        inputs = struct('folder', char(rec.folder), 'file_data', char(rec.file_name), 'organoid_label', char(rec.organoid_label));
        if strlength(rec.crop_excel) > 0
            inputs.option = 1;      inputs.file_excel = char(rec.crop_excel);
        else
            inputs.option = 2;      inputs.ellipse_cropped = rec.crop_ellipse;
        end
    case "IntegrateDataSets"
        if rec.source == "3Brain"                                           % a 3Brain dataset alone
            inputs = struct('folder', char(rec.folder), 'file_data_Kilosort', [], 'file_data_3Brain', char(base), 'include_Kilosort', 0, ...
                'include_3Brain', 1, 'file_output', char(base));
        else                                                                % Kilosort, alone or with a 3Brain dataset
            with_3Brain = strlength(rec.file_name_3Brain) > 0;
            data_3Brain = [];
            if with_3Brain
                data_3Brain = char(rec.file_name_3Brain);
            end
            inputs = struct('folder', char(rec.folder), 'file_data_Kilosort', char(base), 'file_data_3Brain', data_3Brain, 'include_Kilosort', 1, ...
                'include_3Brain', double(with_3Brain), 'file_output', char(base));
        end
    case "Analyze_spikeTrain"
        inputs = struct('folders', string(rec.folder), 'files_data', analyzed_name(rec));
        if ~verbose
            inputs.opt_checking = 0;                                        % (it would wait for a key press after each condition)
        end
    otherwise
        error("The stage '%s' is not run for each recording", stage);
end
inputs.system = 0;                                                          % the MATLAB path is already set (the stages do not modify it)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ME, elapsed] = run_stage(stage, inputs, verbose, file_log)

% function to run one stage: with verbose = 0, its messages (and warnings) are captured and appended to the log file

t0 = tic;
if verbose
    fprintf("      >> %s (started %s)\n", stage, now_text('HH:mm'));
    ME = call_stage(stage, inputs);
    log_text = "(messages printed on screen)" + newline;
else
    fprintf("      %-28s ", stage);
    log_text = string(evalc('ME = call_stage(stage, inputs);'));            % the messages of the stage are captured
end
elapsed = toc(t0);
close all;                                                                  % the figures of the stage (they were saved by the stage, if requested)

status = "done";
if ~isempty(ME)
    status = "ERROR";
    log_text = log_text + newline + string(getReport(ME, 'extended', 'hyperlinks', 'off')) + newline;
end
log_text = regexprep(log_text, '</?(a|strong|span)[^>]*>', '');            % hyperlinks of the warnings
fid = fopen(file_log, 'a');
if fid > 0
    fprintf(fid, '\n==== %s - %s - %s (%s) ====\n', stage, now_text('yyyy-MM-dd HH:mm'), status, time_text(elapsed));
    fprintf(fid, '%s', log_text);
    fclose(fid);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ME = call_stage(stage, inputs)

% function to run a stage and return its error (if any), so that the messages captured before the error are kept

ME = [];
try
    feval(stage, inputs);
catch ME
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function txt = done_text(elapsed, verbose, stage)

if verbose
    txt = sprintf("      << %s done (%s)", stage, time_text(elapsed));
else
    txt = sprintf("done (%s)", time_text(elapsed));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function txt = error_text(ME, elapsed, verbose, stage, file_log)

txt = sprintf("ERROR after %s: %s\n      (details in %s)", time_text(elapsed), first_line(ME.message), file_log);
if verbose
    txt = sprintf("      << %s ", stage) + txt;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function record = start_record(file_params, P, stages, source_text, changes_text, changes, code_ver)

% function to add a new entry to the history of the parameter file of a recording (or a group)

record = struct('parameters', P, 'parameters_stage', struct(), 'history', []);
if isfile(file_params)
    old = load(file_params);
    if isfield(old, 'parameters_stage'), record.parameters_stage = old.parameters_stage; end
    if isfield(old, 'history'), record.history = old.history; end
end
entry = struct('date', string(now_text('yyyy-MM-dd HH:mm:ss')), 'stages', string(stages), 'stages_done', strings(1,0), ...
    'parameter_file', string(source_text), 'changes', {changes}, 'changes_text', string(changes_text), 'code_version', string(code_ver), ...
    'status', "running");
if isempty(record.history)
    record.history = entry;
else
    record.history(end+1) = entry;
end
save_record(file_params, record);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_record(file_params, record)

parameters = record.parameters;
parameters_stage = record.parameters_stage;
history = record.history;
save(file_params, 'parameters', 'parameters_stage', 'history');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = load_saved(file_params)

if ~isfile(file_params)
    error("There is no saved parameter set (%s): run it first with a parameter file", file_params);
end
S = load(file_params, 'parameters');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function collect_analyzed(rec, project)

% function to copy the analyzed file of the recording into the folder of its group (input of "Compile_MEAproject")

group_dir = fullfile(project, rec.group);
if ~isfolder(group_dir)
    mkdir(group_dir);
end
copyfile(fullfile(rec.folder, analyzed_name(rec) + "_Analyzed.mat"), group_dir);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function base = dataset_name(rec)

% name of the dataset after "SelectElectrodes" (the cropped dataset has the label of the organoid)

base = rec.file_name;
if has_crop(rec)
    base = rec.file_name + "_" + rec.organoid_label;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function name = analyzed_name(rec)

% name of the input of "Analyze_spikeTrain" (output of "IntegrateDataSets", without "_QCmanual.mat")

if rec.source ~= "3Brain" && strlength(rec.file_name_3Brain) > 0             % Kilosort + 3Brain
    name = dataset_name(rec) + "_integrated";
else
    name = dataset_name(rec) + "_curated";
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = has_crop(rec)

ok = ~isempty(rec.crop_ellipse) || strlength(rec.crop_excel) > 0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function recs = read_recordings(file, base_dir)

% function to read the recordings table (all the columns are read as text, and converted here)

opts = detectImportOptions(file, 'TextType', 'string', 'VariableNamingRule', 'preserve');
opts = setvartype(opts, opts.VariableNames, 'string');
T = readtable(file, opts);
names = string(T.Properties.VariableNames);
required = ["folder", "file_name", "group"];
missing_col = required(~ismember(required, names));
if ~isempty(missing_col)
    error("The recordings table must have the columns %s (missing: %s)", strjoin(required, ', '), strjoin(missing_col, ', '));
end

recs = struct('row', {}, 'folder', {}, 'file_name', {}, 'group', {}, 'to_discard', {}, 'organoid_label', {}, 'crop_ellipse', {}, ...
    'crop_excel', {}, 'file_name_3Brain', {}, 'source', {}, 'include', {});
for r=1:height(T)
    rec.row = r;
    rec.folder = cell_text(T, r, "folder");
    if strlength(rec.folder) == 0 && strlength(cell_text(T, r, "file_name")) == 0
        continue                                                            % empty row
    end
    if ~is_absolute(char(rec.folder))
        rec.folder = string(fullfile(base_dir, rec.folder));
    end
    rec.file_name = cell_text(T, r, "file_name");
    rec.group = cell_text(T, r, "group");
    if strlength(rec.file_name) == 0 || strlength(rec.group) == 0
        error("Row %d of the recordings table: the columns file_name and group cannot be empty", r);
    end
    rec.to_discard = numbers(cell_text(T, r, "to_discard"));
    rec.organoid_label = cell_text(T, r, "organoid_label");
    rec.crop_ellipse = numbers(cell_text(T, r, "crop_ellipse"));
    rec.crop_excel = cell_text(T, r, "crop_excel");
    if ~isempty(rec.crop_ellipse) && numel(rec.crop_ellipse) ~= 4
        error("Row %d of the recordings table: crop_ellipse must have 4 numbers (x0 y0 radiusx radiusy)", r);
    end
    if has_crop(rec) && strlength(rec.organoid_label) == 0
        error("Row %d of the recordings table: the cropping needs an organoid_label", r);
    end
    rec.file_name_3Brain = cell_text(T, r, "file_name_3Brain");
    rec.source = "Kilosort";
    if strcmpi(cell_text(T, r, "source"), "3Brain")
        rec.source = "3Brain";
    end
    include = lower(cell_text(T, r, "include"));
    rec.include = ~ismember(include, ["0", "no", "false"]);
    recs(end+1) = rec;                                                      %#ok<AGROW>
end
recs = recs([recs.include]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function value = cell_text(T, r, name)

% text of the column "name" in the row r ("" if the column does not exist, or the cell is empty)

value = "";
if ismember(name, string(T.Properties.VariableNames))
    value = T.(name)(r);
    if ismissing(value)
        value = "";
    end
    value = strtrim(value);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = numbers(value)

% numbers written in a cell of the table (e.g., "3 17 25", "[3, 17, 25]"); [] if empty

x = [];
value = regexprep(char(value), '[\[\],;]', ' ');
if ~isempty(strtrim(value))
    x = sscanf(value, '%f')';
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function recs = select_rows(recs, rows)

% recordings selected with the option 'rows' (rows of the table, or file names)

if isempty(rows)
    return
end
if isnumeric(rows)
    recs = recs(ismember([recs.row], rows));
else
    recs = recs(ismember([recs.file_name], string(rows)));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = is_absolute(path_name)

ok = ~isempty(regexp(char(path_name), '^(/|~|[A-Za-z]:[\\/]|\\\\)', 'once'));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function txt = changes_description(changes)

% short description of the changes, for the messages and the history

txt = "";
if isempty(changes)
    return
end
if isstruct(changes)
    txt = "; with changes (structure)";
    return
end
list = changes;
if ~(size(changes,2) == 2 && size(changes,1) > 1)
    list = reshape(changes(:), 2, [])';
end
items = strings(size(list,1),1);
for k=1:size(list,1)
    value = list{k,2};
    if isnumeric(value) || islogical(value)
        value = mat2str(value);
    end
    items(k) = string(list{k,1}) + " = " + string(value);
end
txt = "; changes: " + strjoin(items, ", ");

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function code_ver = code_version(code_dir)

% version of the code (git commit), saved in the history of each recording

code_ver = "unknown";
[status, out] = system(sprintf('git -C "%s" rev-parse --short HEAD', code_dir));
if status == 0
    code_ver = string(strtrim(out));
    [status, out] = system(sprintf('git -C "%s" status --porcelain --untracked-files=no', code_dir));
    if status == 0 && ~isempty(strtrim(out))
        code_ver = code_ver + " (with uncommitted changes)";
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function txt = first_line(message)

lines = splitlines(string(message));
txt = lines(1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function txt = time_text(t)

if t < 60
    txt = sprintf('%.0f s', t);
elseif t < 3600
    txt = sprintf('%.1f min', t/60);
else
    txt = sprintf('%.1f h', t/3600);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function txt = now_text(format)

txt = string(datetime('now', 'Format', format));

end
