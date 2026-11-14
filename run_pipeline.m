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
%   run_pipeline("recordings.csv", 'parameters', "saved", 'renew', "cluster_neuronal_events")   % the same, but with the parameters of a
%                                                                                       % function revised since then taken from the file
%   run_pipeline("recordings.csv", 'changes', {'NetworkBurstDetection.method3.THRESHOLD', 5.0})   % a few values changed
%   run_pipeline("recordings.csv", 'rows', [2 3])                                       % only some recordings of the table

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   recordings:             path of the recordings table (.csv or .xlsx; see "recordings_template.csv"), one row per recording, with the columns
%       folder              folder of the recording (absolute, or relative to the folder of the table)
%       file_name           common name of the files of the recording (e.g., "recording_01"); for a 3Brain recording,
%                           the name of its files without the suffix (e.g., "XXX_3Brain")
%       group               group of the recording: the analyzed files of a group are collected in the folder <project>/<group>,
%                           and compiled together by "Compile_MEAproject"
%       
%       (optional columns)
%       to_discard          clusters discarded by "QualityControl_Manual" (e.g., "3 17 25"; the indexes are defined after eliminating
%                           the silent clusters); empty: none
%       organoid_label      label of the organoid, with crop_ellipse ("x0 y0 radiusx radiusy", in pixels) or crop_excel (file of the
%                           3Brain electrode map). With the photo of the chip: the labels of its organoids, separated by ";" (e.g.,
%                           "56; 57; 58"), given to them in the order of the photo (from the top down) and kept at the next runs (below);
%                           their number must be the number of organoids found (else the run does not start); empty: org1, org2...
%       picture             (November 2026) photo of the chip (file in the folder of the recording, or full path): "SelectElectrodes"
%                           keeps the electrodes under the organoid found in the photo ("OrganoidsFromPhoto", option 3, the default;
%                           crop_ellipse or crop_excel, when given, are used instead). Empty: the photo in the folder of the recording,
%                           when there is one (none: every electrode is kept; several: the run does not start, write the one to use);
%                           "none": no cropping
%       picture_organoid    (optional) number of the organoid in the photo (as numbered in the figure of OrganoidsFromPhoto(picture), from
%                           the top of the array down): only that organoid, with one label. Empty: every organoid of the photo is
%                           cropped, each as its own dataset "<file_name>_<label>" (one row of the table, run once per organoid); the
%                           labels are saved with the place of each organoid on the electrode map in "<file_name>_organoids.mat" (folder
%                           of the recording), and the next runs give them to the same organoids
%       picture_3D_side     (optional) side of the photo where the "3D" text of the chip is (right, top, left, bottom); empty: found
%                           from the photo (when given, it is checked against the pattern of the chip around the array)
%       picture_n_organoids (optional) number of organoids on the chip; empty: found from the photo
%       file_name_3Brain    name of a 3Brain dataset to integrate with the Kilosort dataset ("IntegrateDataSets")
%       source              "Kilosort" (default) or "3Brain" (step by step only)
%       include             =0, skip the recording (default: 1)
%
%   options (name-value pairs):
%       'parameters'        parameter file (default: "pipeline_parameters.m" of this folder), or "saved" (the set saved for each
%                           recording by its last run)
%       'renew'             (new, November 2026) with 'parameters', "saved": parts of each saved set taken from "pipeline_parameters.m"
%                           of this folder instead, e.g., the parameters of a function revised after the last run (a saved set may hold
%                           parameters that the function no longer has, which stops the run): a whole function ("Function": the
%                           parameters it no longer has are removed) or single values ("Function.PARAMETER"); the rest of each saved
%                           set (e.g., the settings of a recording) is kept. They are renewed before the 'changes' (default: none)
%       'changes'           values changed for this run, on top of the parameter set: {'Function.PARAMETER', value, ...}
%       'stages'            stages to run (they always run in the order of the pipeline) (default: all, end to end)
%       'from'              first stage to run; the following ones also run
%       'rows'              recordings to run: rows of the table (1 = first recording) or file names (default: all, with include ~= 0)
%       'project'           folder for the group folders and the compiled files (default: the folder of the table)
%       'verbose'           =0, only the progress is printed (the messages of the stages go to the log file of each recording)
%                           =1, the messages of the stages are also printed (default: 1 for a single stage, 0 otherwise)
%       'stop_on_error'     =1, stop at the first error; =0, continue with the next recording (default: 0)
%       'extra_path'        folders to add to the MATLAB path (e.g., npy-matlab, needed by "ShapeSpikeTrainData")
%       'progress'          (new, November 2026) function called at each step of the run with one structure (field "type": "start",
%                           "stage", "recording", "group" or "end"; see "report_progress"), e.g., by the window "CellCircuit"; when it
%                           returns true, the run stops after the current stage (default: none)

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
    'stop_on_error', 0, 'extra_path', {{}}, 'renew', "", 'progress', []);  % (revised November 2026: 'renew', 'progress')
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
recs = read_recordings(table_file, table_dir);                              % local function
recs = select_rows(recs, opts.rows);                                        % local function
if isempty(recs)
    error("No recording to run (check the column 'include' and the option 'rows')");
end
[recs, notes] = assign_organoids(recs, ismember("SelectElectrodes", stages));   % (November 2026) photos found, organoids labeled
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
% (new, November 2026) parts of the saved sets taken from the parameter file (checked here, before any recording runs)
renew = string(opts.renew);
renew = renew(strlength(renew) > 0);
P_renew = [];
if ~isempty(renew)
    if ~use_saved
        error("The option 'renew' is used with 'parameters', ""saved"" (with a parameter file, every value already comes from the file)");
    end
    renew_file = string(fullfile(code_dir, 'pipeline_parameters.m'));
    P_renew = load_parameters(renew_file);
    check_renew(P_renew, renew);                                            % local function
    source_text = source_text + ", with " + strjoin(renew, ", ") + " from " + renew_file;
end
changes_text = changes_description(opts.changes);
if any(stages == "ShapeSpikeTrainData") && any([recs.source] ~= "3Brain") && exist('readNPY', 'file') ~= 2
    error("The function readNPY (npy-matlab), needed by ShapeSpikeTrainData, is not on the MATLAB path: add it, or use the option 'extra_path'");
end
code_ver = code_version(code_dir);
cleanup = onCleanup(@() parameter_set('clear'));                            % the set is deactivated at the end (also after an error)
setappdata(groot, 'run_pipeline_n_recordings', n_rec);                     % (new, November 2026) number of recordings of the run, read by
cleanup_n = onCleanup(@() setappdata(groot, 'run_pipeline_n_recordings', 0));   % "Analyze_spikeTrain" (HIDE_PLOTS_MIN_RECORDINGS: plots saved, not displayed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the stages, recording by recording

t_start = tic;
fprintf("\n=== run_pipeline (%s): %d recording(s); stages: %s\n", now_text('yyyy-MM-dd HH:mm'), n_rec, strjoin(stages, ' > '));
fprintf("    parameters: %s%s\n", source_text, changes_text);
if ~verbose
    fprintf("    the messages of the stages are saved in the log file of each recording (<file_name>_pipeline_log.txt)\n");
end
for note = notes'
    fprintf("    %s\n", note);
end

rec_stages = stages(stages ~= "Compile_MEAproject");
pre_crop = STAGES(1:find(STAGES == "SelectElectrodes") - 1);                % stages before the cropping (the same for every organoid)
done_pre = strings(0,1);                                                    % "row|stage" already run in this run
errors = strings(0,1);
ok_rec = true(n_rec,1);
% (new, November 2026) progress of the run, for the option 'progress' (e.g., the window "CellCircuit"); stopped: the user asked to stop
stopped = report_progress(opts.progress, 'type', "start", 'recordings', string(arrayfun(@display_name, recs, 'UniformOutput', false)), 'groups', [recs.group], 'stages', stages, ...
    'parameters', source_text + changes_text);
for r=1:n_rec
    if isempty(rec_stages) || stopped
        break
    end
    rec = recs(r);
    fprintf("\n[%d/%d] %s (group: %s) - %s\n", r, n_rec, display_name(rec), rec.group, now_text('HH:mm'));
    if end_to_end && rec.source == "3Brain"
        fprintf("      skipped: 3Brain recordings run step by step (option 'stages')\n");
        ok_rec(r) = false;
        stopped = report_progress(opts.progress, 'type', "recording", 'index', r, 'state', "skipped", 'message', "3Brain recordings run step by step");
        continue
    end
    file_params = fullfile(rec.folder, rec.file_name + "_parameters.mat");
    file_log = fullfile(rec.folder, rec.file_name + "_pipeline_log.txt");
    if ~isfolder(rec.folder)
        msg = sprintf("the folder %s does not exist", rec.folder);
        fprintf("      ERROR: %s\n", msg);
        errors(end+1) = sprintf("[%d/%d] %s: %s", r, n_rec, rec.file_name, msg);
        ok_rec(r) = false;
        stopped = report_progress(opts.progress, 'type', "recording", 'index', r, 'state', "error", 'message', string(msg));
        continue
    end

    % parameter set of the recording
    try
        if use_saved
            S = load_saved(file_params);
            P = load_parameters(renew_set(S.parameters, P_renew, renew), opts.changes);   % (revised November 2026: renew_set)
        else
            P = P_file;
        end
    catch ME
        fprintf("      ERROR: %s\n", ME.message);
        errors(end+1) = sprintf("[%d/%d] %s: %s", r, n_rec, rec.file_name, ME.message);
        ok_rec(r) = false;
        stopped = report_progress(opts.progress, 'type', "recording", 'index', r, 'state', "error", 'message', string(ME.message));
        continue
    end
    parameter_set('set', P);
    record = start_record(file_params, P, rec_stages, source_text, changes_text, opts.changes, code_ver);

    for st = rec_stages
        if rec.expanded && ismember(st, pre_crop) && ismember(string(rec.row) + "|" + st, done_pre)
            fprintf("      %-28s already run for this recording (another organoid of its photo)\n", st);
            report_progress(opts.progress, 'type', "stage", 'index', r, 'stage', st, 'state', "skipped", 'message', "already run for this recording");
            continue
        end
        [inputs, skip] = stage_inputs(st, rec, end_to_end, verbose);
        if strlength(skip) > 0
            fprintf("      %-28s skipped (%s)\n", st, skip);
            report_progress(opts.progress, 'type', "stage", 'index', r, 'stage', st, 'state', "skipped", 'message', skip);
            continue                                                        % (a stop asked meanwhile is seen at the end of the next stage)
        end
        report_progress(opts.progress, 'type', "stage", 'index', r, 'stage', st, 'state', "running", 'message', "");
        [ME, elapsed] = run_stage(st, inputs, verbose, file_log);
        if isempty(ME)
            fprintf("%s\n", done_text(elapsed, verbose, st));
            record.parameters_stage.(st) = P;
            record.history(end).stages_done(end+1) = st;
            if ismember(st, pre_crop), done_pre(end+1) = string(rec.row) + "|" + st; end %#ok<AGROW>
            if st == "Analyze_spikeTrain"                                   % collect the analyzed file in the folder of its group
                collect_analyzed(rec, project);
            end
            stopped = report_progress(opts.progress, 'type', "stage", 'index', r, 'stage', st, 'state', "done", 'elapsed', elapsed, 'message', "");
        else
            fprintf("%s\n", error_text(ME, elapsed, verbose, st, file_log));
            record.history(end).status = "error in " + st + ": " + first_line(ME.message);
            errors(end+1) = sprintf("[%d/%d] %s, %s: %s", r, n_rec, rec.file_name, st, first_line(ME.message));
            ok_rec(r) = false;
            stopped = report_progress(opts.progress, 'type', "stage", 'index', r, 'stage', st, 'state', "error", 'elapsed', elapsed, ...
                'message', first_line(ME.message), 'log', string(file_log));
        end
        save_record(file_params, record);
        if ~isempty(ME)
            if opts.stop_on_error
                rethrow(ME);
            end
            break                                                           % the next stages depend on this one
        end
        if stopped && st ~= rec_stages(end)                                 % (new, November 2026) the user asked to stop: the next stages do not run
            fprintf("      stopped by the user after %s\n", st);
            record.history(end).status = "stopped by the user after " + st;
            save_record(file_params, record);
            ok_rec(r) = false;
            break
        end
    end
    if ok_rec(r)
        record.history(end).status = "done";
        save_record(file_params, record);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compile the results, group by group

if any(stages == "Compile_MEAproject") && ~stopped                         % (revised November 2026: not after a stop asked by the user)
    groups = unique([recs.group], 'stable');
    for g = groups
        if stopped
            break
        end
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
                P = load_parameters(renew_set(S.parameters, P_renew, renew), opts.changes);   % (revised November 2026: renew_set)
            else
                P = P_file;
            end
        catch ME
            fprintf("      ERROR: %s\n", ME.message);
            errors(end+1) = sprintf("[group %s]: %s", g, ME.message);      %#ok<AGROW>
            stopped = report_progress(opts.progress, 'type', "group", 'group', g, 'state', "error", 'message', string(ME.message));
            continue
        end
        parameter_set('set', P);
        record = start_record(file_params, P, "Compile_MEAproject", source_text, changes_text, opts.changes, code_ver);
        inputs = struct('folderPath', char(project), 'group', char(g), 'system', 0);
        report_progress(opts.progress, 'type', "group", 'group', g, 'state', "running", 'message', "");
        [ME, elapsed] = run_stage("Compile_MEAproject", inputs, verbose, file_log);
        if isempty(ME)
            fprintf("%s\n", done_text(elapsed, verbose, "Compile_MEAproject"));
            record.parameters_stage.Compile_MEAproject = P;
            record.history(end).stages_done(end+1) = "Compile_MEAproject";
            record.history(end).status = "done";
            stopped = report_progress(opts.progress, 'type', "group", 'group', g, 'state', "done", 'elapsed', elapsed, 'message', "", ...
                'excel', string(fullfile(project, g + "_Compiled.xlsx")));
        else
            fprintf("%s\n", error_text(ME, elapsed, verbose, "Compile_MEAproject", file_log));
            record.history(end).status = "error in Compile_MEAproject: " + first_line(ME.message);
            errors(end+1) = sprintf("[group %s] Compile_MEAproject: %s", g, first_line(ME.message));   %#ok<AGROW>
            stopped = report_progress(opts.progress, 'type', "group", 'group', g, 'state', "error", 'elapsed', elapsed, ...
                'message', first_line(ME.message), 'log', string(file_log));
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
if stopped
    fprintf("    (stopped by the user: the remaining stages, recordings and compilations did not run)\n");
end
report_progress(opts.progress, 'type', "end", 'errors', errors, 'stopped', stopped, 'elapsed', toc(t_start), 'project', project);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%   LOCAL FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stop = report_progress(fn, varargin)

% (new, November 2026) function to report a step of the run to the function given with the option 'progress' (e.g., the window "CellCircuit"),
% as one structure: type "start" (recordings, groups, stages, parameters), "stage" (index of the recording, stage, state: running / done /
% skipped / error, elapsed, message, log), "recording" (index, state: error / skipped, message: the recording could not start), "group"
% (group, state: running / done / error, elapsed, message, excel, log) or "end" (errors, stopped, elapsed, project). It returns true when
% that function returns true (the user asked to stop). An error of that function (e.g., the window was closed) does not stop the run

stop = false;
if isempty(fn)
    return
end
try
    out = fn(struct(varargin{:}));
    stop = isequal(out, true);
catch ME
    fprintf("      (the progress display failed: %s)\n", ME.message);
end
end

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
        elseif ~isempty(rec.crop_ellipse)
            inputs.option = 2;      inputs.ellipse_cropped = rec.crop_ellipse;
        else                                                                % (November 2026) photo of the chip
            inputs.option = 3;      inputs.file_photo = char(rec.picture);  inputs.photo_organoid = rec.picture_organoid;
            inputs.photo_text_side = rec.picture_3D_side;   inputs.photo_n_organoids = rec.picture_n_organoids;
            inputs.photo_position = rec.photo_position;                     % (saved place of the organoid: it keeps its label)
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
function check_renew(P_renew, renew)

% (new, November 2026) function to check the names of the option 'renew': a function of the parameter file ("Function") or one of
% its parameters ("Function.PARAMETER", "Function.method3.PARAMETER"), not a linked parameter (it is copied from its source)

links = cell(0,2);
if isfield(P_renew, 'links')
    links = P_renew.links;
end
for k=1:numel(renew)
    parts = strsplit(renew(k), '.');
    S = P_renew;
    found = parts(1) ~= "links";
    for j=1:numel(parts)
        if ~found || ~isstruct(S) || ~isfield(S, char(parts(j)))
            found = false;
            break
        end
        S = S.(char(parts(j)));
    end
    if ~found
        error("'%s' (option 'renew') is not in the parameter file (check the name: 'Function' or 'Function.PARAMETER')", renew(k));
    end
    if any(strcmp(links(:,1), renew(k)))
        source_name = links{strcmp(links(:,1), renew(k)), 2};
        error("The parameter '%s' (option 'renew') is linked to '%s' (the same value): renew '%s' instead", renew(k), source_name, source_name);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = renew_set(P, P_renew, renew)

% (new, November 2026) function to replace parts of a saved set by those of the parameter file (option 'renew'): a whole function
% (its parameters that are not in the file are removed) or single values; with no part to renew, the saved set is unchanged

for k=1:numel(renew)
    parts = cellstr(strsplit(renew(k), '.'));
    P = setfield(P, parts{:}, getfield(P_renew, parts{:}));
end

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

ok = ~isempty(rec.crop_ellipse) || strlength(rec.crop_excel) > 0 || (strlength(rec.picture) > 0 && rec.picture ~= "none" ...
    && ~isempty(rec.picture_organoid));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out, notes] = assign_organoids(recs, cropping)

% (new, November 2026) photo of the chip and organoids of each recording cropped from its photo (no crop_ellipse, no crop_excel):
%   - picture empty: the photo in the folder of the recording ("recording_picture"; the one of the saved labels, if any); none: no
%     cropping (every electrode); several, when the cropping runs: the run does not start (write the one to use). "none": no cropping
%   - picture_organoid empty: every organoid of the photo is cropped, each as its own dataset (one recording per organoid, all from the
%     same row of the table), with the labels of organoid_label ("56; 57; 58", given to the organoids in the order of the photo, from the
%     top down; their number must be the number of organoids found, else the run does not start), or org1, org2... when it is empty. The
%     labels are saved in "<file_name>_organoids.mat" with the place of each organoid on the electrode map, and given again to the same
%     organoids at the next runs: when the photo or its settings (picture_3D_side, picture_n_organoids) change, the organoids are found
%     again and each one takes the label of the saved organoid at its place; when the labels of the row change, they are given again in
%     the order of the photo. Without the cropping (later stages only), the saved labels give the cropped datasets; none saved: the
%     recording is not cropped
%   - picture_organoid given: that organoid only, with one label (organoid_label, or org<number>)
% cropping: =true when "SelectElectrodes" runs. notes: what was done, for the messages of the run

out = recs([]);     notes = strings(0, 1);
for i = 1:numel(recs)
    rec = recs(i);
    given = split_labels(rec.organoid_label);                               % labels written in the row
    photo = isempty(rec.crop_ellipse) && strlength(rec.crop_excel) == 0 && rec.picture ~= "none";
    if ~photo || ~isempty(rec.picture_organoid)
        if numel(given) > 1
            error(['Row %d of the recordings table (%s): several labels in organoid_label (%s) are given to the organoids found in the ' ...
                'photo of the chip; with crop_ellipse, crop_excel or picture_organoid, write one label (one row per organoid)'], rec.row, ...
                rec.file_name, rec.organoid_label);
        end
        if photo && isempty(given)
            rec.organoid_label = "org" + rec.picture_organoid;
        end
        out(end+1) = rec; %#ok<AGROW>
        continue
    end
    file_saved = fullfile(rec.folder, rec.file_name + "_organoids.mat");
    saved = [];
    if isfile(file_saved)
        S = load(file_saved, 'organoids');  saved = S.organoids;
    end
    relabel = ~isempty(given) && ~isempty(saved) && ~same_labels(given, saved.labels);  % the labels of the row changed
    if ~cropping                                                            % later stages: the datasets cropped before (saved labels)
        if isempty(saved)
            out(end+1) = rec; %#ok<AGROW>
            continue
        end
        if relabel
            error(['Row %d of the recordings table (%s): the labels in organoid_label (%s) are not those of the last cropping (%s): run the ' ...
                'cropping again (from SelectElectrodes)'], rec.row, rec.file_name, strjoin(given, '; '), strjoin(saved.labels, '; '));
        end
        pic = saved.picture_path;
    else
        pic = rec.picture;                                                  % the photo
        if strlength(pic) == 0
            if ~isempty(saved) && isfile(saved.picture_path)
                pic = saved.picture_path;
            else
                [p, candidates] = recording_picture(rec.folder, rec.file_name);
                if strlength(p) > 0
                    pic = string(fullfile(rec.folder, p));
                elseif numel(candidates) > 1
                    error(['Row %d of the recordings table (%s): several photos in the folder of the recording (%s): write the one to ' ...
                        'use in the column picture (or "none": no cropping)'], rec.row, rec.file_name, strjoin(candidates, ', '));
                end
            end
        end
        if strlength(pic) == 0                                              % no photo: every electrode
            if numel(given) > 1
                error(['Row %d of the recordings table (%s): several labels in organoid_label (%s), but no photo of the chip in the ' ...
                    'folder of the recording (column picture)'], rec.row, rec.file_name, rec.organoid_label);
            end
            out(end+1) = rec; %#ok<AGROW>
            continue
        end
        [~, pname, pext] = fileparts(pic);  pic_name = string(pname) + string(pext);
        same = ~isempty(saved) && saved.picture == pic_name && saved.picture_3D_side == rec.picture_3D_side ...
            && saved.picture_n_organoids == rec.picture_n_organoids;
        if ~same || relabel
            if same                                                         % the same organoids, new labels (in the order of the photo)
                found = struct('text_side', saved.text_side);
                centroids = saved.centroids;    radius = saved.radius;
            else                                                            % organoids found (again)
                try
                    found = OrganoidsFromPhoto(char(pic), struct('FIGURE', false, 'TEXT_SIDE', rec.picture_3D_side, ...
                        'N_ORGANOIDS', rec.picture_n_organoids));
                catch ME
                    error(['Row %d of the recordings table (%s): the organoids could not be found in the photo %s (%s); give crop_ellipse, ' ...
                        'or "none" in the column picture (no cropping)'], rec.row, rec.file_name, pic_name, ME.message);
                end
                centroids = reshape([found.organoids.centroid], 2, [])';
                radius = reshape([found.organoids.radius], [], 1);
            end
            n = size(centroids, 1);
            if ~isempty(given) && numel(given) ~= n
                error(['Row %d of the recordings table (%s): %d organoid(s) found in the photo %s, but %d label(s) in organoid_label ' ...
                    '(%s). Check the labels; if the photo shows %d organoids, write their number in picture_n_organoids (the photo was ' ...
                    'misread)'], rec.row, rec.file_name, n, pic_name, numel(given), strjoin(given, '; '), numel(given));
            end
            keep_saved = saved;
            if relabel, keep_saved = []; end                                % (new labels: in the order of the photo)
            labels = labels_by_place(centroids, radius, keep_saved, given);
            organoids = struct('picture', pic_name, 'picture_path', string(pic), 'picture_3D_side', rec.picture_3D_side, ...
                'picture_n_organoids', rec.picture_n_organoids, 'text_side', found.text_side, 'labels', labels, ...
                'numbers', (1:n)', 'centroids', centroids, 'radius', radius, 'date', string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
            save(file_saved, 'organoids');
            notes(end+1) = sprintf('%s: %d organoid(s) in the photo %s, labeled %s (saved in %s)', rec.file_name, n, pic_name, ...
                strjoin(labels + " (number " + string((1:n)') + ")", ", "), rec.file_name + "_organoids.mat"); %#ok<AGROW>
            saved = organoids;
        end
    end
    for k = 1:numel(saved.labels)                                           % one recording per organoid
        e = rec;
        e.picture = string(pic);
        e.picture_organoid = saved.numbers(k);
        e.organoid_label = saved.labels(k);
        e.photo_position = saved.centroids(k, :);
        e.expanded = true;
        out(end+1) = e; %#ok<AGROW>
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function labels = split_labels(text)

% labels written in organoid_label ("56; 57; 58", or with commas); none: empty

labels = strtrim(split(string(text), [";", ","]));
labels = labels(strlength(labels) > 0);
labels = labels(:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = same_labels(a, b)

% =true if a and b hold the same labels (in any order)

ok = numel(a) == numel(b) && all(ismember(a, b)) && all(ismember(b, a));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function labels = labels_by_place(centroids, radius, saved, pool)

% labels of the organoids found (centroids [x y], radius, in electrodes): each one takes the label of the saved organoid at its place
% (closest pairs first, within the radius of either); the others take, in the order of the photo, the labels of pool (labels written in
% the row) not used yet, or new labels when pool is empty (org<n>, after the largest number used)

n = size(centroids, 1);     labels = strings(n, 1);     used = strings(0, 1);
if ~isempty(saved)
    used = saved.labels(:);
    D = hypot(centroids(:,1) - saved.centroids(:,1)', centroids(:,2) - saved.centroids(:,2)');
    while any(isfinite(D(:)))
        [d, j] = min(D(:));     [a, b] = ind2sub(size(D), j);
        if d <= max(radius(a), saved.radius(b)) && (isempty(pool) || ismember(saved.labels(b), pool))
            labels(a) = saved.labels(b);
        end
        D(a, :) = inf;  D(:, b) = inf;
    end
end
if ~isempty(pool)
    free = pool(~ismember(pool, labels));
    for a = find(labels == "")'
        if isempty(free), break, end
        labels(a) = free(1);    free(1) = [];
    end
else
    num = str2double(regexp([used; labels(labels ~= "")], '\d+$', 'match', 'once'));
    next = max([0; num(~isnan(num))]) + 1;
    for a = find(labels == "")'
        labels(a) = "org" + next;   next = next + 1;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function name = display_name(rec)

% name of a recording in the messages (with the label of its organoid when one row of the table gives several organoids)

name = rec.file_name;
if rec.expanded
    name = rec.file_name + " (" + rec.organoid_label + ")";
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function recs = read_recordings(file, base_dir)

% function to read the recordings table (all the columns are read as text, and converted here)
% (revised November 2026) a .csv table is always read with commas as the delimiter: with file names full of "_" and empty optional
% columns, MATLAB guessed "_" as the delimiter and no column was found
% (revised November 2026) the rows that are not included (include = 0, no or false) are skipped before they are checked, as the window
% "CellCircuit" does: it saves every row, also the unticked ones, which may have no group (e.g., a long list where only the few recordings
% to run were given a group). Before, an unticked row without a group or file name stopped the run

if endsWith(lower(string(file)), ".csv")
    opts = detectImportOptions(file, 'TextType', 'string', 'VariableNamingRule', 'preserve', 'FileType', 'text', 'Delimiter', ',');
else
    opts = detectImportOptions(file, 'TextType', 'string', 'VariableNamingRule', 'preserve');
end
opts = setvartype(opts, opts.VariableNames, 'string');
T = readtable(file, opts);
names = string(T.Properties.VariableNames);
required = ["folder", "file_name", "group"];
missing_col = required(~ismember(required, names));
if ~isempty(missing_col)
    error("The recordings table must have the columns %s (missing: %s)", strjoin(required, ', '), strjoin(missing_col, ', '));
end

recs = struct('row', {}, 'folder', {}, 'file_name', {}, 'group', {}, 'to_discard', {}, 'organoid_label', {}, 'crop_ellipse', {}, ...
    'crop_excel', {}, 'picture', {}, 'picture_organoid', {}, 'picture_3D_side', {}, 'picture_n_organoids', {}, 'file_name_3Brain', {}, ...
    'source', {}, 'include', {}, 'photo_position', {}, 'expanded', {});
for r=1:height(T)
    if ismember(lower(cell_text(T, r, "include")), ["0", "no", "false"])
        continue                                                            % a row that is not included (not checked; see above)
    end
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
    rec.picture = cell_text(T, r, "picture");                               % (November 2026) photo of the chip
    if lower(rec.picture) == "none"
        rec.picture = "none";                                               % (no cropping)
    elseif strlength(rec.picture) > 0 && ~is_absolute(char(rec.picture))
        rec.picture = string(fullfile(rec.folder, rec.picture));            % (a file in the folder of the recording)
    end
    rec.picture_organoid = numbers(cell_text(T, r, "picture_organoid"));
    if ~isempty(rec.picture_organoid) && (~isscalar(rec.picture_organoid) || rec.picture_organoid < 1 || rec.picture_organoid ~= round(rec.picture_organoid))
        error("Row %d of the recordings table: picture_organoid must be one number (the organoid in the photo: 1, 2...)", r);
    end
    rec.picture_3D_side = lower(cell_text(T, r, "picture_3D_side"));
    if ~ismember(rec.picture_3D_side, ["", "right", "top", "left", "bottom"])
        error("Row %d of the recordings table: picture_3D_side must be right, top, left or bottom (or empty: found from the photo)", r);
    end
    rec.picture_n_organoids = numbers(cell_text(T, r, "picture_n_organoids"));
    if isempty(rec.picture_n_organoids)
        rec.picture_n_organoids = 0;                                        % (found from the photo)
    elseif ~isscalar(rec.picture_n_organoids) || rec.picture_n_organoids < 1 || rec.picture_n_organoids ~= round(rec.picture_n_organoids)
        error("Row %d of the recordings table: picture_n_organoids must be one number (or empty: found from the photo)", r);
    end
    if strlength(rec.picture) > 0 && rec.picture ~= "none" && ~isfile(rec.picture)
        error("Row %d of the recordings table: the picture %s does not exist", r, rec.picture);
    end
    if (~isempty(rec.crop_ellipse) || strlength(rec.crop_excel) > 0) && strlength(rec.organoid_label) == 0
        error("Row %d of the recordings table: the cropping needs an organoid_label", r);
    end
    rec.photo_position = [];    rec.expanded = false;                       % (set by "assign_organoids")
    rec.file_name_3Brain = cell_text(T, r, "file_name_3Brain");
    rec.source = "Kilosort";
    if strcmpi(cell_text(T, r, "source"), "3Brain")
        rec.source = "3Brain";
    end
    rec.include = true;                                                     % (the rows not included were skipped above)
    recs(end+1) = rec;                                                      %#ok<AGROW>
end

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
elseif isfile(fullfile(code_dir, 'CellCircuit_version.txt'))               % (new, November 2026) installed toolbox (no git): the version written
    code_ver = string(strtrim(fileread(fullfile(code_dir, 'CellCircuit_version.txt'))));   % by "tools/package_CellCircuit"
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
