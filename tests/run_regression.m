function summary = run_regression(stage, varargin)

% Regression test for the pipeline: run the OLD and the NEW version of the code on the same input data (a reference recording),
% compare all the output variables and measure the running time
% Each version runs in a fresh MATLAB process, on its own copy of the data (the reference recording is never modified)
% (revised October 2026)
%
% Example:      make_baseline();                                     % once: snapshot of the current code, used as the OLD version
%               run_regression("QualityControl_General");            % after modifying the code
%               run_regression("all");                               % the whole pipeline (stage by stage)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stage:                name of the pipeline function to test (see the list STAGES below), or "all"
% name-value pairs (optional):
%   'old':              folder with the OLD version of the code (default: the baseline saved by "make_baseline")
%   'new':              folder with the NEW version of the code (default: the repository folder, i.e., the current version)
%   'profile':          =1, run the MATLAB profiler (default: 1); the profiler slows down the code, but shows where the time goes
%                       =0, no profiler (more accurate timing); ="old" or ="new", only for one of the versions
%   'only':             "old" or "new", to run only one of the versions (no comparison); default: "" (both)
%   any field of "regression_config" (e.g., 'n_iteration_hs', 'work_root')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% summary:              structure with the running time (old / new, in sec), the speed-up, and the result of the comparison
%                       the full report (+ profiler summary) is saved in <work_root>/<stage>/report.txt and summary.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

STAGES = ["ShapeSpikeTrainData", "QualityControl_General", "QualityControl_HyperSynch", "QualityControl_Manual", "IntegrateDataSets", "Analyze_spikeTrain", "Compile_MEAproject"];
N_PROFILE_ROWS = 25;                                                        % number of functions listed in the profiler summary

cfg = regression_config();
opts = struct('old', "", 'new', cfg.repo_dir, 'profile', 1, 'only', "");
for k=1:2:numel(varargin)                                                   % name-value pairs: options or fields of the configuration
    name = char(varargin{k});
    if isfield(opts, name)
        opts.(name) = varargin{k+1};
    elseif isfield(cfg, name)
        cfg.(name) = varargin{k+1};
    else
        error('Unknown option "%s"', name);
    end
end
if strlength(opts.old) == 0
    opts.old = fullfile(cfg.work_root, 'code_baseline');
end
if ~isfolder(opts.old)
    error('The OLD version of the code was not found in %s (run "make_baseline" first)', opts.old);
end

if strcmpi(stage, "all")                                                    % the whole pipeline, stage by stage
    for s=1:numel(STAGES)
        summary(s) = run_regression(STAGES(s), varargin{:});               %#ok<AGROW>
    end
    return
end
if ~ismember(stage, STAGES)
    error('Unknown stage "%s"', stage);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare the input data (the same for both versions)

prepare_chain(stage, cfg, opts.old);                                        % the inputs of the stage are generated with the OLD version (if they do not exist yet)
stage_dir = fullfile(cfg.work_root, stage);
versions = ["old", "new"];
if strlength(opts.only) > 0
    versions = string(opts.only);
end
code_dirs = struct('old', opts.old, 'new', opts.new);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run each version in a fresh MATLAB process

for v=1:numel(versions)
    version = versions(v);
    data_dir = fullfile(stage_dir, version);
    do_profile = isequal(opts.profile, 1) || strcmpi(string(opts.profile), version);    % the profiler can be used for both versions or only for one of them
    [~, out_files] = run_one(stage, code_dirs.(version), data_dir, cfg, do_profile);
    res = load(fullfile(data_dir, 'result.mat'), 'result');
    results.(version) = res.result;
    outputs.(version) = out_files;
    fprintf('%s (%s): %.1f sec\n', stage, version, results.(version).elapsed);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare the outputs and write the report

summary = struct('stage', string(stage), 'elapsed_old', NaN, 'elapsed_new', NaN, 'speedup', NaN, 'same', false, 'report', strings(0,1));
lines = sprintf("=== Regression test: %s (%s)", stage, string(datetime('now')));
lines(end+1) = sprintf("OLD code: %s", opts.old);
lines(end+1) = sprintf("NEW code: %s", opts.new);
for v=1:numel(versions)
    version = versions(v);
    summary.(sprintf('elapsed_%s', version)) = results.(version).elapsed;
    if strlength(results.(version).error) > 0
        lines(end+1) = sprintf("!!! The %s version FAILED:", version);   %#ok<AGROW>
        lines(end+1) = results.(version).error;                            %#ok<AGROW>
    end
end

if numel(versions) == 2
    summary.speedup = summary.elapsed_old / summary.elapsed_new;
    lines(end+1) = sprintf("Elapsed time: old %.1f sec, new %.1f sec (speed-up x%.2f)", summary.elapsed_old, summary.elapsed_new, summary.speedup);
    all_same = strlength(results.old.error) == 0 && strlength(results.new.error) == 0;
    for f=1:numel(outputs.old)                                              % across all the output files of the stage
        file_old = fullfile(stage_dir, "old", outputs.old(f));
        file_new = fullfile(stage_dir, "new", outputs.old(f));
        if ~isfile(file_old) || ~isfile(file_new)
            lines(end+1) = sprintf("Output %s: MISSING (old: %d, new: %d)", outputs.old(f), isfile(file_old), isfile(file_new));  %#ok<AGROW>
            all_same = false;
            continue
        end
        [same, report, n_diff] = compare_outputs(file_old, file_new, cfg.tolerance, cfg.max_report);
        if same
            lines(end+1) = sprintf("Output %s: IDENTICAL", outputs.old(f));                 %#ok<AGROW>
        else
            lines(end+1) = sprintf("Output %s: %g DIFFERENCES", outputs.old(f), n_diff);    %#ok<AGROW>
            lines = [lines, "    " + report'];                                              %#ok<AGROW>
        end
        all_same = all_same && same;
    end
    summary.same = all_same;
else
    lines(end+1) = sprintf("Elapsed time (%s): %.1f sec", versions, summary.(sprintf('elapsed_%s', versions)));
end

for v=1:numel(versions)                                                     % profiler summary, for each version
    version = versions(v);
    if ~isempty(results.(version).profile)
        lines(end+1) = sprintf("--- Profiler (%s version): top %g functions by self time (sec)", version, N_PROFILE_ROWS);   %#ok<AGROW>
        T = results.(version).profile;
        for r=1:min(N_PROFILE_ROWS, height(T))
            lines(end+1) = sprintf("%-45s %-4s calls %9g   total %9.2f   self %9.2f", T.name(r), string(repmat('*', 1, T.pipeline(r))), T.calls(r), T.total(r), T.self(r));  %#ok<AGROW>
        end
        H = results.(version).hot_lines;
        lines(end+1) = sprintf("--- Profiler (%s version): slowest lines of the pipeline functions", version);   %#ok<AGROW>
        for r=1:min(15, height(H))
            lines(end+1) = sprintf("%-35s line %5g   calls %9g   time %9.2f", H.name(r), H.line(r), H.calls(r), H.time(r));   %#ok<AGROW>
        end
    end
end
summary.report = lines';
writelines(lines', fullfile(stage_dir, 'report.txt'));
save(fullfile(stage_dir, 'summary.mat'), 'summary', 'results');
fprintf('%s\n', lines);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data_dir, out_files] = run_one(stage, code_dir, data_dir, cfg, do_profile)

% function to run one stage with one version of the code: copy the input data, write the inputs, launch a fresh MATLAB process

[inputs, in_files, out_files, sub_dir] = stage_spec(stage, cfg, data_dir);
if isfolder(data_dir), rmdir(data_dir, 's'); end                            % start from a clean copy of the data
mkdir(fullfile(data_dir, sub_dir));
for f=1:numel(in_files)                                                     % copy the input files
    copyfile(in_files(f), fullfile(data_dir, sub_dir));
end
extra_path = {cfg.npy_folder};                                              %#ok<NASGU> other libraries needed by the pipeline
inputs_file = fullfile(data_dir, 'inputs.mat');
save(inputs_file, 'inputs', 'extra_path');

tests_dir = fileparts(mfilename('fullpath'));
result_file = fullfile(data_dir, 'result.mat');
log_file = fullfile(data_dir, 'log.txt');
command = sprintf('"%s" -batch "addpath(''%s''); run_stage(''%s'', ''%s'', ''%s'', ''%s'', %d)" > "%s" 2>&1', ...
    cfg.matlab_exe, tests_dir, code_dir, stage, inputs_file, result_file, do_profile, log_file);
system(command);                                                            % run the stage in a fresh MATLAB process
if ~isfile(result_file)
    error('The MATLAB process did not finish properly; see the log file %s', log_file);
end
out_files = fullfile(sub_dir, out_files);                                   % output files, relative to "data_dir"

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inputs, in_files, out_files, sub_dir] = stage_spec(stage, cfg, data_dir)

% function to define, for each stage: the inputs (to overwrite the values in the code), the input files and the output files
% The input files come from the reference recording or from the chain of intermediate results (generated with the OLD version)

chain = fullfile(cfg.work_root, 'chain');
base = cfg.base_name;
sub_dir = "";
switch stage
    case "ShapeSpikeTrainData"
        inputs = struct('folders', string(data_dir), 'files_output', base);
        in_files = [fullfile(chain, "spike_times.npy"), fullfile(chain, "spike_clusters.npy"), fullfile(cfg.ref_folder, "metadata.csv"), fullfile(cfg.ref_folder, "recording_data.csv")];
        out_files = base + "_raw.mat";

    case "QualityControl_General"
        inputs = struct('folder', char(data_dir), 'condition', -1, 'option_plot', 0, 'condition_3Brain', 0, 'option_picture', 0, 'option_pause', 0);
        in_files = fullfile(chain, base + "_raw.mat");
        out_files = base + "_cleaned.mat";

    case "QualityControl_HyperSynch"
        inputs = struct('folders', string(data_dir), 'files_data', base, 'N_ITERATION', cfg.n_iteration_hs);
        in_files = [fullfile(chain, base + "_raw.mat"), fullfile(chain, base + "_cleaned.mat")];
        out_files = base + "_QC.mat";

    case "QualityControl_Manual"
        inputs = struct('folder', char(data_dir), 'file_data', char(base), 'option_main', 3, 'to_discard', [], 'option_saveComplete', 1);
        in_files = [fullfile(chain, base + "_cleaned.mat"), fullfile(chain, base + "_QC.mat")];
        out_files = base + "_QCmanual.mat";

    case "IntegrateDataSets"                                                % Kilosort dataset only (the reference recording has no 3Brain dataset)
        inputs = struct('folder', char(data_dir), 'file_data_Kilosort', char(base), 'file_data_3Brain', [], 'include_Kilosort', 1, 'include_3Brain', 0, ...
            'option_plot', 0, 'option_print_list', 0);
        in_files = [fullfile(chain, base + "_QCmanual.mat"), fullfile(cfg.ref_folder, "metadata.csv")];
        out_files = base + "_curated_QCmanual.mat";

    case "Analyze_spikeTrain"                                               % it uses the curated file of the reference recording (generated by "IntegrateDataSets"; fixed input)
        inputs = struct('folders', string(data_dir), 'files_data', base + "_curated", 'condition_analysis', cfg.analyze_conditions, ...
            'opt_rerun', 0, 'opt_checking', 0, 'print_Plots', 0, 'plot_option', 0, 'plot_network_burst', 0);
        in_files = fullfile(cfg.ref_folder, base + "_curated_QCmanual.mat");
        out_files = base + "_curated_Analyzed.mat";

    case "Compile_MEAproject"                                               % the analyzed file must be inside the folder of the group
        sub_dir = "group1";
        inputs = struct('folderPath', char(data_dir), 'group', 'group1');
        in_files = fullfile(chain, base + "_curated_Analyzed.mat");
        out_files = "../group1_Compiled.mat";
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prepare_chain(stage, cfg, old_dir)

% function to generate (only once) the input files of the stage, running the upstream stages with the OLD version of the code

chain = fullfile(cfg.work_root, 'chain');
if ~isfolder(chain), mkdir(chain); end
base = cfg.base_name;

switch stage                                                                % upstream stage and the file it produces
    case "ShapeSpikeTrainData"
        if ~isfile(fullfile(chain, "spike_times.npy"))                      % convert the spike-sorting output (.mat) into the format read by "ShapeSpikeTrainData" (.npy)
            addpath(cfg.npy_folder);
            S = load(fullfile(cfg.ref_folder, "spike_times.mat"));      writeNPY(S.spike_times, fullfile(chain, "spike_times.npy"));
            S = load(fullfile(cfg.ref_folder, "spike_clusters.mat"));   writeNPY(S.spike_clusters, fullfile(chain, "spike_clusters.npy"));
        end
        return
    case "QualityControl_General",      upstream = "ShapeSpikeTrainData";       product = base + "_raw.mat";
    case "QualityControl_HyperSynch",   upstream = "QualityControl_General";    product = base + "_cleaned.mat";
    case "QualityControl_Manual",       upstream = "QualityControl_HyperSynch"; product = base + "_QC.mat";
    case "IntegrateDataSets",           upstream = "QualityControl_Manual";     product = base + "_QCmanual.mat";
    case "Analyze_spikeTrain",          return                              % fixed input from the reference recording
    case "Compile_MEAproject",          upstream = "Analyze_spikeTrain";        product = base + "_curated_Analyzed.mat";
end

if ~isfile(fullfile(chain, product))                                        % the input of the stage was not generated yet
    prepare_chain(upstream, cfg, old_dir);                                  % first, the inputs of the upstream stage
    data_dir = fullfile(cfg.work_root, upstream, 'old');                    % a previous regression test of the upstream stage may have generated it already
    reuse = isfile(fullfile(data_dir, product)) && isfile(fullfile(data_dir, 'result.mat'));
    if reuse
        res = load(fullfile(data_dir, 'result.mat'), 'result');
        reuse = strcmp(res.result.code_dir, old_dir);                       % only if it was generated with the same OLD version
    end
    if reuse
        fprintf('Using %s generated by the last regression test of %s (OLD version)\n', product, upstream);
    else
        fprintf('Generating %s with the OLD version (%s)...\n', product, upstream);
        data_dir = fullfile(cfg.work_root, 'chain_tmp');
        run_one(upstream, old_dir, data_dir, cfg, 0);
    end
    res = load(fullfile(data_dir, 'result.mat'), 'result');
    if strlength(res.result.error) > 0
        error('The stage %s failed with the OLD version:\n%s', upstream, res.result.error);
    end
    copyfile(fullfile(data_dir, product), chain);                           % keep the output, as input for the next stage
    copyfile(fullfile(data_dir, 'result.mat'), fullfile(chain, sprintf('result_%s.mat', upstream)));
    if upstream == "QualityControl_General" && cfg.hs_n_clusters > 0        % reduce the number of clusters, to keep "QualityControl_HyperSynch" within a reasonable time
        reduce_cleaned(fullfile(chain, product), cfg.hs_n_clusters);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function reduce_cleaned(file, n_keep)

% function to keep only the first "n_keep" clusters in the output of "QualityControl_General" (all its variables are kept consistent)

S = load(file);
keep = 1:min(n_keep, numel(S.spike_plot));                                  % clusters to keep (indexes before discarding)
n_good = sum(S.discarded_cluster(keep) == 0);                               % "spikeTimestamp" only has the good clusters, in the same order
S.spike_plot = S.spike_plot(keep);
S.spikeTimestamp = S.spikeTimestamp(1:n_good);
for c=1:numel(S.Indiv_Cond)
    S.Indiv_Cond{c} = S.Indiv_Cond{c}(keep);
end
per_cluster = {'removedSpikes', 'echo_cluster', 'silent_cluster', 'noisy_cluster', 'ttx_cluster', 'discarded_cluster', 'picture_cluster'};
for k=1:numel(per_cluster)
    S.(per_cluster{k}) = S.(per_cluster{k})(keep);
end
S.n_discarded = sum(S.discarded_cluster > 0);
S.n_cluster0 = numel(keep);
save(file, '-struct', 'S');
fprintf('Test data reduced to the first %g clusters (%g good clusters)\n', numel(keep), n_good);

end
