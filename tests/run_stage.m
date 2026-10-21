function run_stage(code_dir, stage, inputs_file, result_file, do_profile)

% function to run ONE stage of the pipeline with ONE version of the code (called by "run_regression", in a fresh MATLAB process)
% (revised October 2026)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dir:             folder with the version of the code to be used (only this folder provides the pipeline functions)
% stage:                name of the pipeline function to run (e.g., "QualityControl_HyperSynch")
% inputs_file:          .mat file with the variables "inputs" (structure for "override_inputs") and "extra_path" (other libraries, e.g., npy-matlab)
% result_file:          .mat file where the variable "result" is saved
% do_profile:           =1, run the MATLAB profiler; =0, only measure the elapsed time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% result.elapsed:       running time of the stage, in sec
% result.error:         error report (empty if the stage finished without error)
% result.profile:       table with the profiler summary (one row per function: calls, total time and self time, in sec)
% result.hot_lines:     table with the slowest lines of the pipeline functions (function, line, calls, time in sec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up the MATLAB path: the pipeline functions must come ONLY from "code_dir"

load(inputs_file, 'inputs', 'extra_path');
code_dir = char(code_dir);

sentinels = {'BurstDetection', 'QualityControl_General', 'override_inputs'};  % functions that identify a copy of the pipeline
for s=1:numel(sentinels)                                                    % remove any other copy of the pipeline from the MATLAB path
    hits = which('-all', sentinels{s});
    for h=1:numel(hits)
        folder_hit = fileparts(hits{h});
        if ismember(folder_hit, strsplit(path, pathsep))
            rmpath(folder_hit);
        end
    end
end
for k=1:numel(extra_path)                                                   % other libraries needed by the pipeline (e.g., npy-matlab)
    addpath(char(extra_path{k}));
end
addpath(code_dir);                                                          % version of the code under test
work_dir = tempname;    mkdir(work_dir);    cd(work_dir);                   % neutral current folder (no copy of the code can be found there)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the stage

inputs.system = 0;                                                          % the pipeline function must NOT modify the MATLAB path
result = struct('stage', string(stage), 'code_dir', string(code_dir), 'elapsed', NaN, 'error', "", 'profile', [], 'hot_lines', []);
if do_profile
    profile clear;  profile on;                                             % start the profiler
end
t0 = tic;                                                                   % start the timer
try
    feval(stage, inputs);                                                   % run the pipeline function
catch ME
    result.error = string(getReport(ME, 'extended', 'hyperlinks', 'off'));  % keep the error report (the comparison will fail)
end
result.elapsed = toc(t0);                                                   % running time, in sec
if do_profile
    profile off;
    [result.profile, result.hot_lines] = summarize_profile(profile('info'), code_dir);
end
close all force;
save(result_file, 'result');
if strlength(result.error) > 0
    fprintf(2, '%s\n', result.error);
end
fprintf('Elapsed time: %.1f sec\n', result.elapsed);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T, hot] = summarize_profile(info, code_dir)

% function to summarize the profiler output: one row per function (sorted by self time), and the slowest lines of the pipeline functions

N_HOT_LINES = 5;                                                            % number of slowest lines kept for each pipeline function

ft = info.FunctionTable;
n = numel(ft);
name = strings(n,1);    calls = zeros(n,1);     total = zeros(n,1);     self = zeros(n,1);  pipeline = false(n,1);
hot_fun = strings(0,1); hot_line = zeros(0,1);  hot_calls = zeros(0,1); hot_time = zeros(0,1);
for k=1:n                                                                   % across all profiled functions
    name(k) = ft(k).FunctionName;
    calls(k) = ft(k).NumCalls;
    total(k) = ft(k).TotalTime;
    child_time = 0;
    if ~isempty(ft(k).Children)
        child_time = sum([ft(k).Children.TotalTime]);
    end
    self(k) = total(k) - child_time;                                        % time spent in the function itself (without the functions it calls)
    pipeline(k) = strcmp(fileparts(ft(k).FileName), code_dir);              % it is a function of the pipeline (not a MATLAB or a test function)
    if pipeline(k) && ~isempty(ft(k).ExecutedLines)                         % slowest lines of the pipeline function
        lines = sortrows(ft(k).ExecutedLines, -3);                          % columns: line number, number of calls, time
        lines = lines(1:min(N_HOT_LINES, size(lines,1)), :);
        hot_fun = [hot_fun; repmat(name(k), size(lines,1), 1)];             %#ok<AGROW>
        hot_line = [hot_line; lines(:,1)];                                  %#ok<AGROW>
        hot_calls = [hot_calls; lines(:,2)];                                %#ok<AGROW>
        hot_time = [hot_time; lines(:,3)];                                  %#ok<AGROW>
    end
end
T = sortrows(table(name, pipeline, calls, total, self), 'self', 'descend');
hot = sortrows(table(hot_fun, hot_line, hot_calls, hot_time, 'VariableNames', {'name','line','calls','time'}), 'time', 'descend');

end
