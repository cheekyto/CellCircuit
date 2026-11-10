function passed = test_baseline_condition(code_dir)

% Test of the baseline identified by its file name ("baselineCondition", used by "Analyze_spikeTrain"; new, November 2026), without data:
%   1. "baselineCondition" on lists of file names: usual order, TTX listed first, three conditions in every order, the extension .bxr,
%      a single condition, a cell array, two identical names, and names that do not fit the rule (3Brain names, no common start: the
%      first condition, as before)
%   2. "Analyze_spikeTrain" on a synthetic recording (8 units) whose TTX condition is listed first: the baseline is saved as "Baseline"
%      (with its own window) and the TTX condition as "+TTX"; with the usual order, the labels are the same as before
%
% Example (from a terminal):
%   matlab -batch "addpath('<repo>/tests'); test_baseline_condition('<repo>')"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dir:             folder with the code to be tested

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% passed:               =true if all the checks passed (each check is printed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(char(code_dir));
passed = true;

%% 1. baselineCondition
three = ["rec_D3.brw", "rec_D3_01.brw", "rec_D3_ttx_00.brw"];
orders = perms(1:3);
ok_three = true;
for i=1:size(orders,1)
    [b, by_name] = baselineCondition(three(orders(i,:)));
    ok_three = ok_three && by_name && three(orders(i,b)) == "rec_D3.brw";
end
cases = {
    % file names                                                    baseline    by name     description
    ["rec_D3.brw", "rec_D3_ttx_00.brw"],                            1,          true,       'usual order';
    ["rec_D3_ttx_00.brw", "rec_D3.brw"],                            2,          true,       'TTX listed first';
    ["rec_D3_ttx__00.brw", "rec_D3.brw"],                           2,          true,       'TTX listed first (other spelling)';
    ["rec_D3.brw", "rec_D3_baseline2.brw", "rec_D3_ttx_00.brw"],    1,          true,       'second baseline file';
    ["s_org.6_D4_ttx_00.brw", "s_org.6_D4.brw"],                    2,          true,       'dots in the name';
    ["rec_ttx.bxr", "rec.bxr"],                                     2,          true,       'extension .bxr';
    "rec_D3.brw",                                                   1,          true,       'single condition';
    {'rec_D3_ttx_00.brw'; 'rec_D3.brw'},                            2,          true,       'cell array (column)';
    ["rec_D3.brw", "rec_D3.brw"],                                   1,          true,       'identical names: the first one';
    ["rec_Baseline", "rec_TTX"],                                    1,          false,      '3Brain names: the first condition';
    ["rec_TTX", "rec_Baseline"],                                    1,          false,      '3Brain names, other order: the first condition (as before)';
    ["recA.brw", "recB_ttx.brw"],                                   1,          false,      'no common start: the first condition'};
ok_cases = true;
for i=1:size(cases,1)
    [b, by_name] = baselineCondition(cases{i,1});
    if b ~= cases{i,2} || by_name ~= cases{i,3}
        ok_cases = false;
        fprintf('   failed: %s (baseline %d, by name %d)\n', cases{i,4}, b, by_name);
    end
end
passed = report(passed, ok_three && ok_cases, sprintf('1. baselineCondition: %d lists of file names and the 6 orders of three conditions', size(cases,1)));

%% 2. Analyze_spikeTrain on a synthetic recording, TTX condition listed first
folder = tempname;   mkdir(folder);   cleanup = onCleanup(@() rmdir(folder, 's'));
rng(1);
n_units = 8;
win_ttx = [0 300];   win_base = [300 1200];                                 % in sec: the TTX segment comes first in the concatenated data
spikeTimestamp = struct('timeStamp', cell(1,n_units), 'nUnits', [], 'recLength', [], 'window', [], 'cluster', [], 'position', [], ...
    't_step', [], 'positionXY', []);
for i=1:n_units
    t_ttx = sort(win_ttx(1) + rand(round(0.05*diff(win_ttx)),1)*diff(win_ttx));         % ~0.05 Hz during TTX
    t_base = sort(win_base(1) + rand(round(2*diff(win_base)),1)*diff(win_base));        % ~2 Hz during the baseline
    spikeTimestamp(i).timeStamp = [t_ttx; t_base]*1000;                     % in msec
    spikeTimestamp(i).nUnits = numel(spikeTimestamp(i).timeStamp);
    spikeTimestamp(i).recLength = win_base(2) - win_ttx(1);
    spikeTimestamp(i).window = [win_ttx(1) win_base(2)];
    spikeTimestamp(i).cluster = i;
    spikeTimestamp(i).position = i;
    spikeTimestamp(i).t_step = 0.05;
    spikeTimestamp(i).positionXY = [i 1];
end
n_clusterOriginal = n_units;   n_cluster_Kilosort = n_units;   n_cluster_3Brain = 0;
labels = strings(2,2);   files = strings(2,2);   windows = cell(2,2);
for variant = 1:2                                                           % 1: TTX listed first; 2: usual order
    rec_conditions = {'rec_D3_ttx_00.brw', win_ttx; 'rec_D3.brw', win_base};
    if variant == 2
        rec_conditions = rec_conditions([2 1], :);
    end
    name = sprintf('synthetic%d_curated', variant);
    save(fullfile(folder, [name '_QCmanual.mat']), 'spikeTimestamp', 'rec_conditions', 'n_clusterOriginal', 'n_cluster_Kilosort', 'n_cluster_3Brain');
    inputs = struct('folders', string(folder), 'files_data', string(name), 'condition_analysis', -1, 'window', 600, 'opt_rerun', 0, ...
        'print_Plots', 0, 'plot_option', 0, 'plot_network_burst', 0, 'opt_checking', 0, 'system', 0);
    evalc('Analyze_spikeTrain(inputs)');
    S = load(fullfile(folder, [name '_Analyzed.mat']), 'PerCondition');
    for k=1:2
        labels(variant,k) = string(S.PerCondition(k).condition);
        files(variant,k) = string(S.PerCondition(k).file);
        windows{variant,k} = S.PerCondition(k).window;
    end
end
ok_first = isequal(labels(1,:), ["+TTX", "Baseline"]) && files(1,2) == "rec_D3.brw" && isequal(windows{1,2}, [600 1200]) && isequal(windows{1,1}, [0 300]);
ok_usual = isequal(labels(2,:), ["Baseline", "+TTX"]) && files(2,1) == "rec_D3.brw" && isequal(windows{2,1}, [600 1200]);
passed = report(passed, ok_first && ok_usual, sprintf('2. Analyze_spikeTrain: TTX listed first -> %s / %s; usual order -> %s / %s', labels(1,1), labels(1,2), labels(2,1), labels(2,2)));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function passed = report(passed, ok, text)
if ok
    fprintf('PASS  %s\n', text);
else
    fprintf('FAIL  %s\n', text);
    passed = false;
end
end
