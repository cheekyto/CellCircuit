function passed = test_plotFeatures(code_dir)

% Function-level test of "plotFeatures" (November 2026) on synthetic compiled files (3 groups, tables Firing and NETbursting0, a baseline
% and a gabazine condition, known values). It checks that
%   1. "data": the groups (file names without "_Compiled"), the tables found (sheet names of the Excel file), their measures (numeric
%      columns, without ttx_qc, T_start, T_end), their conditions (condition and dose; the baseline first) and the units of the measures
%      (first row of the Excel sheet next to the .mat file)
%   2. bar: one panel per measure, the values of each group and condition (NaN left out and counted as missing), the summary = direct
%      computation (n, mean, SD, SEM, median, percentiles), the bar at the mean and an error bar of one SEM (or SD), one point per value,
%      the names of the groups under the plot, no legend with one condition; groups chosen by name or number, in the order given, each
%      with its own color (the same when another group is left out); new names for the groups
%   3. several conditions: one block of groups per condition, the names of the conditions under the plot, a legend of the groups
%   4. violin and box: a violin (or box) per group, the median, the mean, the whiskers at the percentiles chosen (min-max with [0 100])
%   5. pages: at most "per_page" panels per figure; the PDF has one page per figure; a .png per page; the preview in a panel (first page)
%   6. wrong choices stop with a clear message
% No reference recording is needed
%
% Example (from a terminal):
%   matlab -batch "addpath('<repo>/tests'); test_plotFeatures('<repo>')"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dir:             folder with the code to be tested

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% passed:               =true if all the checks passed (each check is printed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(char(code_dir));
passed = true;
root = tempname;    mkdir(root);
cleanup = onCleanup(@() rmdir(root, 's'));

% synthetic compiled files: group A (4 organoids), B (3 organoids, one NaN rate at baseline), C (2 organoids, no NETbursting0 table)
rate = {[1 2 3 4], [2 NaN 6], [5 7]};                                       % Firing_Rate at baseline
rate_gbz = {[0.5 1 1.5 2], [1 2 3], [2 4]};                                 % Firing_Rate with gabazine
names = ["A", "B", "C"];
files = strings(3, 1);
for g=1:3
    n = numel(rate{g});
    org = "org" + (1:n)';
    Organoid = cellstr([org; org]);     File = cellstr([org + "_base"; org + "_gbz"]);
    Condition = [repmat({'Baseline'}, n, 1); repmat({'+Gabazine'}, n, 1)];
    dose = [repmat({''}, n, 1); repmat({'20uM'}, n, 1)];
    ttx_qc = ones(2*n, 1);  T_start = zeros(2*n, 1);    T_end = 600*ones(2*n, 1);
    Firing_Rate = [rate{g}(:); rate_gbz{g}(:)];
    CV_ISI = (1:2*n)'/10 + g;
    Firing = table(Organoid, ttx_qc, File, Condition, dose, T_start, T_end, Firing_Rate, CV_ISI);
    files(g) = string(fullfile(root, names(g) + "_Compiled.mat"));
    if g < 3
        Rate_NETburst = [10*g + (1:n)'; 5*g + (1:n)'];
        NETbursting0 = table(Organoid, File, ttx_qc, Condition, dose, Rate_NETburst);
        save(files(g), 'Firing', 'NETbursting0');
    else
        save(files(g), 'Firing');
    end
end
xls = fullfile(root, "A_Compiled.xlsx");                                    % units in the first row, as "Compile_MEAproject"
writetable(Firing, xls, 'Sheet', 'Firing', 'Range', 'A2');
writecell({'', '', '', '', '', 'sec', 'sec', 'Hz', 'ad'}, xls, 'Sheet', 'Firing', 'Range', 'A1');

% 1. data
[f0, D] = plotFeatures(files, "data");
S1 = D.sheets(1);   S2 = D.sheets(2);
ok = isempty(f0) && isequal(D.groups, names') && isequal([D.sheets.name], ["Firing", "Network Bursting Pooled"]) && ...
    isequal([D.sheets.var], ["Firing", "NETbursting0"]) && isequal(S1.measures, ["Firing_Rate", "CV_ISI"]) && ...
    isequal(S1.conditions, ["Baseline", "+Gabazine 20uM"]) && isequal(S1.units, ["Hz", "ad"]) && isequal(S2.measures, "Rate_NETburst") && ...
    isequal(S2.units, "");
passed = report("data: groups, tables (sheet names), measures (not ttx_qc / T_start / T_end), conditions (baseline first), units from the Excel sheet", ok) && passed;

% 2. bar, one condition
[f, R] = plotFeatures(D, "bar", struct('visible', false, 'measures', "Firing_Rate"));
ax = findobj(f, 'Type', 'axes');
sm = R.summary;
v = rate{2}(isfinite(rate{2}));
rowB = sm(sm.Group == "B", :);
bars = findobj(ax, 'Type', 'patch');
pts = findobj(ax, 'Type', 'scatter');
ok = isscalar(ax) && height(sm) == 3 && isequal(sm.N', [4 2 2]) && isequal(sm.N_missing', [0 1 0]) && ...
    abs(rowB.Mean - mean(v)) < 1e-12 && abs(rowB.SD - std(v)) < 1e-12 && abs(rowB.SEM - std(v)/sqrt(2)) < 1e-12 && ...
    abs(rowB.Median - median(v)) < 1e-12 && abs(sm.P25(1) - prctile(rate{1}, 25)) < 1e-12 && abs(sm.Whisker_high(1) - prctile(rate{1}, 95)) < 1e-12 && ...
    height(R.values) == 8 && isequal(sort(R.values.Value(R.values.Group == "A"))', rate{1}) && numel(bars) == 3 && ...
    sum(arrayfun(@(p) numel(p.XData), pts)) == 8 && isequal(string(ax.XTickLabel(:))', names) && isempty(findobj(f, 'Type', 'legend')) && ...
    string(ax.YLabel.String) == "Firing_Rate (Hz)";
heights = sort(arrayfun(@(p) max(p.YData), bars))';
ok = ok && max(abs(heights - sort([mean(rate{1}), mean(v), mean(rate{3})]))) < 1e-12;
err_lines = findobj(ax, 'Type', 'line');
tops = arrayfun(@(h) max(h.YData), err_lines);
ok = ok && any(abs(tops - (mean(rate{1}) + std(rate{1})/2)) < 1e-12);    % (SEM of A: std/sqrt(4))
close(f);
passed = report("bar: values per group (NaN left out, counted), summary = direct computation, bars at the means, SEM error bars, one point per value, group names", ok) && passed;
[f, R] = plotFeatures(D, "bar", struct('visible', false, 'measures', "Firing_Rate", 'error_bar', "sd", 'groups', ["C", "A"]));
ax = findobj(f, 'Type', 'axes');
tops = arrayfun(@(h) max(h.YData), findobj(ax, 'Type', 'line'));
colA = get(findobj(ax, 'Type', 'patch', '-and', 'XData', 2 + 0.36*[-1; 1; 1; -1]), 'FaceColor');
close(f);
[f, ~] = plotFeatures(D, "bar", struct('visible', false, 'measures', "Firing_Rate"));
colA_all = get(findobj(findobj(f, 'Type', 'axes'), 'Type', 'patch', '-and', 'XData', 1 + 0.36*[-1; 1; 1; -1]), 'FaceColor');
close(f);
[f, R2] = plotFeatures(D, "bar", struct('visible', false, 'measures', "Firing_Rate", 'groups', [2 1], 'labels', ["WT", "KO", "Het"]));
lab2 = string(findobj(f, 'Type', 'axes').XTickLabel(:))';
close(f);
ok = isequal(R.groups', ["C", "A"]) && any(abs(tops - (mean(rate{1}) + std(rate{1}))) < 1e-12) && isequal(colA, colA_all) && ...
    isequal(lab2, ["KO", "WT"]) && isequal(unique(R2.values.Group)', ["KO", "WT"]);
passed = report("bar: SD error bars; groups by name or number in the order given; a group keeps its color when another is left out; new names", ok) && passed;

% 3. several conditions
[f, R] = plotFeatures(D, "bar", struct('visible', false, 'measures', ["Firing_Rate", "CV_ISI"], 'conditions', "all"));
ax = findobj(f, 'Type', 'axes');
lg = findobj(f, 'Type', 'legend');
gb = R.summary(R.summary.Measure == "Firing_Rate" & R.summary.Condition == "+Gabazine 20uM" & R.summary.Group == "B", :);
ok = numel(ax) == 2 && isequal(R.conditions, ["Baseline", "+Gabazine 20uM"]) && isscalar(lg) && isequal(string(lg.String), names) && ...
    all(arrayfun(@(a) isequal(string(a.XTickLabel(:))', ["Baseline", "+Gabazine 20uM"]), ax)) && abs(gb.Mean - 2) < 1e-12 && ...
    height(R.summary) == 2*2*3;
close(f);
passed = report("several conditions: one block of groups per condition (names under the plot), a legend of the groups, summary per condition", ok) && passed;

% 4. violin and box
[f, R] = plotFeatures(D, "violin", struct('visible', false, 'measures', "Firing_Rate"));
ax = findobj(f, 'Type', 'axes');
n_violin = numel(findobj(ax, 'Type', 'patch'));
diamonds = findobj(ax, 'Type', 'line', 'Marker', 'diamond');
dots = findobj(ax, 'Type', 'line', 'Marker', 'o');
ok = n_violin == 3 && numel(diamonds) == 3 && numel(dots) == 3 && ...
    any(arrayfun(@(h) abs(h.YData - mean(rate{1})) < 1e-12, diamonds)) && any(arrayfun(@(h) abs(h.YData - median(rate{1})) < 1e-12, dots));
close(f);
[f, R] = plotFeatures(D, "box", struct('visible', false, 'measures', "Firing_Rate", 'whiskers', [0 100], 'groups', "A"));
ax = findobj(f, 'Type', 'axes');
box_patch = findobj(ax, 'Type', 'patch');
ylines = arrayfun(@(h) h.YData(:)', findobj(ax, 'Type', 'line'), 'UniformOutput', false);
ys = unique(round([ylines{:}], 12));
q = prctile(rate{1}, [25 75]);
ok = ok && isscalar(box_patch) && abs(min(box_patch.YData) - q(1)) < 1e-12 && abs(max(box_patch.YData) - q(2)) < 1e-12 && ...
    any(abs(ys - 1) < 1e-12) && any(abs(ys - 4) < 1e-12) && any(abs(ys - median(rate{1})) < 1e-12) && R.summary.Whisker_low == 1 && R.summary.Whisker_high == 4;
close(f);
passed = report("violin and box: one per group, median and mean marked, box from the 25th to the 75th percentile, whiskers at min-max with [0 100]", ok) && passed;

% 5. pages, files, preview
many = struct('visible', false, 'measures', ["Firing_Rate", "CV_ISI"], 'per_page', 1, 'output', fullfile(root, "feat.pdf"));
[f, R] = plotFeatures(D, "bar", many);
n_fig = numel(f);   close(f);
many.output = fullfile(root, "feat.png");
[f, ~] = plotFeatures(D, "bar", many);     close(f);
p = uipanel(figure('Visible', 'off'));
[fp, Rp] = plotFeatures(D, "box", struct('measures', ["Firing_Rate", "CV_ISI"], 'per_page', 1, 'parent', p));
ok = n_fig == 2 && R.n_pages == 2 && isfile(fullfile(root, "feat.pdf")) && isfile(fullfile(root, "feat.png")) && isfile(fullfile(root, "feat_2.png")) && ...
    fp == p && Rp.n_pages == 1 && isscalar(findobj(p, 'Type', 'axes'));
close(ancestor(p, 'figure'));
passed = report("pages: one figure per 'per_page' panels, PDF and one .png per page; the preview in a panel (first page)", ok) && passed;

% 6. wrong choices
ok = throws(@() plotFeatures(D, "pie"), "the plot must be") && ...
    throws(@() plotFeatures(D, "bar", struct('sheet', "Spikes")), "not in these compiled files") && ...
    throws(@() plotFeatures(D, "bar", struct('measures', "Rate")), "not a measure") && ...
    throws(@() plotFeatures(D, "bar", struct('conditions', "+TTX")), "condition not found") && ...
    throws(@() plotFeatures(D, "bar", struct('groups', "Z")), "group not found") && ...
    throws(@() plotFeatures(D, "box", struct('whiskers', [30 70])), "whiskers must be") && ...
    throws(@() plotFeatures(D, "bar", struct('error_bar', "ci")), "error_bar must be") && ...
    throws(@() plotFeatures(fullfile(root, "none_Compiled.mat"), "data"), "file not found");
passed = report("wrong choices stop with a clear message", ok) && passed;

if passed
    fprintf('test_plotFeatures: all checks passed\n');
else
    fprintf('test_plotFeatures: SOME CHECKS FAILED\n');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = throws(fn, text)
% true when fn stops with an error whose message contains text
try
    fn();
    ok = false;
catch err
    ok = contains(err.message, text);
    if ~ok
        fprintf('        (message: %s)\n', err.message);
    end
end
end

function ok = report(name, ok)
if ok
    fprintf("  PASS  %s\n", name);
else
    fprintf("  FAIL  %s\n", name);
end
end
