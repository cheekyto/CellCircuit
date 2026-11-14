function passed = test_compare_neurons(code_dir)

% Function-level test of "CompareNeurons" on a synthetic project (new, November 2026)
% Three groups of 8 organoids (A, B and C; C only helps to fit the axes), with known effects: the units of B are broader (higher waveform
% score 1), the rate of the network bursts follows the mean of that score of each organoid (so the group difference in rate is explained
% by the composition), and the broader units start the network bursts earlier. One organoid of C has no saved neuron scores and no
% templates (scored here, firing pattern only), the TTX condition is listed before the baseline, and the recordings table has a column
% "batch". It checks that
%   1. the axes are fitted on the 24 organoids and saved; the files of the recordings are not changed; the baseline row of the compiled
%      tables is used (not the TTX row)
%   2. composition: B has a higher mean waveform score 1 than A
%   3. activity: the composition explains the group difference in the rate of the network bursts (effect of the composition > 0, group
%      effect smaller once the composition is in the model), with the random intercept of "batch"; the yes/no measure gets a logistic model;
%      "joint" puts the three scores in one model (with the similarity to the QC noise shape)
%   4. roles: broader units start earlier (onset ~ waveform score 1 < 0), across the organoids, in each group, and their difference
%   5. a custom list of activity measures replaces the default one; the defaults of the code are those of "comparison_settings.m"
%   6. nothing runs with fewer than 20 organoids, a compared group outside FIT_GROUPS, an unknown group, an unknown setting or a wrong
%      activity measure
% No reference recording is needed
%
% Example (from a terminal):
%   matlab -batch "addpath('<repo>/tests'); test_compare_neurons('<repo>')"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dir:             folder with the code to be tested

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% passed:               =true if all the checks passed (each check is printed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(char(code_dir));
parameter_set('clear');
passed = true;
rng(20261113);
tmp = tempname;
mkdir(tmp);
cleanup = onCleanup(@() cleanup_test(tmp));
project = fullfile(tmp, 'project');
out = fullfile(tmp, 'results');
info = make_compare_project(project);

C.TABLES = string(fullfile(project, 'recordings.csv'));
C.COMPARISONS = {["A", "B"]};
C.RANDOM = "batch";
C.OUTPUT = string(out);
C.NAME = "test";
[~, Res] = evalc('CompareNeurons(C)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Axes, files, baseline

O = Res.Organoids;
v = whos('-file', info.no_scores_file);
k = find(O.Organoid == info.no_scores_name);
passed = report(sprintf("axes fitted on the 24 organoids of A, B and C and saved (agreement with the fixed axes %.3f, %.3f, %.3f); " + ...
    "results written", Res.axesAgreement), height(O) == 24 && Res.axes.fit.nGroups == 3 && Res.axes.fit.nOrganoidsWaveform == 23 && ...
    isfile(Res.files.axes) && isfile(Res.files.excel) && isfile(Res.files.mat) && contains(Res.axes.version, "test") && ...
    all(Res.axesAgreement > 0.9)) && passed;
passed = report("an organoid without saved scores is scored here (firing pattern only), and its file is not changed", ...
    O.ScoresSource(k) == "computed here, firing pattern only" && ~any(strcmp({v.name}, 'NeuronScore')) && isnan(O.MeanWaveformScore1(k)) && ...
    sum(O.ScoresSource == "saved") == 23) && passed;
i = find(O.Organoid == info.names(3));
passed = report("the baseline row of the compiled tables is used (the TTX condition is listed first)", ...
    abs(O.("Activity_Network_burstRate")(i) - info.rate(3)) < 1e-12 && abs(O.Activity_FiringRate(i) - info.firing(3)) < 1e-12) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Composition

P = Res.Composition;
r = P(P.Score == "MeanWaveformScore1", :);
passed = report(sprintf("composition: B has a higher mean waveform score 1 than A (%.2f SD of the organoid means, p = %.2g)", ...
    r.DifferenceInSD, r.p), height(P) == 3 && r.Difference > 0 && r.p < 0.05 && r.Organoids == 16 && contains(r.Model, "batch") && ...
    P.Covariates(P.Score == "MeanWaveformScore2") == "MeanNoiseSimilarity") && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Activity

A = Res.Activity;
r = A(A.Activity == "Network-burst rate" & A.CompositionTerm == "MeanWaveformScore1", :);
passed = report(sprintf("activity: the composition explains the group difference in network-burst rate (group effect %.2f -> %.2f, " + ...
    "composition %.2f, p = %.2g; %d bursting organoids)", r.GroupEffect, r.GroupEffectAdjusted, r.CompositionEffect, r.CompositionP, r.N), ...
    height(A) == 15 && r.CompositionEffect > 0 && r.CompositionP < 0.01 && r.GroupP < 0.05 && abs(r.GroupEffectAdjusted) < abs(r.GroupEffect) / 2 && ...
    r.GroupChangePercent > 50 && r.CompositionTestP < 0.01 && r.N == info.n_bursting && contains(r.Model, "linear mixed model") && ...
    contains(r.Model, "batch") && r.Organoids == "bursting") && passed;
b = A(A.Activity == "Bursting", :);
passed = report("the yes/no measure (bursting) gets a logistic model on all the organoids", all(startsWith(b.Model, "logistic")) && ...
    all(b.N == 16)) && passed;

C2 = C;  C2.COMPOSITION_MODE = "joint";  C2.RANDOM = [];  C2.ROLES = [];  C2.NAME = "joint";
[~, R2] = evalc('CompareNeurons(C2)');
j = R2.Activity(R2.Activity.Activity == "Network-burst rate", :);
passed = report("joint: the three mean scores in one model, with the similarity to the QC noise shape; ordinary regression without " + ...
    "random intercepts", height(j) == 3 && all(j.Composition == "MeanWaveformScore1 + MeanWaveformScore2 + MeanFiringPatternScore") && ...
    all(j.Covariates == "MeanNoiseSimilarity") && all(j.Model == "linear regression") && isempty(R2.Roles)) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Roles

L = Res.Roles;
on = L(L.Role == "onset" & L.Score == "waveformScore1", :);
all_groups = on(on.Groups == "all groups", :);
passed = report(sprintf("roles: broader units start earlier (onset ~ waveform score 1: %.2f, p = %.2g, %d organoids), in A and in B, " + ...
    "and the difference B - A is estimated", all_groups.Estimate, all_groups.p, all_groups.Organoids), all_groups.Estimate < 0 && ...
    all_groups.p < 0.01 && all(on.Estimate(on.Groups == "A" | on.Groups == "B") < 0) && any(on.Groups == "B - A") && ...
    ~isnan(on.Estimate(on.Groups == "B - A")) && all_groups.Organoids == info.n_role_organoids && ...
    numel(unique(L.Role)) == 2) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Custom activity measures, defaults

C3 = C;  C3.ACTIVITY = ["Burst rate", "NETbursting1.Rate_NETburst", "log", "bursting"];  C3.ROLES = [];  C3.NAME = "custom";
[~, R3] = evalc('CompareNeurons(C3)');
D = CompareNeurons("defaults");
T = template_settings(fullfile(char(code_dir), 'comparison_settings.m'));
same = true;
for f=string(fieldnames(D))'
    if ismember(f, ["TABLES", "COMPARISONS"]), continue, end
    same = same && isfield(T, f) && isequal(string(T.(f)(:)), string(D.(f)(:)));
end
passed = report("a custom list of activity measures replaces the default one; the defaults of the code are those of comparison_settings.m", ...
    all(R3.Activity.Activity == "Burst rate") && height(R3.Activity) == 3 && same && isequal(sort(fieldnames(T)), sort(fieldnames(D)))) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Refusals

few = C;  few.FIT_GROUPS = ["A", "B"];
outside = C;  outside.FIT_GROUPS = ["A", "B", "C"];  outside.COMPARISONS = {["A", "D"]};
helper = C;  helper.FIT_GROUPS = ["A", "C"];
typo = C;  typo.COVARIATE = "LogUnits";
bad = C;  bad.ACTIVITY = ["Rate", "Rate_NETburst", "log", "all"];
refused = [fails(@() quiet(few), "at least 20"), fails(@() quiet(outside), "no organoids for the group"), ...
    fails(@() quiet(helper), "among FIT_GROUPS"), fails(@() quiet(typo), "unknown setting"), fails(@() quiet(bad), "<table>.<column>")];
passed = report(sprintf("nothing runs with fewer than 20 organoids, an unknown group, a compared group outside FIT_GROUPS, an unknown " + ...
    "setting or a wrong activity measure (%d of 5 refused)", sum(refused)), all(refused)) && passed;

if passed
    fprintf("test_compare_neurons: all checks passed\n");
else
    fprintf("test_compare_neurons: SOME CHECKS FAILED\n");
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quiet(S)
% CompareNeurons without its report
evalc('CompareNeurons(S)');
end

function S = template_settings(file)
C = struct();
run(file);
S = C;
end

function failed = fails(f, message)
% =true if f stops with an error whose message contains the text given
failed = false;
try
    f();
catch err
    failed = contains(err.message, message);
    if ~failed
        fprintf("        (other error: %s)\n", err.message);
    end
end
end

function cleanup_test(folder)
parameter_set('clear');
if isfolder(folder)
    rmdir(folder, 's');
end
end

function ok = report(name, ok)
if ok
    fprintf("  PASS  %s\n", name);
else
    fprintf("  FAIL  %s\n", name);
end
end
