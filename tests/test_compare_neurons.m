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
info = make_project(project);

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
function info = make_project(project)
% three groups of 8 organoids: recordings table, group folders with the "_Analyzed.mat" files, compiled tables
mkdir(project);
R = NeuronScores(struct('relativeRate', 0, 'lv', 1, 'cv', 1, 'burstIndex', 0), [], [], [], NaN).reference;
D.a1 = R.WF_AXES(:,1) / norm(R.WF_AXES(:,1));
a2 = R.WF_AXES(:,2) - (D.a1' * R.WF_AXES(:,2)) * D.a1;
D.a2 = a2 / norm(a2);
D.b = R.FP_AXIS / norm(R.FP_AXIS);
groups = ["A", "B", "C"];
shift = [0 1 0];                                                            % B: broader units
not_bursting = ["A_org1", "A_org2", "B_org1"];
rows = strings(0, 5);
info.names = strings(0, 1);  info.rate = [];  info.firing = [];
info.n_bursting = 0;  info.n_role_organoids = 0;
for g=1:3
    mkdir(fullfile(project, groups(g)));
    Firing = table();  NETbursting1 = table();
    for k=1:8
        name = groups(g) + "_org" + k;
        rec = name + "_curated";
        n = 60 + randi(40);
        z1 = shift(g) + 2 * randn(n, 1);
        [NS, ~] = synth(R, D, z1);
        nb = 20 * ~ismember(name, not_bursting);
        P = conditions(name, n, nb, z1);
        spikeTimestamp = struct('cluster', num2cell(1:n));                  %#ok<NASGU>
        file = fullfile(project, groups(g), rec + "_Analyzed.mat");
        if name == "C_org8"                                                 % no saved scores, no templates
            P(2).spikeTrain = arrayfun(@(u) struct('timeStamp', sort(rand(1, 300) * 600e3)), 1:n);
            P(2).window = [0 600];
            PerCondition = P;                                               %#ok<NASGU>
            save(file, 'PerCondition', 'spikeTimestamp');
            info.no_scores_file = file;  info.no_scores_name = rec;
        else
            NS.cluster = (1:n)';  NS.condition = 'Baseline';
            NeuronScore = NS;  PerCondition = P;                            %#ok<NASGU>
            save(file, 'PerCondition', 'spikeTimestamp', 'NeuronScore');
        end
        rate = NaN;  part = NaN;  dur = NaN;
        if nb > 0
            rate = 10^(0.3 + 0.5 * mean(z1) + 0.03 * randn);                % follows the composition
            part = 0.5 + 0.1 * rand;  dur = 0.5 + rand;
        end
        firing = 10^(0.2 * randn);
        if groups(g) ~= "C" && nb > 0
            info.n_bursting = info.n_bursting + 1;
            info.n_role_organoids = info.n_role_organoids + 1;
        end
        info.names(end+1) = rec;  info.rate(end+1) = rate;  info.firing(end+1) = firing;
        Firing = [Firing; table({char(rec); char(rec)}, {char(name + "_ttx_00.brw"); char(name + ".brw")}, {'+TTX'; 'Baseline'}, ...
            [0.01; firing], 'VariableNames', {'Organoid', 'File', 'Condition', 'Firing_Rate'})];   %#ok<AGROW>
        NETbursting1 = [NETbursting1; table({char(rec); char(rec)}, {char(name + "_ttx_00.brw"); char(name + ".brw")}, ...
            {'+TTX'; 'Baseline'}, [0; nb], [NaN; rate], [NaN; part], [NaN; dur], 'VariableNames', {'Organoid', 'File', 'Condition', ...
            'Number_NETburst_Events', 'Rate_NETburst', 'Participation_NETburst', 'Duration_NETburst'})];   %#ok<AGROW>
        rows(end+1, :) = [fullfile(project, "raw", name), name, groups(g), "b" + mod(k, 3), "1"];   %#ok<AGROW>
    end
    save(fullfile(project, groups(g) + "_Compiled.mat"), 'Firing', 'NETbursting1');
end
writetable(array2table(rows, 'VariableNames', {'folder', 'file_name', 'group', 'batch', 'include'}), fullfile(project, 'recordings.csv'));
end

function P = conditions(name, n, nb, z1)
% TTX (listed first) and baseline; primary network bursts of the baseline: 70% of the units take part, the broader ones start earlier
P = struct('file', {char(name + "_ttx_00.brw"), char(name + ".brw")}, 'condition', {'+TTX', 'Baseline'}, ...
    'NetworkBurstDetected', {struct('nBurst', 0), []}, 'NetworkBurstStructure', {[], []}, 'spikeTrain', {[], []}, 'window', {[], []});
N.nBurst = nb;
if nb > 0
    part = find(rand(n, 1) < 0.7);
    N.Channels = part;
    [~, order] = sort(-z1(part) + 0.5 * randn(numel(part), 1));            % broader units (higher z1) first
    rank = zeros(numel(part), 1);  rank(order) = 1:numel(part);
    N.timeStamp = cell(nb, n);
    for u=part'
        N.timeStamp(rand(nb, 1) < 0.8, u) = {1};
    end
    P(2).NetworkBurstStructure = struct('structure', [rank, rank * 5]);
else
    N.nBurst = NaN;                                                         % (rejected by the final check: no bursts)
    N.timeStamp = cell(0, n);
end
P(2).NetworkBurstDetected = N;
end

function [NS, z1] = synth(R, D, z1)
% scores of a synthetic organoid: waveform measures along D.a1 (latent z1) and D.a2, firing pattern along D.b
n = numel(z1);
T = R.WF_MEAN + R.WF_SD .* ([z1, 1.5 * randn(n, 1)] * [D.a1 D.a2]' + 0.3 * randn(n, numel(R.WF_MEAN)) + 0.1 * randn(1, numel(R.WF_MEAN)));
for m=1:numel(R.WF_MEASURES)
    switch R.WF_TRANSFORM{m}
        case 'log',      W.(R.WF_MEASURES{m}) = exp(T(:,m));
        case 'log(-x)',  W.(R.WF_MEASURES{m}) = -exp(T(:,m));
        otherwise,       W.(R.WF_MEASURES{m}) = T(:,m);
    end
end
t = linspace(-1, 1, numel(R.NOISE_SHAPE));
W.waveform = -exp(-(t / 0.2).^2) + 0.3 * sin(3 * t) + 0.02 * randn(n, numel(t));
W.nonSomatic = false(n, 1);
W.parameters.samplingRate = R.REF_SAMPLING_RATE;
Ft = R.FP_MEAN + R.FP_SD .* ((2 * randn(n, 1)) * D.b' + 0.5 * randn(n, 4));
FP.relativeRate = nan(n, 1);  FP.lv = Ft(:,2);  FP.cv = exp(Ft(:,3));  FP.burstIndex = Ft(:,4);
FP.rateOutside = 10.^(Ft(:,1) + randn);
FP.enoughSpikes = rand(n, 1) > 0.1;
FP.nSpikesOutside = 50 + 100 * FP.enoughSpikes;
FP.parameters.MIN_REFERENCE = 5;
FP = FiringPatternFeatures(FP, [], [], []);
NS = NeuronScores(FP, [], [], W, NaN);
end

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
