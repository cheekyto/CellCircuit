function passed = test_Analyze_cNE(code_dir)

% Function-level test of "Analyze_cNE", the cNE step of "Analyze_spikeTrain" (November 2026), on a synthetic recording with five
% conditions, in this order: a baseline (0-300 s), a "drug" (300-600 s) in which one of the two planted ensembles no longer fires
% together, a second drug (600-900 s, the same ensembles as under the first drug), a condition with too few units, and TTX. It checks that
%   1. each condition gets its cNEs and their description (the same as "cNEDetection" + "cNEAnalysis" called directly); each drug is
%      compared with the baseline and the second drug with the first one (the same as "cNECompare" detecting both conditions itself); the
%      2nd half of the baseline with its 1st half; the condition with too few units and its comparisons are left out, with their status;
%      the TTX condition is left out of every step ("excluded (TTX)") and is no one's previous condition ("previousCondition")
%   2. a rerun keeps the results of the conditions that were not analyzed again, and computes again the comparisons that involve the ones
%      that were (all of them when the baseline is analyzed again, or when the baseline is another condition)
%   3. a step that fails (a reversed window) is reported with its status and a warning; the other steps go on
%   4. the figures are saved (one per condition, per comparison and for the baseline halves; the index of the condition is added when two
%      conditions have the same name), hidden figures are closed, and the default visibility of the figures is restored
%
% Example (from a terminal):
%   matlab -batch "addpath('<repo>/tests'); test_Analyze_cNE('<repo>')"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dir:             folder with the code to be tested

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% passed:               =true if all the checks passed (each check is printed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(char(code_dir));
passed = true;
N = 30;     BASE = [0 300];     DRUG = [300 600];   DRUG2 = [600 900];  OTHER = [900 1200];     TTX = [1200 1500];
G = {1:6, 10:15};                                                           % planted ensembles
FAST = struct('detection', struct('N_RESAMPLE', 0, 'N_CHANCE', 0), 'analysis', struct('N_PERM', 200));   % fast settings for the test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Synthetic recording

rng(2);
trains = background(N, TTX(2));
trains = plant(trains, G{1}, events(0.3, [0 DRUG2(2)]), 0.7);              % kept under both drugs
trains = plant(trains, G{2}, events(0.3, BASE), 0.7);                       % coordinated in the baseline only ...
for i = G{2}                                                                % ... same rate under the drugs, independent spikes
    trains{i} = [trains{i}; DRUG(1)*1000 + rand(round(0.7*0.3*(DRUG2(2) - DRUG(1))),1)*(DRUG2(2) - DRUG(1))*1000];
end
trains = plant(trains, G{1}, events(0.3, TTX), 0.7);                        % (an ensemble under TTX: left out anyway)
trains = cellfun(@sort, trains, 'UniformOutput', false);
ev = @(w) (w(1) + 10:20:w(2) - 10)';                                        % large network events (1 s, tails of 1 s)
P = struct('spikeTrain', {}, 'window', {}, 'NetworkBurstDetected', {}, 'condition', {}, 'dose', {}, 'too_few_units', {});
P(1) = condition(trains, BASE, network(ev(BASE)), "Baseline", "", false);
P(2) = condition(trains, DRUG, network(ev(DRUG)), "+Drug", "10uM", false);
P(3) = condition(trains, DRUG2, network(ev(DRUG2)), "+Drug2", "", false);
P(4) = condition(trains(1:3), OTHER, [], "+Other", "", true);
P(5) = condition(trains, TTX, network(ev(TTX)), "+TTX", "2uM", false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Each condition, against the baseline and the condition before, the baseline halves

tic;
R = Analyze_cNE(P, 1, 1:5, [], FAST);
fprintf("  (first run: %.0f sec)\n", toc);
passed = report(sprintf("status of the conditions: %s", strjoin([R.condition.status], " | ")), ...
    isequal([R.condition.status], ["ok" "ok" "ok" "too few units" "excluded (TTX)"]) && ...
    isequal([R.condition.name], ["Baseline" "+Drug10uM" "+Drug2" "+Other" "+TTX2uM"]) && isempty(R.condition(5).cNE)) && passed;
passed = report(sprintf("against the baseline: %s | halves: %s", strjoin([R.comparison.status], " | "), R.halves.status), ...
    isequal([R.comparison.status], ["not compared" "ok" "ok" "no cNEs for the condition (too few units)" "excluded (TTX)"]) && ...
    R.halves.status == "ok" && R.baseline == 1) && passed;
passed = report(sprintf("against the condition before: %s (condition before: %s)", strjoin([R.previous.status], " | "), ...
    mat2str([R.previous.reference])), isequal([R.previous.reference], [0 0 2 3 0]) && isequal([R.previous.status], ...
    ["not compared" "not compared" "ok" "no cNEs for the condition (too few units)" "not compared"])) && passed;
passed = report("previousCondition: in the order of the windows, TTX left out, not the baseline", ...
    isequal(previousCondition({P.window}, 1, [false false false false true]), [0 0 2 3 0]) && ...
    isequal(previousCondition({P([3 1 2]).window}, 2), [3 0 0]) && isequal(previousCondition([0 10; NaN NaN; 20 30], 1), [0 0 0])) && passed;
C2 = cNEDetection(P(2).spikeTrain, DRUG, FAST.detection);
A2 = cNEAnalysis(C2, P(2).NetworkBurstDetected, P(2).spikeTrain, FAST.analysis);
passed = report(sprintf("drug condition: the same as cNEDetection + cNEAnalysis (%d cNEs)", C2.nCNE), ...
    isequaln(R.condition(2).cNE, C2) && isequaln(R.condition(2).report.summary, A2.summary)) && passed;
in = struct('DETECTION', FAST.detection, 'ANALYSIS', FAST.analysis);
K = cNECompare(input_of(P(1), BASE, "Baseline"), input_of(P(2), DRUG, "+Drug10uM"), in);
S = R.comparison(2).result.summary;
passed = report(sprintf("drug vs baseline: the same as cNECompare (%d cNEs -> %d; found %d, lost %d, new %d)", S.nCNE, S.nCNETest, ...
    S.nFound, S.nLost, S.nNew), isequaln(S, K.summary) && S.nCNE == 2 && S.nFound == 1 && S.nLost == 1) && passed;
K3 = cNECompare(input_of(P(2), DRUG, "+Drug10uM"), input_of(P(3), DRUG2, "+Drug2"), in);
S3 = R.previous(3).result.summary;
passed = report(sprintf("drug 2 vs drug 1: the same as cNECompare (%d cNEs -> %d; found %d, lost %d)", S3.nCNE, S3.nCNETest, S3.nFound, ...
    S3.nLost), isequaln(S3, K3.summary) && R.previous(3).name == "+Drug10uM vs +Drug2" && S3.nFound == S3.nCNE) && passed;
H = cNECompare(input_of(P(1), [0 150], "Baseline 1st half"), input_of(P(1), [150 300], "Baseline 2nd half"), in);
passed = report(sprintf("baseline halves: the same as cNECompare (found %d of %d)", R.halves.result.summary.nFound, ...
    R.halves.result.summary.nCNE), isequaln(R.halves.result.summary, H.summary) && R.halves.result.summary.nFound == 2) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Reruns

M = R;                                                                      % markers: kept results keep them, recomputed ones lose them
M.condition(1).report.summary.nCNE = 99;    M.condition(2).report.summary.nCNE = 99;    M.condition(3).report.summary.nCNE = 99;
M.comparison(2).result.summary.nFound = 99;  M.comparison(3).result.summary.nFound = 99;
M.previous(3).result.summary.nFound = 99;    M.halves.result.summary.nFound = 99;
M.condition(5) = R.condition(2);                                            % (TTX with cNEs, as if from an earlier run: removed)
R2 = Analyze_cNE(P, 1, 2, M, FAST);                                         % only the first drug analyzed again
passed = report("rerun of drug 1: baseline, drug 2, its baseline comparison and the halves kept; drug 1 and its comparisons computed again", ...
    R2.condition(1).report.summary.nCNE == 99 && R2.condition(3).report.summary.nCNE == 99 && ...
    R2.comparison(3).result.summary.nFound == 99 && R2.halves.result.summary.nFound == 99 && ...
    isequaln(R2.condition(2).report.summary, A2.summary) && isequaln(R2.comparison(2).result.summary, K.summary) && ...
    R2.previous(3).result.summary.nFound ~= 99 && R2.condition(5).status == "excluded (TTX)" && isempty(R2.condition(5).cNE)) && passed;
R1 = Analyze_cNE(P, 1, 1, M, FAST);                                         % only the baseline analyzed again
passed = report("rerun of the baseline: drugs and drug 2 vs drug 1 kept; their baseline comparisons and the halves computed again", ...
    R1.condition(2).report.summary.nCNE == 99 && R1.previous(3).result.summary.nFound == 99 && ...
    isequaln(R1.comparison(2).result.summary, K.summary) && R1.comparison(3).result.summary.nFound ~= 99 && ...
    isequaln(R1.halves.result.summary, H.summary)) && passed;
Rb = Analyze_cNE(P, 2, [], M, FAST);                                        % another baseline: everything computed again
passed = report(sprintf("another baseline: everything computed again (against it: %s; before: %s)", strjoin([Rb.comparison.status], " | "), ...
    mat2str([Rb.previous.reference])), Rb.baseline == 2 && isequaln(Rb.condition(2).report.summary, A2.summary) && ...
    Rb.condition(1).report.summary.nCNE ~= 99 && Rb.comparison(1).status == "ok" && Rb.comparison(2).status == "not compared" && ...
    Rb.halves.result.summary.nFound ~= 99 && isequal([Rb.previous.reference], [0 0 0 3 0])) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. A step that fails

Pbad = P;   Pbad(2).window = [450 300];                                     % a reversed window: cNEDetection stops with an error
lastwarn('');
warning('off', 'Analyze_cNE:error');
Re = Analyze_cNE(Pbad, 1, 1:5, [], FAST);
warning('on', 'Analyze_cNE:error');
[~, wid] = lastwarn;
short = @(s) extractBefore(s + blanks(40), 41);                             % first 40 characters of a status
passed = report(sprintf('failed step: drug 1 "%s...", its comparisons "%s...", "%s..."; the rest done', short(Re.condition(2).status), ...
    short(Re.comparison(2).status), short(Re.previous(3).status)), startsWith(Re.condition(2).status, "error: ") && ...
    startsWith(Re.comparison(2).status, "no cNEs for the condition (error: ") && ...
    startsWith(Re.previous(3).status, "no cNEs for the previous condition (error: ") && Re.condition(1).status == "ok" && ...
    Re.comparison(3).status == "ok" && Re.halves.status == "ok" && strcmp(wid, 'Analyze_cNE:error')) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Figures

folder = tempname;  mkdir(folder);
cleanup = onCleanup(@() rmdir(folder, 's'));
close all;
before = get(groot, 'DefaultFigureVisible');
plots = FAST;   plots.plots = true;     plots.hide = true;  plots.folder = folder;  plots.file_data = "rec";
Analyze_cNE(P, 1, 1:5, [], plots);
d = dir(fullfile(folder, '*.fig'));     files = string({d.name});
expected = ["rec_cNE_Baseline.fig" "rec_cNE_+Drug10uM.fig" "rec_cNE_+Drug2.fig" "rec_cNECompare_+Drug10uM.fig" ...
    "rec_cNECompare_+Drug2.fig" "rec_cNECompare_+Drug2_vs_+Drug10uM.fig" "rec_cNECompare_BaselineHalves.fig"];
passed = report(sprintf("figures saved (%d, none for TTX); none left open; default visibility restored", numel(files)), ...
    isempty(setxor(files, expected)) && isempty(findall(groot, 'Type', 'figure')) && isequal(get(groot, 'DefaultFigureVisible'), before)) && passed;
P2 = P(1:2);    P2(3) = P(2);                                               % two conditions with the same name (and window)
delete(fullfile(folder, '*.fig'));
Analyze_cNE(P2, 1, 1:3, [], plots);
d = dir(fullfile(folder, '*.fig'));     files = string({d.name});
expected = ["rec_cNE_Baseline.fig" "rec_cNE_+Drug10uM_2.fig" "rec_cNE_+Drug10uM_3.fig" "rec_cNECompare_+Drug10uM_2.fig" ...
    "rec_cNECompare_+Drug10uM_3.fig" "rec_cNECompare_BaselineHalves.fig"];
passed = report(sprintf("two conditions with the same name: their figures get the index of the condition (%d figures)", numel(files)), ...
    isempty(setxor(files, expected))) && passed;
f = openfig(fullfile(folder, "rec_cNE_Baseline.fig"));
passed = report("a saved hidden figure opens visible", strcmp(f.Visible, 'on')) && passed;
close(f);

if passed
    fprintf("Analyze_cNE: all checks passed\n");
else
    fprintf("Analyze_cNE: SOME CHECKS FAILED\n");
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = report(name, ok)
if ok
    fprintf("  PASS  %s\n", name);
else
    fprintf("  FAIL  %s\n", name);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = condition(trains, window, NB, name, dose, too_few)
% one condition, as "Analyze_spikeTrain" saves it in PerCondition (spike trains: structure array with .timeStamp, in msec)
C = struct('spikeTrain', struct('timeStamp', trains(:), 'cluster', num2cell((1:numel(trains))')), 'window', window, ...
    'NetworkBurstDetected', NB, 'condition', name, 'dose', dose, 'too_few_units', too_few);
end

function S = input_of(C, window, name)
% one condition, as "cNECompare" takes it (cNEs detected by "cNECompare")
S = struct('spikeTimestamp', {C.spikeTrain}, 'tRange', window, 'NetworkBurstDetected', {C.NetworkBurstDetected}, 'name', name);
end

function trains = background(n_units, T_sec)
% independent Poisson units (1-4 Hz), spike times in msec
rates = 1 + 3*rand(n_units,1);
trains = cell(n_units,1);
for i=1:n_units
    trains{i} = rand(poissrnd(rates(i)*T_sec),1)*T_sec*1000;
end
end

function t = events(rate, range)
% Poisson event times (sec) in the range
t = sort(range(1) + 0.01 + rand(poissrnd(rate*diff(range)),1)*(diff(range) - 0.02));
end

function trains = plant(trains, group, t, p)
% at each event time t (sec), each member fires with probability p, within +/- 2 msec
for i=group
    part = rand(numel(t),1) < p;
    trains{i} = [trains{i}; t(part)*1000 + (2*rand(sum(part),1) - 1)*2];
end
end

function NB = network(s)
% output of "NetworkBurstDetection" (option 5) with large events [s, s + 1] and tails of 1 s
E = struct('start', num2cell(s), 'last', num2cell(s + 1), 'peak', num2cell(s + 0.3), 'level', num2cell(ones(numel(s),1)));
NB = struct('nBurst', numel(E), 'IndivBurst', E, 'closedZones', [s, s + 1, s + 2]);
end
