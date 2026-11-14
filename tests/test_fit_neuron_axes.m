function passed = test_fit_neuron_axes(code_dir)

% Function-level test of "FitNeuronAxes" and AXES_FILE of "NeuronScores" on synthetic organoids (new, November 2026)
% It checks that
%   1. "FiringPatternFeatures" given a firing pattern already measured computes its relative rate again: the same values against the same
%      reference units, the values of a full measurement against other reference units
%   2. the axes fitted on organoids of 3 groups are those that generated the data (|cosine| > 0.95), with their signs; placed on them with
%      AXES_FILE, the units used have scores of weighted mean 0, weighted SD 1 and no weighted correlation (every organoid weighted equally)
%   3. every organoid counts the same (its units given twice: the same axes); look-alikes and non-somatic units are left out (extreme values
%      added for them: the same axes); eigenvalues and stability between halves of the organoids are reported
%   4. nothing is fitted with too few organoids, a single group, a file given twice, a file without scores, templates of different kinds,
%      an existing axes file or the version of the axes already used
%   5. AXES_FILE: the axes of section REFERENCE written in a file give the same scores; a missing file or wrong axes stop with an error
%   6. templates of another sampling rate: no waveform scores on the fixed axes but their measures are kept; axes fitted on them (no QC noise
%      shape) give them waveform scores, and the relative rate is then measured against their somatic units
% No reference recording is needed
%
% Example (from a terminal):
%   matlab -batch "addpath('<repo>/tests'); test_fit_neuron_axes('<repo>')"

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
TOL = 1e-9;
tmp = tempname;
mkdir(tmp);
cleanup = onCleanup(@() cleanup_test(tmp));

% the fixed axes, and the directions that generate the synthetic data
R = NeuronScores(fp_struct(zeros(1, 4)), [], [], [], NaN).reference;
D.a1 = R.WF_AXES(:,1) / norm(R.WF_AXES(:,1));
a2 = R.WF_AXES(:,2) - (D.a1' * R.WF_AXES(:,2)) * D.a1;
D.a2 = a2 / norm(a2);
D.b = R.FP_AXIS / norm(R.FP_AXIS);
D.samples = numel(R.NOISE_SHAPE);   D.rate = R.REF_SAMPLING_RATE;
groups = repelem(["A"; "B"; "C"], 8);
NS = arrayfun(@(k) synth(R, D, 20 + randi(100)), (1:24)', 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Relative rate of a firing pattern already measured

n = 12;
for i=n:-1:1
    trains(i).timeStamp = sort(rand(1, 200 + 40 * i) * 600e3);              % msec
end
somatic = [true(9, 1); false(3, 1)];
FPall = FiringPatternFeatures(trains, [0 600], [], []);
FPsom = FiringPatternFeatures(trains, [0 600], [], somatic);
again = FiringPatternFeatures(FPall, [], [], []);
other = FiringPatternFeatures(FPall, [], [], somatic);
passed = report("a firing pattern already measured: the same relative rate against the same units, the values of a full measurement " + ...
    "against other units", isequaln(again.relativeRate, FPall.relativeRate) && max(abs(other.relativeRate - FPsom.relativeRate)) < 1e-12 && ...
    isequaln(rmfield(other, 'relativeRate'), rmfield(FPall, 'relativeRate'))) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Axes fitted on 24 organoids of 3 groups

[A, out] = evalc_fit(NS, groups, fullfile(tmp, 'axes.mat'), "test axes");
cos_wf = abs(sum(A.WF_AXES ./ vecnorm(A.WF_AXES) .* [D.a1 D.a2]));
cos_fp = abs(A.FP_AXIS' * D.b) / norm(A.FP_AXIS);
passed = report(sprintf("the axes that generated the data are found (|cos| %.3f, %.3f and %.3f)", cos_wf, cos_fp), ...
    all(cos_wf > 0.95) && cos_fp > 0.95 && contains(out, "test axes") && isfile(fullfile(tmp, 'axes.mat'))) && passed;

S = cellfun(@(x) score(x, fullfile(tmp, 'axes.mat')), NS, 'UniformOutput', false);
[m1, s1, r12] = wstats(S, 'waveformScore1', 'waveformScore2', 'waveformScored');
[m2, s2] = wstats(S, 'waveformScore2', 'waveformScore1', 'waveformScored');
[m3, s3] = wstats(S, 'firingPatternScore', 'firingPatternScore', 'firingScored');
passed = report("with AXES_FILE, the units used get scores of weighted mean 0, weighted SD 1, and the two waveform scores are uncorrelated", ...
    max(abs([m1 m2 m3 r12])) < TOL && max(abs([s1 s2 s3] - 1)) < TOL && S{1}.reference.version == "test axes" && ...
    S{1}.reference.file == string(fullfile(tmp, 'axes.mat'))) && passed;

hw = cellfun(@(x) x.WaveformFeature.halfWidth(x.waveformScored), S, 'UniformOutput', false);
as = cellfun(@(x) x.WaveformFeature.asymmetry(x.waveformScored), S, 'UniformOutput', false);
cv = cellfun(@(x) x.FiringPattern.cv(x.firingScored), S, 'UniformOutput', false);
w1 = cellfun(@(x) x.waveformScore1(x.waveformScored), S, 'UniformOutput', false);
w2 = cellfun(@(x) x.waveformScore2(x.waveformScored), S, 'UniformOutput', false);
f1 = cellfun(@(x) x.firingPatternScore(x.firingScored), S, 'UniformOutput', false);
passed = report("signs: waveform score 1 grows with the half-width, score 2 with the asymmetry, the firing-pattern score with the CV", ...
    pearson(vertcat(w1{:}), vertcat(hw{:})) > 0 && pearson(vertcat(w2{:}), vertcat(as{:})) > 0 && pearson(vertcat(f1{:}), vertcat(cv{:})) > 0) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Equal weight per organoid, units left out, report

NS2 = NS;
NS2{5} = twice(NS{5});
A2 = evalc_fit(NS2, groups, "", "test axes");
passed = report("every organoid counts the same: the units of one organoid given twice leave the axes unchanged", ...
    max(abs([A2.WF_AXES(:) - A.WF_AXES(:); A2.FP_AXIS - A.FP_AXIS; A2.WF_MEAN(:) - A.WF_MEAN(:)])) < 1e-10) && passed;

NS3 = NS;
NS3{3} = with_outliers(NS{3}, R);
A3 = evalc_fit(NS3, groups, "", "test axes");
passed = report("look-alikes of the QC noise shape and non-somatic units (with extreme measures) are left out of the axes", ...
    max(abs([A3.WF_AXES(:) - A.WF_AXES(:); A3.FP_AXIS - A.FP_AXIS; A3.WF_SD(:) - A.WF_SD(:)])) < 1e-10 && ...
    sum(NS3{3}.lookAlike) == 5) && passed;

F = A.fit;
passed = report(sprintf("eigenvalues and stability reported (halves of the organoids: median |cos| %.3f, %.3f, %.3f)", ...
    median(F.stabilityWaveform), median(F.stabilityFiring)), abs(sum(F.eigenvaluesWaveform) - 10) < 1e-9 && ...
    abs(sum(F.eigenvaluesFiring) - 4) < 1e-9 && issorted(F.eigenvaluesWaveform, 'descend') && ...
    isequal(A.WF_EIGENVALUES, F.eigenvaluesWaveform(1:2)) && isequal(size(F.stabilityWaveform), [50 2]) && ...
    all(F.stabilityWaveform(:) <= 1 + 1e-12) && all(median(F.stabilityWaveform) > 0.9) && median(F.stabilityFiring) > 0.9 && ...
    F.nOrganoidsWaveform == 24 && F.nGroups == 3 && height(F.organoids) == 24) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Refusals

files = strings(4, 1);
for k=1:4
    files(k) = fullfile(tmp, sprintf('org%d_curated_Analyzed.mat', k));
    NeuronScore = NS{k};  %#ok<NASGU>
    save(files(k), 'NeuronScore');
end
no_scores = fullfile(tmp, 'other_curated_Analyzed.mat');
x = 1;  save(no_scores, 'x');
NS_rate = NS;
NS_rate{2}.WaveformFeature.parameters.samplingRate = 2 * D.rate;
refused = [fails(@() evalc_fit(NS(1:19), groups(1:19), "", ""), "at least 20"), ...                % 19 organoids
           fails(@() evalc_fit(NS, repmat("A", 24, 1), "", ""), "1 group(s)"), ...                  % a single group
           fails(@() evalc_fit(files([1 2 2]), ["A"; "B"; "B"], "", ""), "given twice"), ...        % a file given twice
           fails(@() evalc_fit([files(1:2); no_scores], ["A"; "B"; "B"], "", ""), "no neuron scores"), ...   % a file without scores
           fails(@() evalc_fit(NS_rate, groups, "", ""), "templates of the sources differ"), ...    % templates of different kinds
           fails(@() evalc_fit(NS, groups, fullfile(tmp, 'axes.mat'), ""), "already exists"), ...   % the axes file exists
           fails(@() evalc_fit(NS, groups, "", R.version), "already used")];                        % the version of the fixed axes
passed = report(sprintf("nothing fitted with too few organoids, one group, a file twice, a file without scores, templates of different " + ...
    "kinds, an existing file or a used version (%d of 7 refused)", sum(refused)), all(refused)) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. AXES_FILE with the fixed axes; wrong files

NeuronAxes = rmfield(R, 'file');
NeuronAxes.version = "copy of the fixed axes";
save(fullfile(tmp, 'fixed.mat'), 'NeuronAxes');
S0 = score(NS{1}, "");
S1 = score(NS{1}, fullfile(tmp, 'fixed.mat'));
passed = report("the fixed axes written in a file give the same scores", isequaln([S0.waveformScore1 S0.waveformScore2 S0.firingPatternScore ...
    S0.noiseSimilarity S0.weightedScores], [S1.waveformScore1 S1.waveformScore2 S1.firingPatternScore S1.noiseSimilarity S1.weightedScores]) && ...
    S1.reference.version == "copy of the fixed axes" && S0.reference.version == R.version) && passed;

bad_size = NeuronAxes;  bad_size.WF_AXES = bad_size.WF_AXES(:, 1);
bad_measures = NeuronAxes;  bad_measures.WF_MEASURES{1} = 'width';
save(fullfile(tmp, 'bad_size.mat'), '-struct', 'bad_size');                 % (not in a variable NeuronAxes)
NeuronAxes = bad_size;      save(fullfile(tmp, 'bad_size2.mat'), 'NeuronAxes');
NeuronAxes = bad_measures;  save(fullfile(tmp, 'bad_measures.mat'), 'NeuronAxes');
refused = [fails(@() score(NS{1}, fullfile(tmp, 'missing.mat')), "does not exist"), ...
    fails(@() score(NS{1}, fullfile(tmp, 'bad_size.mat')), "has no axes"), ...
    fails(@() score(NS{1}, fullfile(tmp, 'bad_size2.mat')), "expected sizes"), ...
    fails(@() score(NS{1}, fullfile(tmp, 'bad_measures.mat')), "other measures")];
passed = report("AXES_FILE: a missing file, a file without NeuronAxes, axes of the wrong size or of other measures stop with an error", ...
    all(refused)) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Templates of another sampling rate

warnings = warning('off', 'all');                                           % (the templates are not like the reference ones)
D2 = D;  D2.samples = 50;  D2.rate = 17855;
NSo = arrayfun(@(k) synth(R, D2, 20 + randi(100), true), (1:24)', 'UniformOutput', false);
kept = all(cellfun(@(x) ~x.templatesUsed && ~isempty(x.WaveformFeature) && all(isnan(x.waveformScore1)) && any(x.firingScored), NSo));
passed = report("templates of another sampling rate: no waveform scores on the fixed axes, their measures kept", kept) && passed;

Ao = evalc_fit(NSo, groups, fullfile(tmp, 'axes_other.mat'), "");
So = cellfun(@(x) score(x, fullfile(tmp, 'axes_other.mat')), NSo, 'UniformOutput', false);
[m1, s1] = wstats(So, 'waveformScore1', 'waveformScore1', 'waveformScored');
som = ~NSo{1}.WaveformFeature.nonSomatic;
rr = FiringPatternFeatures(NSo{1}.FiringPattern, [], [], som).relativeRate;
passed = report(sprintf("axes fitted on them (%d samples at %.0f Hz, no QC noise shape): waveform scores (weighted mean %.1g, SD %.6f), " + ...
    "no look-alikes, relative rate against their somatic units", Ao.REF_SAMPLES, Ao.REF_SAMPLING_RATE, m1, s1), ...
    Ao.REF_SAMPLES == 50 && abs(Ao.REF_SAMPLING_RATE - 17855) < 1 && isempty(Ao.NOISE_SHAPE) && abs(m1) < TOL && abs(s1 - 1) < TOL && ...
    all(cellfun(@(x) x.templatesUsed && ~any(x.lookAlike) && all(isnan(x.noiseSimilarity)), So)) && ...
    isequaln(So{1}.FiringPattern.relativeRate, rr) && ~isequaln(rr, NSo{1}.FiringPattern.relativeRate) && ...
    contains(Ao.version, "24 organoids of 3 groups")) && passed;
warning(warnings);

if passed
    fprintf("test_fit_neuron_axes: all checks passed\n");
else
    fprintf("test_fit_neuron_axes: SOME CHECKS FAILED\n");
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NS = synth(R, D, n, nonsomatic)
% scores of a synthetic organoid: 10 waveform measures from two latent variables along D.a1 and D.a2 (SD 2 and 1.5, noise 0.3), and 4
% firing-pattern features from one along D.b (SD 2, noise 0.5), with an offset of the organoid; with nonsomatic, 20% of the units are
% non-somatic with another firing rate
if nargin < 4
    nonsomatic = false;
end
shift = 0.2 * randn(1, numel(R.WF_MEAN));
T = R.WF_MEAN + R.WF_SD .* ((randn(n, 2) .* [2 1.5]) * [D.a1 D.a2]' + 0.3 * randn(n, numel(R.WF_MEAN)) + shift);
for m=1:numel(R.WF_MEASURES)
    switch R.WF_TRANSFORM{m}
        case 'log',      W.(R.WF_MEASURES{m}) = exp(T(:,m));
        case 'log(-x)',  W.(R.WF_MEASURES{m}) = -exp(T(:,m));
        otherwise,       W.(R.WF_MEASURES{m}) = T(:,m);
    end
end
t = linspace(-1, 1, D.samples);
W.waveform = -exp(-(t / 0.2).^2) + 0.3 * sin(3 * t) + 0.02 * randn(n, D.samples);
W.nonSomatic = false(n, 1);
if nonsomatic
    W.nonSomatic(1:round(n / 5)) = true;
end
W.parameters.samplingRate = D.rate;
Fz = (2 * randn(n, 1)) * D.b' + 0.5 * randn(n, 4);
Ft = R.FP_MEAN + R.FP_SD .* Fz;
FP.relativeRate = nan(n, 1);
FP.lv = Ft(:,2);  FP.cv = exp(Ft(:,3));  FP.burstIndex = Ft(:,4);
FP.rateOutside = 10.^(Ft(:,1) + randn + W.nonSomatic);
FP.enoughSpikes = rand(n, 1) > 0.2;
FP.lv(~FP.enoughSpikes) = NaN;  FP.cv(~FP.enoughSpikes) = NaN;  FP.burstIndex(~FP.enoughSpikes) = NaN;
FP.nSpikesOutside = 50 + 100 * FP.enoughSpikes;
FP.parameters.MIN_REFERENCE = 5;
FP = FiringPatternFeatures(FP, [], [], []);
NS = NeuronScores(FP, [], [], W, NaN);
end

function NS = twice(NS)
% the same organoid with each unit given twice
W = NS.WaveformFeature;  FP = NS.FiringPattern;
for f=string(fieldnames(W))'
    if f ~= "parameters", W.(f) = [W.(f); W.(f)]; end
end
for f=["relativeRate", "lv", "cv", "burstIndex", "rateOutside", "enoughSpikes", "nSpikesOutside"]
    FP.(f) = [FP.(f); FP.(f)];
end
NS = NeuronScores(FP, [], [], W, NaN);
end

function NS = with_outliers(NS, R)
% 5 look-alikes (template = QC noise shape) and 5 non-somatic units added, all with extreme measures and features
W = NS.WaveformFeature;  FP = NS.FiringPattern;
for m=1:numel(R.WF_MEASURES)
    W.(R.WF_MEASURES{m}) = [W.(R.WF_MEASURES{m}); repmat(10 * W.(R.WF_MEASURES{m})(1), 10, 1)];
end
W.waveform = [W.waveform; repmat(R.NOISE_SHAPE, 5, 1); W.waveform(1:5,:)];
W.nonSomatic = [W.nonSomatic; false(5, 1); true(5, 1)];
FP.lv = [FP.lv; 50 * ones(10, 1)];  FP.cv = [FP.cv; 50 * ones(10, 1)];  FP.burstIndex = [FP.burstIndex; 50 * ones(10, 1)];
FP.rateOutside = [FP.rateOutside; 1e3 * ones(10, 1)];
FP.enoughSpikes = [FP.enoughSpikes; true(10, 1)];  FP.nSpikesOutside = [FP.nSpikesOutside; 500 * ones(10, 1)];
FP.relativeRate = [FP.relativeRate; nan(10, 1)];
NS = NeuronScores(FP, [], [], W, NaN);
end

function S = score(NS, axes_file)
% scores of the saved measures and features with AXES_FILE
parameter_set('set', struct('NeuronScores', struct('AXES_FILE', string(axes_file))));
c = onCleanup(@() parameter_set('clear'));
S = NeuronScores(NS.FiringPattern, [], [], NS.WaveformFeature, NaN);
end

function [m, s, r] = wstats(S, name, other, scored)
% weighted mean and SD of a score over the scored units (every organoid weighted equally), and its weighted correlation with another one
x = [];  y = [];  w = [];
for k=1:numel(S)
    u = S{k}.(scored);
    x = [x; S{k}.(name)(u)];  y = [y; S{k}.(other)(u)];  w = [w; ones(sum(u), 1) / sum(u)];   %#ok<AGROW>
end
w = w / sum(w);
m = w' * x;
s = sqrt(w' * (x - m).^2);
my = w' * y;
r = (w' * ((x - m) .* (y - my))) / (s * sqrt(w' * (y - my).^2));
end

function [A, out] = evalc_fit(sources, groups, file, version)
% FitNeuronAxes, its printed report in out
[out, A] = evalc('FitNeuronAxes(sources, groups, file, version)');
end

function F = fp_struct(V)
% firing pattern with the features V (one row per unit; transformed values)
F.relativeRate = V(:,1);  F.lv = V(:,2);  F.cv = exp(V(:,3));  F.burstIndex = V(:,4);
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

function r = pearson(a, b)
% Pearson correlation
c = corrcoef(a, b);
r = c(1, 2);
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
