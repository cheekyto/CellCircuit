function passed = test_neuron_scores(code_dir)

% Function-level test of "NeuronScores" on synthetic inputs (new, November 2026)
% Like "test_firing_pattern_features", it does not compare two versions of the code: it checks that
%   1. the scores are the projections of the standardized measures on the reference axes: 0 at the reference means, one SD of a measure
%      moves each score by its loading, and two real units of the reference set get the scores they had when the axes were computed
%      (and the axes carry their version)
%   2. the similarity to the QC noise shape is a correlation (1 for the shape itself, also scaled and shifted; -1 for the flipped shape); the
%      look-alikes, the non-somatic units and the units with a missing measure get no scores; "clean" = similarity <= CLEAN_SIMILARITY
%   3. the firing-pattern score: the reference of the relative rate is the somatic units without look-alikes (all units without templates)
%   4. templates given as a matrix = their "WaveformFeatures" given directly; templates of another length or sampling rate are not used
%   5. weightedScores follow WAVEFORM_WEIGHT (also when changed through a parameter set) and the scores do not depend on it
%   6. composition: means of the scored units, NaN below MIN_COMPOSITION_UNITS; edge cases (no templates, wrong number of units)
% No reference recording is needed
%
% Example (from a terminal):
%   matlab -batch "addpath('<repo>/tests'); test_neuron_scores('<repo>')"

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
rng(20261111);
TOL = 1e-9;

% the reference axes, read from the output of a call
R = NeuronScores(fp_struct(zeros(1, 4)), [], [], [], NaN).reference;
n_wf = numel(R.WF_MEASURES);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Projection on the reference axes

T = [R.WF_MEAN; R.WF_MEAN + R.WF_SD .* eye(n_wf)];                         % reference means, then +1 SD of each measure
W = wf_struct(T, R);
F = [R.FP_MEAN; R.FP_MEAN + R.FP_SD .* eye(4)];
F = [F; repmat(R.FP_MEAN, n_wf + 1 - size(F, 1), 1)];
NS = NeuronScores(fp_struct(F), [], [], W, NaN);
passed = report("scores = 0 at the reference means", abs(NS.waveformScore1(1)) < TOL && abs(NS.waveformScore2(1)) < TOL && ...
    abs(NS.firingPatternScore(1)) < TOL) && passed;
passed = report("+1 SD of each waveform measure moves the scores by its loadings", ...
    max(abs([NS.waveformScore1(2:end) NS.waveformScore2(2:end)] - R.WF_AXES), [], 'all') < TOL) && passed;
passed = report("+1 SD of each firing-pattern feature moves the score by its loading", max(abs(NS.firingPatternScore(2:5) - R.FP_AXIS)) < TOL) && passed;

% two units of the reference set (measures and features as computed in November 2026) and their scores on the axes of 13 November 2026
% (equal weight per organoid; computed independently of this code)
measures = [0.977028422079 0.268011865187 0.300146697515 0.0657932520832 0.640415034459 -6.46296038611 4.31089283819 1.49921620154 0.212617584079 31.9272512133;
            0.632790428807 0.237392269663 0.346183362334 0.0048256498399 0.972504125692 -6.21159302809 5.67721766044 1.09412627798 0.202492937218 35.2345023797];
features = [0.321035930291 0.793201111263 1.90340978473 -1.00259798072;
            2.14334356218 0.333278900714 4.68028619508 -1.59900923982];
expected = [-0.0429133029274 -1.16537448208 0.882995135945;
            -0.412549920607 0.925959581015 3.38484285773];
W = wf_struct(measures, R, true);
NS = NeuronScores(fp_struct(features, true), [], [], W, NaN);
passed = report(sprintf("two reference units get their scores on the axes (%s)", R.version), ...
    max(abs([NS.waveformScore1 NS.waveformScore2 NS.firingPatternScore] - expected), [], 'all') < 1e-8) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Similarity to the QC noise shape, units left out

T = repmat(R.WF_MEAN, 7, 1);
W = wf_struct(T, R);
shape = R.NOISE_SHAPE;
W.waveform(1,:) = shape;                                                    % the shape itself
W.waveform(2,:) = 0.3 * shape - 0.1;                                        % scaled and shifted
W.waveform(3,:) = -shape;                                                   % flipped
W.waveform(4,:) = 0.9 * shape + 0.6 * randn(size(shape)) * std(shape);      % noisy copy
W.nonSomatic(5) = true;                                                     % non-somatic
W.troughToPeak(6) = NaN;                                                    % a measure missing
NS = NeuronScores(fp_struct(repmat(R.FP_MEAN, 7, 1)), [], [], W, NaN);
r4 = corr(W.waveform(4,:)', shape');
passed = report("similarity: 1 for the shape, also scaled and shifted; -1 flipped; Pearson r for a noisy copy", ...
    abs(NS.noiseSimilarity(1) - 1) < TOL && abs(NS.noiseSimilarity(2) - 1) < TOL && abs(NS.noiseSimilarity(3) + 1) < TOL && ...
    abs(NS.noiseSimilarity(4) - r4) < TOL) && passed;
passed = report("look-alikes, non-somatic units and units with a missing measure get no scores", ...
    isequal(NS.lookAlike(1:3)', [true true false]) && all(isnan(NS.waveformScore1([1 2 5 6]))) && ...
    all(isnan(NS.firingPatternScore([1 2 5]))) && ~any(NS.waveformScored([1 2 5 6])) && all(NS.waveformScored([3 4 7])) && ...
    NS.nonSomatic(5)) && passed;
passed = report("clean = similarity <= CLEAN_SIMILARITY", isequal(NS.clean, NS.noiseSimilarity <= NS.parameters.CLEAN_SIMILARITY)) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Firing-pattern score from spike trains: reference of the relative rate

T0 = 300;   T1 = 900;
n = 12;
rates = [1 2 3 4 5 6 7 8 30 40 50 60];                                      % the last 4 units: very active look-alikes / non-somatic
trains = struct('timeStamp', cell(1, n));
for i=1:n
    trains(i).timeStamp = sort(T0*1000 + (T1-T0)*1000 * rand(round(rates(i)*(T1-T0)), 1));
end
W = wf_struct(repmat(R.WF_MEAN, n, 1), R);
W.waveform(9:10,:) = repmat(shape, 2, 1);                                   % look-alikes
W.nonSomatic(11:12) = true;                                                 % non-somatic
NS = NeuronScores(trains, [T0 T1], [], W, NaN);
FP = FiringPatternFeatures(trains, [T0 T1], [], [true(8,1); false(4,1)]);
passed = report("relative rate against the somatic units without look-alikes", ...
    max(abs(NS.FiringPattern.relativeRate - FP.relativeRate)) < TOL && all(isfinite(NS.firingPatternScore(1:8))) && ...
    all(isnan(NS.firingPatternScore(9:12)))) && passed;
NS0 = NeuronScores(trains, [T0 T1], [], [], NaN);
FP0 = FiringPatternFeatures(trains, [T0 T1], []);
passed = report("without templates: no waveform scores, relative rate against all units, every unit with the features scored", ...
    ~NS0.templatesUsed && all(isnan(NS0.waveformScore1)) && all(isnan(NS0.noiseSimilarity)) && ...
    max(abs(NS0.FiringPattern.relativeRate - FP0.relativeRate)) < TOL && all(NS0.firingScored) && ...
    isnan(NS0.composition.nSomatic) && isempty(NS0.WaveformFeature)) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Templates: matrix = "WaveformFeatures" given directly; other lengths or sampling rates not used

fs = R.REF_SAMPLING_RATE;
t = 1:61;
templates = zeros(n, 82);
for i=1:n
    templates(i, 22:end) = -exp(-(t - 21).^2 / (4 + i)) + (0.2 + 0.03*i) * exp(-(t - 30 - i/2).^2 / 40);
end
A = NeuronScores(trains, [T0 T1], [], templates, fs);
B = NeuronScores(trains, [T0 T1], [], WaveformFeatures(templates, fs), NaN);
passed = report("templates as a matrix = their WaveformFeatures", A.templatesUsed && ...
    isequaln([A.waveformScore1 A.waveformScore2 A.firingPatternScore], [B.waveformScore1 B.waveformScore2 B.firingPatternScore]) && ...
    any(isfinite(A.waveformScore1))) && passed;
warning('off', 'all');
C1 = NeuronScores(trains, [T0 T1], [], templates(:, 1:72), fs);            % 51 samples after the padding
C2 = NeuronScores(trains, [T0 T1], [], templates, 30000);                  % another sampling rate
warning('on', 'all');
passed = report("templates of another length or sampling rate: not used", ~C1.templatesUsed && ~C2.templatesUsed && ...
    all(isnan(C1.waveformScore1)) && all(isnan(C2.waveformScore2))) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Weighted scores

T = repmat(R.WF_MEAN, 20, 1) + randn(20, n_wf) .* R.WF_SD;
F = repmat(R.FP_MEAN, 20, 1) + randn(20, 4) .* R.FP_SD;
W = wf_struct(T, R);
NS = NeuronScores(fp_struct(F), [], [], W, NaN);
lam = R.WF_EIGENVALUES;
S3 = [NS.waveformScore1 NS.waveformScore2 NS.firingPatternScore];
passed = report("weightedScores with the default weight 2/3", ...
    max(abs(NS.weightedScores - S3 .* sqrt([4/3, 4/3*lam(2)/lam(1), 2/3])), [], 'all') < TOL) && passed;
parameter_set('set', struct('NeuronScores', struct('WAVEFORM_WEIGHT', 0.5)));
NSh = NeuronScores(fp_struct(F), [], [], W, NaN);
parameter_set('clear');
passed = report("WAVEFORM_WEIGHT = 0.5 through a parameter set: equal weights, same scores", NSh.parameters.WAVEFORM_WEIGHT == 0.5 && ...
    max(abs(NSh.weightedScores - S3 .* sqrt([1, lam(2)/lam(1), 1])), [], 'all') < TOL && ...
    isequaln(S3, [NSh.waveformScore1 NSh.waveformScore2 NSh.firingPatternScore])) && passed;
passed = report("with weight 2/3, one SD on waveform axis 1 counts twice one SD on the firing-pattern axis", ...
    abs(NS.weightedScores(1,1)^2 / S3(1,1)^2 / (NS.weightedScores(1,3)^2 / S3(1,3)^2) - 2) < 1e-12) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Composition and edge cases

c = NS.composition;
passed = report("composition: means of the 20 scored units", c.nWaveformScored == 20 && c.nFiringScored == 20 && ...
    abs(c.meanWaveformScore1 - mean(S3(:,1))) < TOL && abs(c.meanFiringPatternScore - mean(S3(:,3))) < TOL && ...
    abs(c.meanNoiseSimilarity - mean(NS.noiseSimilarity)) < TOL) && passed;
W9 = wf_struct(T(1:9,:), R);
NS9 = NeuronScores(fp_struct(F(1:9,:)), [], [], W9, NaN);
passed = report("composition: means NaN below MIN_COMPOSITION_UNITS (9 units)", NS9.composition.nWaveformScored == 9 && ...
    isnan(NS9.composition.meanWaveformScore1) && isnan(NS9.composition.meanFiringPatternScore)) && passed;
W.waveform(1:12,:) = repmat(shape, 12, 1) + 0.05;                           % similarity 1: 12 look-alikes, 8 units left
NSl = NeuronScores(fp_struct(F), [], [], W, NaN);
passed = report("composition: look-alikes counted and left out of the means", NSl.composition.nLookAlike == 12 && ...
    NSl.composition.nWaveformScored == 8 && isnan(NSl.composition.meanWaveformScore2)) && passed;
ok = false;
try
    NeuronScores(fp_struct(F(1:5,:)), [], [], W, NaN);
catch
    ok = true;
end
passed = report("error when the waveforms and the spike trains have different numbers of units", ok) && passed;
E = NeuronScores(struct('timeStamp', cell(1, 0)), [T0 T1], [], [], NaN);
passed = report("no units", isempty(E.waveformScore1) && E.composition.nUnits == 0 && isnan(E.composition.meanFiringPatternScore)) && passed;

if passed
    fprintf("\nAll checks passed\n");
else
    fprintf("\nSOME CHECKS FAILED\n");
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = wf_struct(V, R, raw)
% output of "WaveformFeatures" for units with the transformed measures V (one row per unit, in the order of R.WF_MEASURES), or with the raw
% measures when raw = true; the waveform of every unit is a smooth spike (similarity to the noise shape well below the cutoffs)
if nargin < 3 || ~raw
    for m=1:numel(R.WF_MEASURES)
        switch R.WF_TRANSFORM{m}
            case 'log',     V(:,m) = exp(V(:,m));
            case 'log(-x)', V(:,m) = -exp(V(:,m));
        end
    end
end
n = size(V, 1);
for m=1:numel(R.WF_MEASURES)
    W.(R.WF_MEASURES{m}) = V(:,m);
end
t = 1:numel(R.NOISE_SHAPE);
W.waveform = repmat(-exp(-(t - 21).^2 / 6) + 0.6 * exp(-(t - 26).^2 / 6), n, 1);
W.nonSomatic = false(n, 1);
W.parameters = struct('samplingRate', R.REF_SAMPLING_RATE);
end

function F = fp_struct(V, ~)
% output of "FiringPatternFeatures" with the features V (one row per unit: relative rate, LV, CV, burst index); with one input, V holds the
% transformed values (log of the CV), with two the raw ones
if nargin < 2
    V(:,3) = exp(V(:,3));
end
F.relativeRate = V(:,1);  F.lv = V(:,2);  F.cv = V(:,3);  F.burstIndex = V(:,4);
F.rateOutside = nan(size(V, 1), 1);  F.nSpikesOutside = nan(size(V, 1), 1);
end

function ok = report(name, ok)
if ok
    fprintf("  PASS  %s\n", name);
else
    fprintf("  FAIL  %s\n", name);
end
end
