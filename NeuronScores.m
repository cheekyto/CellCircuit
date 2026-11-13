function NeuronScore = NeuronScores(spikeTrain, tRange, netBursts, waveforms, samplingRate, unitType)

% function to describe each neuron of a recording by three continuous scores, from its spike waveform and its firing pattern at baseline
% (new, November 2026)
% The neurons of the organoids do not form classes: the 10 core waveform measures and the 4 firing-pattern features together form one
% continuous cloud, no more grouped than cluster-free data of the same shape, with two independent axes, one for the waveform and one for
% the firing pattern (5,163 units of 51 organoids, November 2026). Each neuron is therefore described by its position on these axes
% (user's choice, November 2026: continuous scores instead of classes, with the second axis of the waveform added):
%   waveformScore1      first principal component of the 10 core waveform measures (44% of their variance): + broader and slower spikes
%                       (half-width, phase-plot area and steepest descent follow it with |rho| >= 0.9)
%   waveformScore2      second principal component (29%): + larger second peak and faster recovery after the trough, - larger first peak.
%                       Its negative end points toward the stereotyped shape that the quality control removes (the "QC noise shape"): the
%                       more a template resembles that shape, the lower this score (within-organoid rho -0.74). Whenever this score is used,
%                       add noiseSimilarity as a covariate and repeat the result on the units with clean = true (user's choice)
%   firingPatternScore  first principal component of the 4 features of "FiringPatternFeatures" (47%), measured outside the network bursts:
%                       + higher relative rate, higher CV, lower LV and lower burst index (firing in regular episodes separated by pauses),
%                       - irregular, bursty firing. It includes the relative rate, so it partly follows the firing rate
% The three scores are almost unrelated to each other (|r| <= 0.10). Their axes are fixed (section REFERENCE; user's choice, 13 November
% 2026): they were computed once on the 7,963 somatic units of 58 organoids (Kilosort 2.5), without the look-alikes of the QC noise shape,
% with every organoid weighted equally, and every recording is placed on them, so that the scores of recordings analyzed at different
% times, of different groups and of different projects are comparable (recordings can be added at any time; the scores of the others do
% not change). Each score is in standard deviations of the reference units (0 = their mean, every organoid weighted equally). The axes are
% stable: recomputed on halves of the organoids, on the organoids of one dataset, or with every unit counting once, the loadings agree with
% |cosine| >= 0.92. Fitted on a single group of 5-10 organoids, they would move the differences between groups by up to 0.3-0.65 SD of the
% organoid means: axes must not be refitted on one of the groups being compared. NeuronScore.reference.version names the axes used
% The waveform axes only apply to templates like those of the reference (Kilosort 2.5 export: REF_SAMPLES samples after the padding, at
% about REF_SAMPLING_RATE Hz); for other templates (or without templates, e.g., 3Brain units) the waveform scores are NaN
% Units scored: the waveform scores for the somatic units (not "nonSomatic" in "WaveformFeatures") that are not look-alikes
% (noiseSimilarity > LOOKALIKE_SIMILARITY; user's choice: they are kept in the network burst analysis, but their template says nothing
% about the cell type) and have the 10 measures; the firing-pattern score for these units when they have the 4 features (at least MIN_SPIKES
% spikes outside the bursts, see "FiringPatternFeatures"). Without usable templates, the firing-pattern score is given to every unit with
% the 4 features, with the relative rate computed against all units
% WAVEFORM_WEIGHT (user's choice) is only used where the scores are combined into one distance ("weightedScores"): it does not change the
% three scores themselves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   spikeTrain:     units of the baseline condition (as PerCondition(cond).spikeTrain of "Analyze_spikeTrain"), or the output of
%                   "FiringPatternFeatures" for these units
%   tRange:         [start end] of the analysis window, in sec (as PerCondition(cond).window); not used if spikeTrain is the output of
%                   "FiringPatternFeatures"
%   netBursts:      primary network bursts of the window (structure NetworkBurstDetected, or [start last] in sec, or []; see
%                   "FiringPatternFeatures"); not used if spikeTrain is the output of "FiringPatternFeatures"
%   waveforms:      Kilosort template of each unit on its peak electrode (as given by "loadKilosortTemplates", same order as spikeTrain),
%                   or the output of "WaveformFeatures" for these units; [] if there are no templates
%   samplingRate:   sampling rate of the recording, in Hz (not used if waveforms is the output of "WaveformFeatures" or [])
%   unitType:       (optional) type of each unit given by bombcell (see "WaveformFeatures")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scores of each unit (column vectors, one value per unit; NaN when the unit is not scored):
% NeuronScore.waveformScore1(i)        position on the first waveform axis, in SD of the reference units
% NeuronScore.waveformScore2(i)        position on the second waveform axis
% NeuronScore.firingPatternScore(i)    position on the firing-pattern axis
% NeuronScore.weightedScores(i,:)      the three scores scaled for a combined distance with WAVEFORM_WEIGHT: the waveform axes get
%                                      2*WAVEFORM_WEIGHT of variance along the first one (and the same share of their own variance along the
%                                      second one), the firing pattern 2*(1-WAVEFORM_WEIGHT); the Euclidean distance between two rows weighs
%                                      the waveform WAVEFORM_WEIGHT / (1-WAVEFORM_WEIGHT) times as much as the firing pattern
% Information about each unit:
% NeuronScore.noiseSimilarity(i)       correlation of the template (scaled to trough = -1) with the QC noise shape (NaN without template)
% NeuronScore.lookAlike(i)             =true if noiseSimilarity > LOOKALIKE_SIMILARITY (not scored)
% NeuronScore.clean(i)                 =true if noiseSimilarity <= CLEAN_SIMILARITY (the units on which to repeat any result of waveformScore2)
% NeuronScore.nonSomatic(i)            =true for the non-somatic units of "WaveformFeatures" (not scored)
% NeuronScore.waveformScored(i)        =true if the unit has the waveform scores
% NeuronScore.firingScored(i)          =true if the unit has the firing-pattern score
% NeuronScore.templatesUsed            =true if the templates of the recording could be placed on the reference axes
% Composition of the recording (means of the scored units; NaN with less than MIN_COMPOSITION_UNITS of them):
% NeuronScore.composition              structure with nUnits, nSomatic, nLookAlike, nWaveformScored, nFiringScored, nClean,
%                                      meanWaveformScore1, meanWaveformScore2, meanFiringPatternScore, meanNoiseSimilarity (of the scored
%                                      units), meanWaveformScore2Clean (of the scored units with clean = true)
% Details:
% NeuronScore.WaveformFeature          output of "WaveformFeatures" ([] without templates)
% NeuronScore.FiringPattern            output of "FiringPatternFeatures"
% NeuronScore.scores                   names of the three scores, in the order of the columns of weightedScores
% NeuronScore.reference                the reference axes used (section REFERENCE), with their version (.version)
% NeuronScore.parameters               parameters used

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WAVEFORM_WEIGHT = 2/3;              % (user's choice) weight of the waveform in a combined distance ("weightedScores"); the firing pattern gets
                                    % 1 - WAVEFORM_WEIGHT (2/3: the waveform counts twice as much); it does not change the three scores
LOOKALIKE_SIMILARITY = 0.95;        % (user's choice) units whose template correlates more than this with the QC noise shape are not scored
CLEAN_SIMILARITY = 0.7;             % (user's choice) units with a similarity up to this value are "clean": waveformScore2 barely depends on the
                                    % similarity among them (within-organoid rho -0.15, against -0.74 for all units)
MIN_COMPOSITION_UNITS = 10;         % minimum number of scored units for the means of "composition"
apply_parameters('NeuronScores');                                           % values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   REFERENCE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Axes computed once (13 November 2026) on the 7,963 somatic units of 58 organoids without the look-alikes (waveform), and on the 5,180 of
% them with the 4 firing-pattern features (outside the bursts, with the verdict of the pipeline for rejected bursts). Every organoid has the
% same weight (each of its n units counts 1/n; user's choice: large recordings do not dominate the axes). Each measure is transformed
% (WF_TRANSFORM, FP_TRANSFORM), standardized with the weighted mean and SD of the reference units, and projected on the principal axes of
% their weighted correlation matrix (scaled so that the scores of the reference units have weighted mean 0 and SD 1). With every unit
% counting once instead, the axes agree with |cosine| >= 0.99. They are not parameters: changing them would make the scores of different
% recordings incomparable. A new set of axes must get a new REFERENCE_VERSION
REFERENCE_VERSION = "2026-11-13, 58 organoids, equal weight per organoid";
WF_MEASURES = {'troughToPeak', 'halfWidth', 'peakAfterRatio', 'peakBeforeRatio', 'asymmetry', 'maxDescent', 'maxRecovery', ...
               'descentRecoveryRatio', 'slopeTroughToPeak', 'phasePlotArea'};       % the 10 core measures of "WaveformFeatures"
WF_TRANSFORM = {'log', 'log', 'log', '', '', 'log(-x)', '', 'log', 'log', 'log'};    % natural log of the skewed measures (maxDescent < 0)
WF_MEAN = [-0.3554245095 -1.325475078 -1.100911239 0.04769946803 0.7582337047 1.756306768 5.096373602 0.1406314711 -1.520315925 3.449709662];
WF_SD = [0.2094981082 0.1715143043 0.2312192704 0.05026080092 0.2060075887 0.234363527 0.8019693358 0.2396345612 0.2003764028 0.3361374948];
WF_AXES = [[0.08813886702 0.2188018137 0.1445013584 0.001331194543 0.02623917823 -0.2022885356 -0.1329343182 -0.1057277324 0.1788318354 -0.2222603561]' ...
           [-0.1965547372 -0.06446337537 0.1596232199 -0.2782000278 0.2985458415 -0.07652324176 0.2303102045 -0.23132724 -0.0591309757 -0.009686549157]'];
WF_EIGENVALUES = [4.378961801 2.850668345];         % variance of the first two principal components of the 10 standardized measures
FP_FEATURES = {'relativeRate', 'lv', 'cv', 'burstIndex'};                   % the 4 features of "FiringPatternFeatures"
FP_TRANSFORM = {'', '', 'log', ''};
FP_MEAN = [0.2450965648 0.8588435806 0.5690304987 -0.3053691768];
FP_SD = [0.5499248034 0.3903916109 0.6119461611 0.4191643091];
FP_AXIS = [0.2816459799 -0.3519558472 0.3532207892 -0.445992107]';
NOISE_SHAPE = [0.1565210108 0.141434632 0.1496622743 0.1551954213 0.1577640361 0.1576077911 0.1555531746 0.1526028806 0.1500158857 ...
    0.1469970728 0.1426464707 0.1320495421 0.1022909299 0.02969489272 -0.1134240641 -0.3225864009 -0.5494611123 -0.7473685369 ...
    -0.8978422241 -0.9914604965 -1 -0.921634841 -0.7997240356 -0.670705226 -0.5555430857 -0.4584712872 -0.3775868807 -0.3089220346 ...
    -0.2494651749 -0.1981466887 -0.1496778542 -0.1073891362 -0.06722832606 -0.03177653935 -0.0002807467258 0.02616722814 0.05132831673 ...
    0.07552434948 0.09881153656 0.1207062826 0.1393413713 0.1547758061 0.1679229793 0.1784345433 0.1891114909 0.1981114333 0.2073784233 ...
    0.2136400125 0.2174219565 0.2194454891 0.2201458556 0.2205321304 0.2211390991 0.2201139108 0.2178944097 0.2153748024 0.2118629746 ...
    0.2083753986 0.2031948251 0.1973172795 0.1830521783];       % QC noise shape: median template of the densest group of units removed by
                                                                % the quality control (scaled to trough = -1, without the padding)
REF_SAMPLES = numel(NOISE_SHAPE);   % samples of the reference templates after removing the padding
REF_SAMPLING_RATE = 19753.8;        % sampling rate of the reference recordings, in Hz (templates within 1% of it are accepted)

SCORES = {'waveformScore1', 'waveformScore2', 'firingPatternScore'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Waveform: measures, similarity to the QC noise shape, scores

if nargin < 6
    unitType = [];
end
if isstruct(spikeTrain) && isfield(spikeTrain, 'lv') && isfield(spikeTrain, 'relativeRate')
    n_units = numel(spikeTrain.lv);                                         % output of "FiringPatternFeatures"
else
    n_units = numel(spikeTrain);
end

if isstruct(waveforms) && isfield(waveforms, 'troughToPeak')
    W = waveforms;                                                          % output of "WaveformFeatures"
elseif isempty(waveforms)
    W = [];
else
    W = WaveformFeatures(waveforms, samplingRate, unitType);
end
templates_used = ~isempty(W);
if templates_used && numel(W.troughToPeak) ~= n_units
    error('NeuronScores: %d units in the waveforms, but %d in the spike trains', numel(W.troughToPeak), n_units);
end
if templates_used && size(W.waveform, 2) ~= REF_SAMPLES
    warning('NeuronScores: the templates have %d samples (without padding), the reference ones %d: no waveform scores', ...
        size(W.waveform, 2), REF_SAMPLES);
    templates_used = false;
end
if templates_used && isfield(W, 'parameters') && isfield(W.parameters, 'samplingRate') && ...
        abs(W.parameters.samplingRate / REF_SAMPLING_RATE - 1) > 0.01
    warning('NeuronScores: sampling rate %.1f Hz, reference %.1f Hz: no waveform scores', W.parameters.samplingRate, REF_SAMPLING_RATE);
    templates_used = false;
end

similarity = nan(n_units, 1);
non_somatic = false(n_units, 1);
wf_scores = nan(n_units, 2);
if templates_used
    X = W.waveform - mean(W.waveform, 2);                                   % Pearson correlation with the QC noise shape
    r0 = NOISE_SHAPE - mean(NOISE_SHAPE);
    similarity = (X * r0') ./ (sqrt(sum(X.^2, 2)) * norm(r0));
    non_somatic = logical(W.nonSomatic(:));
    M = nan(n_units, numel(WF_MEASURES));
    for m=1:numel(WF_MEASURES)
        M(:,m) = transform(W.(WF_MEASURES{m})(:), WF_TRANSFORM{m});
    end
    Z = (M - WF_MEAN) ./ WF_SD;
    wf_scores = Z * WF_AXES;
end
look_alike = similarity > LOOKALIKE_SIMILARITY;
waveform_scored = templates_used & ~non_somatic & ~look_alike & all(isfinite(wf_scores), 2);
wf_scores(~waveform_scored, :) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Firing pattern: features outside the network bursts, score

if templates_used
    reference = ~non_somatic & ~look_alike;                                 % units that set the median rate of the organoid
else
    reference = true(n_units, 1);
end
if isstruct(spikeTrain) && isfield(spikeTrain, 'lv') && isfield(spikeTrain, 'relativeRate')
    FP = spikeTrain;
else
    FP = FiringPatternFeatures(spikeTrain, tRange, netBursts, reference);
end
F = nan(n_units, numel(FP_FEATURES));
for m=1:numel(FP_FEATURES)
    F(:,m) = transform(FP.(FP_FEATURES{m})(:), FP_TRANSFORM{m});
end
fp_score = ((F - FP_MEAN) ./ FP_SD) * FP_AXIS;
firing_scored = reference & isfinite(fp_score);
fp_score(~firing_scored) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output

clean = similarity <= CLEAN_SIMILARITY;                                     % (false when the similarity is NaN)
weights = sqrt([2*WAVEFORM_WEIGHT, 2*WAVEFORM_WEIGHT*WF_EIGENVALUES(2)/WF_EIGENVALUES(1), 2*(1-WAVEFORM_WEIGHT)]);

NeuronScore.waveformScore1 = wf_scores(:,1);
NeuronScore.waveformScore2 = wf_scores(:,2);
NeuronScore.firingPatternScore = fp_score;
NeuronScore.weightedScores = [wf_scores fp_score] .* weights;
NeuronScore.noiseSimilarity = similarity;
NeuronScore.lookAlike = look_alike;
NeuronScore.clean = clean;
NeuronScore.nonSomatic = non_somatic;
NeuronScore.waveformScored = waveform_scored;
NeuronScore.firingScored = firing_scored;
NeuronScore.templatesUsed = templates_used;

C.nUnits = n_units;
C.nSomatic = NaN;   C.nLookAlike = NaN;
if templates_used
    C.nSomatic = sum(~non_somatic);
    C.nLookAlike = sum(look_alike & ~non_somatic);
end
C.nWaveformScored = sum(waveform_scored);
C.nFiringScored = sum(firing_scored);
C.nClean = sum(waveform_scored & clean);
C.meanWaveformScore1 = mean_if(wf_scores(waveform_scored, 1), MIN_COMPOSITION_UNITS);
C.meanWaveformScore2 = mean_if(wf_scores(waveform_scored, 2), MIN_COMPOSITION_UNITS);
C.meanFiringPatternScore = mean_if(fp_score(firing_scored), MIN_COMPOSITION_UNITS);
C.meanNoiseSimilarity = mean_if(similarity(waveform_scored), MIN_COMPOSITION_UNITS);
C.meanWaveformScore2Clean = mean_if(wf_scores(waveform_scored & clean, 2), MIN_COMPOSITION_UNITS);
NeuronScore.composition = C;

if templates_used
    NeuronScore.WaveformFeature = W;
else
    NeuronScore.WaveformFeature = [];
end
NeuronScore.FiringPattern = FP;
NeuronScore.scores = SCORES;
NeuronScore.reference = struct('version', REFERENCE_VERSION, 'WF_MEASURES', {WF_MEASURES}, 'WF_TRANSFORM', {WF_TRANSFORM}, 'WF_MEAN', WF_MEAN, 'WF_SD', WF_SD, ...
    'WF_AXES', WF_AXES, 'WF_EIGENVALUES', WF_EIGENVALUES, 'FP_FEATURES', {FP_FEATURES}, 'FP_TRANSFORM', {FP_TRANSFORM}, 'FP_MEAN', FP_MEAN, ...
    'FP_SD', FP_SD, 'FP_AXIS', FP_AXIS, 'NOISE_SHAPE', NOISE_SHAPE, 'REF_SAMPLING_RATE', REF_SAMPLING_RATE);
NeuronScore.parameters = struct('WAVEFORM_WEIGHT', WAVEFORM_WEIGHT, 'LOOKALIKE_SIMILARITY', LOOKALIKE_SIMILARITY, ...
    'CLEAN_SIMILARITY', CLEAN_SIMILARITY, 'MIN_COMPOSITION_UNITS', MIN_COMPOSITION_UNITS);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = transform(x, kind)
% transformation of a measure before the standardization (values out of the domain of the logarithm -> NaN)
switch kind
    case 'log'
        x(~(x > 0)) = NaN;      y = log(x);
    case 'log(-x)'
        x(~(x < 0)) = NaN;      y = log(-x);
    otherwise
        y = x;
end
end

function m = mean_if(x, n_min)
% mean of x when it has at least n_min values (NaN otherwise)
m = NaN;
if numel(x) >= n_min
    m = mean(x);
end
end
