function NeuronAxes = FitNeuronAxes(sources, groups, axesFile, version)

% function to fit a new set of axes for the neuron scores of "NeuronScores", on the scores saved by many recordings (new, November 2026)
% The axes of "NeuronScores" are fixed, so that every recording, analyzed at any time, is placed on the same axes and the scores of all the
% recordings are comparable. They were fitted once on 58 organoids of all the groups (section REFERENCE of "NeuronScores"). New axes are
% only needed for data of another kind (e.g., templates of another spike sorting or of another sampling rate), and they must be fitted on a
% broad set: many organoids of ALL the groups that will be compared (user's choice, 13 November 2026: axes fitted on a single group of 5-10
% organoids moved the differences between groups by up to 0.3-0.65 SD of the organoid means). It is not a stage of the pipeline:
% "CompareNeurons" runs it on all the groups of a comparison before comparing them (user's choice, 13 November 2026), and places the units
% on the new axes with AXES_FILE of "NeuronScores"; it can also be run by hand on the NeuronScore of analyzed recordings
% The axes are fitted as those of section REFERENCE. Units used: for the waveform axes, the somatic units that are not look-alikes of the QC
% noise shape and have the 10 measures; for the firing-pattern axis, the same units (all the units of a recording without templates) with the
% 4 features, the relative rate measured against them. Every organoid (source) has the same weight: each of its n units counts 1/n. Each
% measure is transformed as in "NeuronScores", standardized with the weighted mean and SD, and projected on the principal axes of the
% weighted correlation matrix, scaled so that the scores of these units have weighted mean 0 and SD 1. Signs: waveform score 1 grows with
% the half-width, waveform score 2 with the asymmetry (larger second peak), the firing-pattern score with the CV
% Checks: at least MIN_ORGANOIDS organoids of at least MIN_GROUPS groups for each set of axes (otherwise nothing is fitted), and the
% stability of each axis between random halves of the organoids (|cosine| of the loadings, N_SPLITS splits; a warning below MIN_STABILITY)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   sources:    "_Analyzed.mat" files with NeuronScore (one per organoid; string array or cell array of paths), or a cell array of
%               NeuronScore structures (output of "NeuronScores")
%   groups:     group of each source (e.g., the column "group" of the recordings table)
%   axesFile:   (optional) .mat file where the axes are written (variable NeuronAxes); an existing file is not replaced; "" or omitted: not
%               written
%   version:    (optional) name of the new axes (saved with the scores of every recording); default "<date>, <n> organoids of <g> groups,
%               equal weight per organoid". It must differ from the version of the axes the sources were scored with

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NeuronAxes:   the axes, with the fields of section REFERENCE of "NeuronScores" (version, WF_MEASURES, WF_TRANSFORM, WF_MEAN, WF_SD,
%               WF_AXES, WF_EIGENVALUES, FP_FEATURES, FP_TRANSFORM, FP_MEAN, FP_SD, FP_AXIS, NOISE_SHAPE, REF_SAMPLES, REF_SAMPLING_RATE;
%               NOISE_SHAPE is [] when the templates are not like those of the QC noise shape: then no unit is a look-alike), and
% NeuronAxes.fit                    how they were fitted:
%   .date, .sources, .groups        date, name and group of each source
%   .organoids                      table, one row per source: its group, its units for each set of axes and whether it was used
%   .nOrganoidsWaveform, .nUnitsWaveform, .nOrganoidsFiring, .nUnitsFiring, .nGroups
%   .eigenvaluesWaveform            the 10 eigenvalues of the weighted correlation matrix of the measures (sum = 10)
%   .eigenvaluesFiring              the 4 of the features (sum = 4)
%   .stabilityWaveform              |cosine| between the loadings of two random halves of the organoids (one row per split, one column per axis)
%   .stabilityFiring                the same for the firing-pattern axis
%   .parameters                     parameters used

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MIN_ORGANOIDS = 20;                 % minimum number of organoids for each set of axes (the axes of section REFERENCE: 58 for the waveform,
                                    % 51 for the firing pattern)
MIN_GROUPS = 2;                     % minimum number of groups among them (axes must not be fitted on one of the groups being compared)
MIN_ORGANOID_UNITS = 1;             % organoids with fewer units for a set of axes are left out of it (1: every organoid with units, as for
                                    % the axes of section REFERENCE, where the smallest had 2 units)
N_SPLITS = 50;                      % random splits of the organoids in two halves (within each group) for the stability of the axes
MIN_STABILITY = 0.9;                % warning when the median |cosine| between the axes of two halves is lower
SEED = 20261113;                    % seed of the random splits (the same splits at every run)
SAME_RATE = 0.01;                   % sampling rates within 1% are the same (as in "NeuronScores")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sources

if nargin < 3
    axesFile = "";
end
if nargin < 4
    version = "";
end
axesFile = string(axesFile);
version = string(version);
[NS, names] = read_sources(sources);
n_src = numel(NS);
groups = string(groups(:));
if numel(groups) ~= n_src
    error("FitNeuronAxes: %d sources but %d groups", n_src, numel(groups));
end
if strlength(axesFile) > 0 && isfile(axesFile)
    error("FitNeuronAxes: %s already exists; choose another file (the scores of the recordings name the axes they used)", axesFile);
end

R = NS{1}.reference;                                                        % measures and transformations of the scores
for k=2:n_src
    if ~isequal(NS{k}.reference.WF_MEASURES, R.WF_MEASURES) || ~isequal(NS{k}.reference.WF_TRANSFORM, R.WF_TRANSFORM) || ...
            ~isequal(NS{k}.reference.FP_FEATURES, R.FP_FEATURES) || ~isequal(NS{k}.reference.FP_TRANSFORM, R.FP_TRANSFORM)
        error("FitNeuronAxes: %s was scored with other measures than %s", names(k), names(1));
    end
end
old_versions = unique(string(cellfun(@(x) char(x.reference.version), NS, 'UniformOutput', false)));
if strlength(version) > 0 && any(version == old_versions)
    error("FitNeuronAxes: the version '%s' is already used; new axes need a new version name", version);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Units of each source

n_wf = numel(R.WF_MEASURES);  n_fp = numel(R.FP_FEATURES);
WF = cell(n_src, 1);  FPF = cell(n_src, 1);
samples = nan(n_src, 1);  rate = nan(n_src, 1);
for k=1:n_src
    W = NS{k}.WaveformFeature;
    FP = NS{k}.FiringPattern;
    n = numel(FP.lv);
    if ~isempty(W)
        reference = ~logical(W.nonSomatic(:)) & ~logical(NS{k}.lookAlike(:));   % somatic units, not look-alikes
        M = nan(n, n_wf);
        for m=1:n_wf
            M(:,m) = transform(W.(R.WF_MEASURES{m})(:), R.WF_TRANSFORM{m});
        end
        WF{k} = M(reference & all(isfinite(M), 2), :);
        samples(k) = size(W.waveform, 2);
        if isfield(W, 'parameters') && isfield(W.parameters, 'samplingRate')
            rate(k) = W.parameters.samplingRate;
        end
    else
        reference = true(n, 1);                                             % no templates: every unit (as "NeuronScores")
        WF{k} = zeros(0, n_wf);
    end
    if isfield(FP, 'rateOutside') && isfield(FP, 'enoughSpikes')
        FP = FiringPatternFeatures(FP, [], [], reference);                  % relative rate against these units
    end
    F = nan(n, n_fp);
    for m=1:n_fp
        F(:,m) = transform(FP.(R.FP_FEATURES{m})(:), R.FP_TRANSFORM{m});
    end
    FPF{k} = F(reference & all(isfinite(F), 2), :);
end

% the templates of all the sources must be alike (same number of samples, same sampling rate)
has_wf = ~isnan(samples);
if ~any(has_wf)
    error("FitNeuronAxes: no source has templates; the waveform axes cannot be fitted");
end
REF_SAMPLES = samples(find(has_wf, 1));
REF_SAMPLING_RATE = median(rate(has_wf), 'omitnan');
if any(samples(has_wf) ~= REF_SAMPLES) || any(abs(rate(has_wf) / REF_SAMPLING_RATE - 1) > SAME_RATE)
    error("FitNeuronAxes: the templates of the sources differ (%s samples, %s Hz); fit separate axes for each kind of data", ...
        strjoin(string(unique(samples(has_wf)))', '/'), strjoin(string(unique(round(rate(has_wf), 1)))', '/'));
end
NOISE_SHAPE = [];                                                           % the QC noise shape, if the templates are like it
for k=find(has_wf)'
    r = NS{k}.reference;
    if numel(r.NOISE_SHAPE) == REF_SAMPLES && isfield(r, 'REF_SAMPLING_RATE') && abs(REF_SAMPLING_RATE / r.REF_SAMPLING_RATE - 1) <= SAME_RATE
        NOISE_SHAPE = r.NOISE_SHAPE(:)';
        break
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Organoids used, checks

n_units_wf = cellfun(@(x) size(x, 1), WF);
n_units_fp = cellfun(@(x) size(x, 1), FPF);
use_wf = n_units_wf >= max(MIN_ORGANOID_UNITS, 1);
use_fp = n_units_fp >= max(MIN_ORGANOID_UNITS, 1);
check_set(use_wf, groups, 'waveform axes', MIN_ORGANOIDS, MIN_GROUPS);
check_set(use_fp, groups, 'firing-pattern axis', MIN_ORGANOIDS, MIN_GROUPS);
n_groups = numel(unique(groups(use_wf | use_fp)));
if strlength(version) == 0
    version = sprintf("%s, %d organoids of %d groups, equal weight per organoid", string(datetime('today', 'Format', 'yyyy-MM-dd')), ...
        sum(use_wf | use_fp), n_groups);
    if any(version == old_versions)
        version = version + " (" + string(datetime('now', 'Format', 'HH:mm:ss')) + ")";
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Axes, stability

orient_wf = [find(strcmp(R.WF_MEASURES, 'halfWidth')), find(strcmp(R.WF_MEASURES, 'asymmetry'))];
orient_fp = find(strcmp(R.FP_FEATURES, 'cv'));
[A_wf, org_wf] = stack(WF, use_wf);
[A_fp, org_fp] = stack(FPF, use_fp);
[WF_MEAN, WF_SD, WF_AXES, lam_wf] = fit_axes(A_wf, org_wf, 2, orient_wf);
[FP_MEAN, FP_SD, FP_AXIS, lam_fp] = fit_axes(A_fp, org_fp, 1, orient_fp);

state = rng;
rng(SEED, 'twister');
stab_wf = nan(N_SPLITS, 2);  stab_fp = nan(N_SPLITS, 1);
for s=1:N_SPLITS
    stab_wf(s,:) = split_half(WF, use_wf, groups, 2, orient_wf);
    stab_fp(s,:) = split_half(FPF, use_fp, groups, 1, orient_fp);
end
rng(state);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output

NeuronAxes = struct('version', version, 'WF_MEASURES', {R.WF_MEASURES}, 'WF_TRANSFORM', {R.WF_TRANSFORM}, 'WF_MEAN', WF_MEAN, ...
    'WF_SD', WF_SD, 'WF_AXES', WF_AXES, 'WF_EIGENVALUES', lam_wf(1:2)', 'FP_FEATURES', {R.FP_FEATURES}, 'FP_TRANSFORM', {R.FP_TRANSFORM}, ...
    'FP_MEAN', FP_MEAN, 'FP_SD', FP_SD, 'FP_AXIS', FP_AXIS, 'NOISE_SHAPE', NOISE_SHAPE, 'REF_SAMPLES', REF_SAMPLES, ...
    'REF_SAMPLING_RATE', REF_SAMPLING_RATE);
fit.date = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
fit.sources = names;
fit.groups = groups;
fit.organoids = table(names, groups, n_units_wf, n_units_fp, use_wf, use_fp, 'VariableNames', ...
    {'Source', 'Group', 'WaveformUnits', 'FiringUnits', 'UsedWaveform', 'UsedFiring'});
fit.nOrganoidsWaveform = sum(use_wf);   fit.nUnitsWaveform = sum(n_units_wf(use_wf));
fit.nOrganoidsFiring = sum(use_fp);     fit.nUnitsFiring = sum(n_units_fp(use_fp));
fit.nGroups = n_groups;
fit.eigenvaluesWaveform = lam_wf';
fit.eigenvaluesFiring = lam_fp';
fit.stabilityWaveform = stab_wf;
fit.stabilityFiring = stab_fp;
fit.parameters = struct('MIN_ORGANOIDS', MIN_ORGANOIDS, 'MIN_GROUPS', MIN_GROUPS, 'MIN_ORGANOID_UNITS', MIN_ORGANOID_UNITS, ...
    'N_SPLITS', N_SPLITS, 'MIN_STABILITY', MIN_STABILITY, 'SEED', SEED);
NeuronAxes.fit = fit;

fprintf("\nFitNeuronAxes: %s\n", version);
fprintf("  waveform axes: %d units of %d organoids; %.0f%% and %.0f%% of the variance (next axis %.0f%%); halves of the organoids: " + ...
    "|cos| median %.3f and %.3f (lowest %.3f and %.3f)\n", fit.nUnitsWaveform, fit.nOrganoidsWaveform, 100 * lam_wf(1:3) / n_wf, ...
    median(stab_wf), min(stab_wf));
fprintf("  firing-pattern axis: %d units of %d organoids; %.0f%% of the variance (next axis %.0f%%); halves: |cos| median %.3f " + ...
    "(lowest %.3f)\n", fit.nUnitsFiring, fit.nOrganoidsFiring, 100 * lam_fp(1:2) / n_fp, median(stab_fp), min(stab_fp));
[g, ~, gi] = unique(groups(use_wf | use_fp));
fprintf("  organoids per group: %s\n", strjoin(g + " " + string(accumarray(gi, 1)), ', '));
if isempty(NOISE_SHAPE)
    fprintf("  the templates (%d samples at %.1f Hz) are not like the QC noise shape: no look-alikes with these axes\n", REF_SAMPLES, REF_SAMPLING_RATE);
end
names_axes = ["waveform score 1", "waveform score 2", "firing-pattern score"];
low = [median(stab_wf), median(stab_fp)] < MIN_STABILITY;
if any(low)
    warning("FitNeuronAxes: %s not stable between halves of the organoids (median |cos| below %g): more organoids are needed", ...
        strjoin(names_axes(low), ', '), MIN_STABILITY);
end
if strlength(axesFile) > 0
    save(axesFile, 'NeuronAxes');
    fprintf("  written in %s (AXES_FILE of NeuronScores)\n", axesFile);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NS, names] = read_sources(sources)
% NeuronScore of each source, and its name (file name without "_Analyzed.mat", or "source <k>")
if iscell(sources) && ~isempty(sources) && all(cellfun(@isstruct, sources))
    NS = sources(:);
    names = "source " + string(1:numel(NS))';
    return
end
files = string(sources);                                                    % (a path as text: one file)
files = files(:);
if numel(unique(files)) < numel(files)
    error("FitNeuronAxes: a file is given twice (every organoid must count once)");
end
NS = cell(numel(files), 1);
names = strings(numel(files), 1);
missing = strings(0, 1);
for k=1:numel(files)
    [~, name] = fileparts(files(k));
    names(k) = regexprep(name, '_Analyzed$', '');
    if ~isfile(files(k)) || isempty(whos('-file', files(k), 'NeuronScore'))
        missing(end+1) = files(k);                                          %#ok<AGROW>
        continue
    end
    NS{k} = load(files(k), 'NeuronScore').NeuronScore;
    if isempty(NS{k})
        missing(end+1) = files(k);                                          %#ok<AGROW>
    end
end
if ~isempty(missing)
    error("FitNeuronAxes: no neuron scores in %s (run Analyze_spikeTrain with NEURON_SCORES = 1 first)", strjoin(missing, ', '));
end
if numel(unique(names)) < numel(names)
    error("FitNeuronAxes: two files have the same name (every organoid must count once)");
end
end

function check_set(use, groups, what, min_organoids, min_groups)
% enough organoids and groups for a set of axes (otherwise nothing is fitted)
if sum(use) < min_organoids
    error("FitNeuronAxes: %d organoids with units for the %s, at least %d are needed (fit the axes on a broad set of all the groups)", ...
        sum(use), what, min_organoids);
end
if numel(unique(groups(use))) < min_groups
    error("FitNeuronAxes: the organoids of the %s come from %d group(s), at least %d are needed (axes must not be fitted on one of " + ...
        "the groups being compared)", what, numel(unique(groups(use))), min_groups);
end
end

function [A, org] = stack(X, use)
% units of the organoids used, with the index of their organoid
idx = find(use);
A = vertcat(X{idx});
org = repelem((1:numel(idx))', cellfun(@(x) size(x, 1), X(idx)));
end

function [mu, sd, ax, lam] = fit_axes(A, org, k, orient)
% weighted (each organoid 1, each of its n units 1/n) standardization and principal axes of the weighted correlation matrix, scaled so that
% the scores have weighted SD 1; the sign of axis j makes its scores grow with the measure orient(j)
w = 1 ./ accumarray(org, 1);
w = w(org);
w = w / sum(w);
mu = w' * A;
sd = sqrt(w' * (A - mu).^2);
Z = (A - mu) ./ sd;
C = (Z .* w)' * Z;
[V, L] = eig((C + C') / 2);
[lam, o] = sort(diag(L), 'descend');
V = V(:, o);
ax = V(:, 1:k) ./ sqrt(lam(1:k))';
for j=1:k
    r = corrcoef(Z * ax(:,j), A(:, orient(j)));
    if r(1,2) < 0
        ax(:,j) = -ax(:,j);
    end
end
end

function c = split_half(X, use, groups, k, orient)
% |cosine| between the axes fitted on two random halves of the organoids (split within each group)
idx = find(use);
half = false(numel(idx), 1);
g = groups(idx);
for u=unique(g)'
    in = find(g == u);
    in = in(randperm(numel(in)));
    first = mod(randi(2), 2);                                               % which half gets the extra organoid of an odd group
    half(in(1 + first:2:end)) = true;
end
use_a = false(size(use));  use_a(idx(half)) = true;
use_b = false(size(use));  use_b(idx(~half)) = true;
[A, org] = stack(X, use_a);
[~, ~, ax_a] = fit_axes(A, org, k, orient);
[B, org] = stack(X, use_b);
[~, ~, ax_b] = fit_axes(B, org, k, orient);
c = abs(sum(ax_a .* ax_b, 1)) ./ (vecnorm(ax_a) .* vecnorm(ax_b));
end

function y = transform(x, kind)
% transformation of a measure before the standardization (the same as in "NeuronScores")
switch kind
    case 'log'
        x(~(x > 0)) = NaN;      y = log(x);
    case 'log(-x)'
        x(~(x < 0)) = NaN;      y = log(-x);
    otherwise
        y = x;
end
end
