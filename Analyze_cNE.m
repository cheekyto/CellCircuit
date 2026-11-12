function cNEResults = Analyze_cNE(PerCondition, base_cond, analyzed, cNEResults, options)

% (new, November 2026) cNE step of "Analyze_spikeTrain": the coordinated neuronal ensembles (cNEs) of each condition of a recording, and
% their comparison with the baseline and with the condition recorded before
%   1. each condition: its cNEs, detected on its own analysis window ("cNEDetection"), and described against its network events
%      ("cNEAnalysis")
%   2. each other condition against the baseline ("cNECompare"): the baseline cNEs (detected once, in 1) followed in the condition, and
%      the cNEs of the condition (also from 1) matched to them
%   3. (user's choice) each condition against the condition recorded just before it, when that one is not the baseline (e.g., with a
%      baseline and three drugs, drug 2 against drug 1 and drug 3 against drug 2; see "previousCondition"), in the same way as 2
%   4. the baseline, its 2nd half against its 1st half ("cNECompare"): the change that time alone brings, to compare the drugs with
% (user's choice) The TTX conditions (label with "TTX") are left out of every step: no cNEs, no comparison, and never the condition before
% another one. The short list of measures compiled by "Compile_MEAproject" (sheets "cNEs", "cNEs vs baseline" and "cNEs vs previous") is
% described in docs/cNE_METRICS.md
% A step that fails is reported with a warning and its status ("error: <message>"); the other steps, and the rest of the pipeline, go on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   PerCondition:   results of "Analyze_spikeTrain" for the recording (spikeTrain, window, NetworkBurstDetected, condition, dose and
%                   too_few_units of each condition)
%   base_cond:      index of the baseline condition
%   analyzed:       indices of the conditions analyzed in this run: their results, and the comparisons that involve them, are computed
%                   again; the results of the other conditions are kept from cNEResults (and computed if they are missing)
%   cNEResults:     (optional) results of a previous run; [] = none
%   options:        (optional) structure with any of the fields
%       plots           true = draw the figures and save them as .fig files (default: false)
%       folder          folder of the figures
%       file_data       name of the recording, at the start of the name of the figures
%       hide            true = the figures are created invisible, saved so that they open visible, and closed (default: false)
%       detection, analysis, compare    inputs of "cNEDetection", "cNEAnalysis" and "cNECompare" (default: none; the comparisons use
%                       the same inputs of "cNEDetection" and "cNEAnalysis")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cNEResults       structure:
%   .condition(c)       each condition: name (condition and dose), window, status ("ok", "excluded (TTX)", or why there is no
%                       result), cNE (output of "cNEDetection"), report (output of "cNEAnalysis")
%   .comparison(c)      each condition against the baseline (empty for the baseline): name, status, result (output of "cNECompare")
%   .previous(c)        each condition against the condition recorded just before it: reference (index of that condition; 0 = no
%                       such comparison, see "previousCondition"), name, status, result (output of "cNECompare")
%   .halves             the 2nd half of the baseline against its 1st half: name, status, result (output of "cNECompare")
%   .baseline           index of the baseline condition
%   .created            date and time of the run
% Figures (plots = true): <file_data>_cNE_<condition><dose>.fig (each condition), <file_data>_cNECompare_<condition><dose>.fig (each
% comparison with the baseline), <file_data>_cNECompare_<condition><dose>_vs_<previous condition><dose>.fig (each comparison with the
% condition before), <file_data>_cNECompare_BaselineHalves.fig; when two conditions have the same name, the index of the condition is added
% (e.g., _cNE_+TTX_2.fig and _cNE_+TTX_3.fig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4, cNEResults = []; end
if nargin < 5, options = struct(); end
opts = struct('plots', false, 'folder', pwd, 'file_data', "recording", 'hide', false, 'detection', struct(), 'analysis', struct(), ...
    'compare', struct());
f = fieldnames(options);
for i=1:numel(f)
    opts.(f{i}) = options.(f{i});
end
cmp_inputs = opts.compare;                                                  % the comparisons detect and describe the cNEs as in 1
if ~isfield(cmp_inputs, 'DETECTION'), cmp_inputs.DETECTION = opts.detection; end
if ~isfield(cmp_inputs, 'ANALYSIS'), cmp_inputs.ANALYSIS = opts.analysis; end
if opts.plots && opts.hide                                                  % figures created invisible (they would take over the screen)
    previous = get(groot, 'DefaultFigureVisible');
    set(groot, 'DefaultFigureVisible', 'off');
    restore = onCleanup(@() set(groot, 'DefaultFigureVisible', previous));
end

if numel(PerCondition) < base_cond                                          % the baseline was never analyzed (an empty condition)
    PerCondition(base_cond).window = [];
end
n = numel(PerCondition);
empty_cond = struct('name', "", 'window', [NaN NaN], 'status', "not analyzed", 'cNE', [], 'report', []);
empty_cmp = struct('name', "", 'status', "not compared", 'result', []);
empty_prev = struct('reference', 0, 'name', "", 'status', "not compared", 'result', []);
old = cNEResults;
same_baseline = isstruct(old) && isfield(old, 'baseline') && isequal(old.baseline, base_cond);

cNEResults = struct();
cNEResults.condition = repmat(empty_cond, 1, n);
cNEResults.comparison = repmat(empty_cmp, 1, n);
cNEResults.previous = repmat(empty_prev, 1, n);
cNEResults.halves = empty_cmp;
cNEResults.baseline = base_cond;
cNEResults.created = datetime('now');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. cNEs of each condition (kept from the previous run when the condition was not analyzed again)

tags = strings(1, n);                                                       % names of the figures: condition and dose, and the index of
for c=1:n                                                                   % the condition when two conditions have the same name
    tags(c) = as_text(PerCondition(c).condition) + as_text(PerCondition(c).dose);
end
dup = arrayfun(@(c) sum(tags == tags(c)) > 1, 1:n);
tags(dup) = tags(dup) + "_" + string(find(dup));

ttx = arrayfun(@(c) contains(as_text(PerCondition(c).condition), "TTX", 'IgnoreCase', true), 1:n);   % (user's choice) left out

changed = false(1, n);
changed(analyzed(analyzed >= 1 & analyzed <= n)) = true;
for c=1:n
    if ttx(c)
        cNEResults.condition(c).name = tags(c);
        cNEResults.condition(c).status = "excluded (TTX)";
        if ~isempty(PerCondition(c).window)
            cNEResults.condition(c).window = PerCondition(c).window;
        end
    elseif ~changed(c) && same_baseline && c <= numel(old.condition) && ismember(old.condition(c).status, ["ok" "too few units"])
        cNEResults.condition(c) = old.condition(c);
    else
        changed(c) = true;
        cNEResults.condition(c) = own_cnes(PerCondition(c), empty_cond, opts, tags(c));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Each other condition against the baseline

B = cNEResults.condition(base_cond);
for c=setdiff(1:n, base_cond)
    if ~ttx(c) && ~changed(c) && ~changed(base_cond) && same_baseline && c <= numel(old.comparison) && ...
            ~startsWith(old.comparison(c).status, "error")
        cNEResults.comparison(c) = old.comparison(c);
        continue
    end
    K = empty_cmp;
    K.name = B.name + " vs " + cNEResults.condition(c).name;
    if ttx(c)
        K.status = "excluded (TTX)";
    elseif B.status ~= "ok"
        K.status = "no baseline cNEs (" + B.status + ")";
    elseif cNEResults.condition(c).status ~= "ok"
        K.status = "no cNEs for the condition (" + cNEResults.condition(c).status + ")";
    else
        ref = condition_input(PerCondition(base_cond), PerCondition(base_cond).window, B.name, B.cNE);
        test = condition_input(PerCondition(c), PerCondition(c).window, cNEResults.condition(c).name, cNEResults.condition(c).cNE);
        K = run_compare(K, ref, test, cmp_inputs, opts, "cNECompare_" + tags(c));
    end
    cNEResults.comparison(c) = K;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. (user's choice) Each condition against the condition recorded just before it (when that one is not the baseline)

prev = previousCondition({PerCondition.window}, base_cond, ttx);
for c=find(prev > 0)
    p = prev(c);
    if ~changed(c) && ~changed(p) && isstruct(old) && isfield(old, 'previous') && c <= numel(old.previous) && ...
            old.previous(c).reference == p && ~startsWith(old.previous(c).status, "error")
        cNEResults.previous(c) = old.previous(c);
        continue
    end
    K = empty_prev;
    K.reference = p;
    K.name = cNEResults.condition(p).name + " vs " + cNEResults.condition(c).name;
    if cNEResults.condition(p).status ~= "ok"
        K.status = "no cNEs for the previous condition (" + cNEResults.condition(p).status + ")";
    elseif cNEResults.condition(c).status ~= "ok"
        K.status = "no cNEs for the condition (" + cNEResults.condition(c).status + ")";
    else
        ref = condition_input(PerCondition(p), PerCondition(p).window, cNEResults.condition(p).name, cNEResults.condition(p).cNE);
        test = condition_input(PerCondition(c), PerCondition(c).window, cNEResults.condition(c).name, cNEResults.condition(c).cNE);
        K = run_compare(K, ref, test, cmp_inputs, opts, "cNECompare_" + tags(c) + "_vs_" + tags(p));
    end
    cNEResults.previous(c) = K;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. The baseline: 2nd half against the 1st half

if ~changed(base_cond) && same_baseline && isfield(old, 'halves') && ~startsWith(old.halves.status, "error")
    cNEResults.halves = old.halves;
else
    K = empty_cmp;
    K.name = B.name + " 1st half vs 2nd half";
    if B.status ~= "ok"
        K.status = "no baseline cNEs (" + B.status + ")";
    else
        w = PerCondition(base_cond).window;     mid = mean(w);
        ref = condition_input(PerCondition(base_cond), [w(1) mid], B.name + " 1st half", []);
        test = condition_input(PerCondition(base_cond), [mid w(2)], B.name + " 2nd half", []);
        K = run_compare(K, ref, test, cmp_inputs, opts, "cNECompare_BaselineHalves");
    end
    cNEResults.halves = K;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions

function R = own_cnes(P, R, opts, tag)
% cNEs of one condition ("cNEDetection") and their description ("cNEAnalysis")
R.name = as_text(P.condition) + as_text(P.dose);
if ~isempty(P.window)
    R.window = P.window;
end
if isfield(P, 'too_few_units') && isequal(P.too_few_units, true)
    R.status = "too few units";
    return
elseif isempty(P.spikeTrain) || numel(P.window) ~= 2 || any(~isfinite(P.window))
    R.status = "not analyzed";
    return
end
t0 = tic;
try
    R.cNE = cNEDetection(P.spikeTrain, P.window, opts.detection);
    if opts.plots
        [R.report, fig] = cNEAnalysis(R.cNE, P.NetworkBurstDetected, P.spikeTrain, opts.analysis);
        save_figure(fig, opts, "cNE_" + tag);
    else
        R.report = cNEAnalysis(R.cNE, P.NetworkBurstDetected, P.spikeTrain, opts.analysis);
    end
    R.status = "ok";
    fprintf("cNEs of %s: %g cNEs (%.0f sec) \n", R.name, R.cNE.nCNE, toc(t0));
catch err
    R.status = "error: " + string(err.message);
    warning('Analyze_cNE:error', 'cNEs of %s not computed: %s', R.name, err.message);
end
end

function K = run_compare(K, ref, test, inputs, opts, figure_name)
% one comparison ("cNECompare")
t0 = tic;
try
    if opts.plots
        [K.result, fig] = cNECompare(ref, test, inputs);
        save_figure(fig, opts, figure_name);
    else
        K.result = cNECompare(ref, test, inputs);
    end
    K.status = "ok";
    s = K.result.summary;
    fprintf("cNEs, %s: %g found, %g lost, %g new (%.0f sec) \n", K.name, s.nFound, s.nLost, s.nNew, toc(t0));
catch err
    K.status = "error: " + string(err.message);
    warning('Analyze_cNE:error', 'cNEs, %s not compared: %s', K.name, err.message);
end
end

function s = as_text(x)
% a label as a string ("" when empty)
s = "";
if ~isempty(x) && ~all(ismissing(string(x)))
    s = strjoin(string(x), "");
end
end

function S = condition_input(P, window, name, cNE)
% one condition, as "cNECompare" takes it (cNE = [] -> detected by "cNECompare")
S = struct('spikeTimestamp', {P.spikeTrain}, 'tRange', window, 'NetworkBurstDetected', {P.NetworkBurstDetected}, 'name', name, ...
    'cNE', {cNE});
end

function save_figure(fig, opts, name)
% saves a figure as a .fig file (as "save_plot" of "Analyze_spikeTrain"): a hidden figure is saved so that it opens visible, and closed
file = fullfile(opts.folder, sprintf('%s_%s.fig', opts.file_data, name));
if opts.hide
    set(fig, 'CreateFcn', 'set(gcbo, ''Visible'', ''on'')');
end
saveas(fig, file);
if opts.hide
    close(fig);
end
end
