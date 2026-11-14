function Results = CompareNeurons(settings)

% function to compare groups of organoids by the neurons they are made of (new, November 2026; user's choices of 13 November 2026)
% One run, on groups that went through the pipeline:
%   1. the axes of the neuron scores are fitted again on all the organoids of FIT_GROUPS (every organoid weighted equally; "FitNeuronAxes";
%      at least 20 organoids of at least 2 groups, otherwise nothing is run), and every unit is placed on them ("NeuronScores"); the files
%      of the recordings are not changed, and the axes are saved with the results. Groups that only help the fitting can be listed
%   2. for each comparison chosen in the settings (2 or more groups; the first one is the reference):
%      - composition: do the groups differ in the mean scores of their organoids?
%      - activity: does the composition explain the activity? For each activity measure, activity ~ group + composition + covariates
%        (+ random intercepts), fitted with and without the composition: the group effect both ways (how much of it the composition
%        accounts for), the effect of the composition within the groups, and a likelihood-ratio test of the composition. Linear model, or
%        logistic for a yes/no measure; with random intercepts (columns of the recordings tables, e.g., "batch") it is a mixed model
%      - roles: within the organoids, do the scores relate to the roles of the units in the network bursts of the baseline (onset,
%        participation)? Every variable is turned into normal scores within each organoid (the organoids differ in size and activity),
%        and a mixed model with a random intercept and random slopes of the scores for each organoid gives the association across all the
%        organoids of the comparison, within each group, and its difference between the groups
% The organoids of a group are the "_Analyzed.mat" files of its folder <project>/<group> (as for "Compile_MEAproject"); the activity measures
% come from <project>/<group>_Compiled.mat (the condition chosen in the settings, the baseline by default). A recording without saved neuron
% scores (analyzed before NEURON_SCORES) gets them here, from its baseline and the Kilosort templates of its folder (column "folder" of
% the recordings table), without changing its file
% The settings are in a copy of "comparison_settings.m" (see there for each of them, with the default set of activity measures)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   settings:   path of a settings file (a copy of "comparison_settings.m"), or a structure with the same fields;
%               "defaults": returns the default settings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Results (also saved in <OUTPUT>/<NAME>_<date and time>.mat, with the Excel file of the same name; the axes in <...>_axes.mat):
%   .settings       the settings used
%   .axes           the axes fitted (output of "FitNeuronAxes"), .axesAgreement: |cosine| with the fixed axes of "NeuronScores"
%   .Organoids      one row per organoid of FIT_GROUPS: group, units, composition (mean scores), activity, covariates (sheet "Organoids")
%   .Composition    group differences of the mean scores, for each comparison (sheet "Composition")
%   .Activity       models of the activity, for each comparison, measure and composition score (sheet "Activity")
%   .Roles          models of the roles in the network bursts, for each comparison, role and score (sheet "Roles")
%   .Units          one row per unit: scores on the new axes and roles (sheet "Units")
%   .files          files written

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MIN_FIT_ORGANOIDS = 20;             % (user's choice) minimum number of organoids in FIT_GROUPS
MIN_ORGANOIDS_PER_TERM = 10;        % a note when a model of the activity has fewer organoids per term (rule of thumb; nothing is blocked)
MIN_RANDOM_LEVELS = 3;              % a note when a random intercept has fewer levels among the organoids of a model
MIN_GROUP_ORGANOIDS = 3;            % role models within a group need at least this many organoids
LEVEL = 0.95;                       % confidence intervals
RATE_OFFSET = 0.01;                 % logRateOutside = log10(rate outside the bursts + RATE_OFFSET)
MIN_UNITS = 5;                      % recordings with fewer units have no neuron scores (as MIN_UNITS of "Analyze_spikeTrain")

if (ischar(settings) || isstring(settings)) && string(settings) == "defaults"
    Results = default_settings();
    return
end
C = read_settings(settings);
stamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
base = fullfile(C.OUTPUT, C.NAME + "_" + stamp);
if ~isfolder(C.OUTPUT)
    mkdir(C.OUTPUT);
end
fprintf("\nCompareNeurons: %s\n", C.NAME);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Organoids of the groups

[O, Tab] = find_organoids(C);
all_groups = unique(O.Group, 'stable');
fit_groups = C.FIT_GROUPS;
if isempty(fit_groups)
    fit_groups = all_groups;
end
compared = unique([C.COMPARISONS{:}], 'stable');
unknown = setdiff([fit_groups(:); compared(:)], all_groups);
if ~isempty(unknown)
    error("CompareNeurons: no organoids for the group(s) %s (groups found: %s)", strjoin(unknown, ', '), strjoin(all_groups, ', '));
end
outside = setdiff(compared, fit_groups);
if ~isempty(outside)
    error("CompareNeurons: the groups compared must be among FIT_GROUPS (the axes are fitted on them): %s", strjoin(outside, ', '));
end
O = O(ismember(O.Group, fit_groups), :);
dup = unique(O.Recording(duplicated_values(O.Recording)));
if ~isempty(dup)
    error("CompareNeurons: the same recording is in more than one group folder (every organoid must count once): %s", strjoin(dup, ', '));
end
if height(O) < MIN_FIT_ORGANOIDS
    error("CompareNeurons: %d organoids in FIT_GROUPS, at least %d are needed to fit the axes", height(O), MIN_FIT_ORGANOIDS);
end
fprintf("  %d organoids of %d groups; comparisons: %s\n", height(O), numel(fit_groups), comparison_labels(C));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Neuron scores and baseline of each organoid

n = height(O);
NS = cell(n, 1);  Base = cell(n, 1);
source = strings(n, 1);  note = strings(n, 1);
for i=1:n
    [NS{i}, Base{i}, source(i), note(i)] = load_organoid(O(i,:), MIN_UNITS);
end
fprintf("  neuron scores: %d saved, %d computed here, %d organoids without units\n", sum(startsWith(source, "saved")), ...
    sum(startsWith(source, "computed")), sum(source == "none"));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Axes fitted on FIT_GROUPS, every unit placed on them

has = ~cellfun(@isempty, NS);
version = sprintf("%s %s, %d organoids of %d groups, equal weight per organoid", C.NAME, stamp, sum(has), numel(unique(O.Group(has))));
axes_file = base + "_axes.mat";
NeuronAxes = FitNeuronAxes(NS(has), O.Group(has), axes_file, version);
fixed = fixed_reference();
agreement = [abs(cosine(NeuronAxes.WF_AXES(:,1), fixed.WF_AXES(:,1))), abs(cosine(NeuronAxes.WF_AXES(:,2), fixed.WF_AXES(:,2))), ...
    abs(cosine(NeuronAxes.FP_AXIS, fixed.FP_AXIS))];
fprintf("  agreement with the fixed axes (|cos|): %.3f, %.3f, %.3f\n", agreement);
for i=find(has)'
    NS{i} = rescore(NS{i}, axes_file);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Units (scores and roles) and organoids (composition, activity, covariates)

Units = units_table(O, NS, Base, RATE_OFFSET);
Org = organoid_table(O, NS, Base, source, note);
compiled = load_compiled(O, C);
[Org, activity_names] = add_activity(Org, O, Base, compiled, C);
[Org, covariate_names, covariate_kind] = add_covariates(Org, O, Base, compiled, Tab, C.COVARIATES, C);
[Org, random_names] = add_covariates(Org, O, Base, compiled, Tab, C.RANDOM, C, true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparisons

Composition = table();  Activity = table();  Roles = table();
opts = struct('LEVEL', LEVEL, 'MIN_ORGANOIDS_PER_TERM', MIN_ORGANOIDS_PER_TERM, 'MIN_RANDOM_LEVELS', MIN_RANDOM_LEVELS, ...
    'MIN_GROUP_ORGANOIDS', MIN_GROUP_ORGANOIDS);
for c=1:numel(C.COMPARISONS)
    groups = C.COMPARISONS{c};
    label = strjoin(groups, " vs ");
    fprintf("  %s\n", label);
    rows = ismember(Org.Group, groups);
    Composition = [Composition; compare_composition(Org(rows,:), groups, label, C, random_names, opts)];   %#ok<AGROW>
    for a=1:size(C.ACTIVITY, 1)
        Activity = [Activity; compare_activity(Org(rows,:), groups, label, C.ACTIVITY(a,:), activity_names(a), C, ...
            covariate_names, covariate_kind, random_names, opts)];         %#ok<AGROW>
    end
    for r=1:numel(C.ROLES)
        Roles = [Roles; compare_roles(Units(ismember(Units.Group, groups) & Units.PrimaryBursts >= C.MIN_BURSTS, :), groups, label, ...
            C.ROLES(r), C, opts)];                                          %#ok<AGROW>
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results

Results.settings = C;
Results.axes = NeuronAxes;
Results.axesAgreement = agreement;
Results.Organoids = Org;
Results.Composition = Composition;
Results.Activity = Activity;
Results.Roles = Roles;
Results.Units = Units;
Results.files = struct('excel', base + ".xlsx", 'mat', base + ".mat", 'axes', axes_file);
save(Results.files.mat, 'Results', '-v7.3');
write_excel(Results, base + ".xlsx");
fprintf("  results: %s (.xlsx and .mat; axes in %s)\n", base, axes_file);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   SETTINGS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = default_settings()
% the defaults (the same as "comparison_settings.m"; TABLES and COMPARISONS must be given)
D.TABLES = strings(0, 1);
D.PROJECTS = strings(0, 1);
D.CONDITION = "baseline";
D.FIT_GROUPS = strings(0, 1);
D.COMPARISONS = {};
D.ACTIVITY = [
    "Bursting",                     "NETbursting1.Number_NETburst_Events",  "yes/no",   "all"
    "Network-burst rate",           "NETbursting1.Rate_NETburst",           "log",      "bursting"
    "Network-burst participation",  "NETbursting1.Participation_NETburst",  "none",     "bursting"
    "Network-burst duration",       "NETbursting1.Duration_NETburst",       "log",      "bursting"
    "Firing rate",                  "Firing.Firing_Rate",                   "log",      "all"
    ];
D.COMPOSITION = ["MeanWaveformScore1", "MeanWaveformScore2", "MeanFiringPatternScore"];
D.COMPOSITION_MODE = "separate";
D.COVARIATES = "MeanNoiseSimilarity";
D.RANDOM = strings(0, 1);
D.ROLES = ["onset", "participation"];
D.ROLE_SCORES = ["waveformScore1", "waveformScore2", "firingPatternScore"];
D.ROLE_COVARIATES = ["noiseSimilarity", "logRateOutside"];
D.MIN_BURSTS = 5;
D.MIN_ROLE_UNITS = 20;
D.OUTPUT = strings(0, 1);
D.NAME = "comparison";
end

function C = read_settings(settings)
% the settings of the file (or structure), with the defaults for the ones not given, checked
D = default_settings();
here = string(pwd);
if isstruct(settings)
    S = settings;
else
    file = string(settings);
    if ~isfile(file)
        error("CompareNeurons: the settings file %s does not exist", file);
    end
    S = run_settings(file);
    here = string(fileparts(char(file)));
    if strlength(here) == 0, here = string(pwd); end
end
unknown = setdiff(fieldnames(S), fieldnames(D));
if ~isempty(unknown)
    error("CompareNeurons: unknown setting(s) %s", strjoin(unknown, ', '));
end
C = D;
for f=string(fieldnames(S))'
    C.(f) = S.(f);
end
for f=["TABLES", "PROJECTS", "CONDITION", "FIT_GROUPS", "COMPOSITION", "COMPOSITION_MODE", "COVARIATES", "RANDOM", "ROLES", ...
        "ROLE_SCORES", "ROLE_COVARIATES", "OUTPUT", "NAME"]
    C.(f) = reshape(string(C.(f)), 1, []);
end
C.ACTIVITY = string(C.ACTIVITY);
if isempty(C.TABLES)
    error("CompareNeurons: TABLES is needed (the recordings tables)");
end
if ~isempty(C.PROJECTS) && numel(C.PROJECTS) ~= numel(C.TABLES)
    error("CompareNeurons: PROJECTS needs one folder per table (or [] for the folders of the tables)");
end
if ~iscell(C.COMPARISONS) || isempty(C.COMPARISONS)
    error("CompareNeurons: COMPARISONS is needed, e.g. {[""groupA"", ""groupB""]}");
end
for k=1:numel(C.COMPARISONS)
    C.COMPARISONS{k} = reshape(string(C.COMPARISONS{k}), 1, []);
    if numel(unique(C.COMPARISONS{k})) < 2 || numel(unique(C.COMPARISONS{k})) < numel(C.COMPARISONS{k})
        error("CompareNeurons: each comparison needs 2 or more different groups (comparison %d)", k);
    end
end
if size(C.ACTIVITY, 2) ~= 4 || isempty(C.ACTIVITY)
    error("CompareNeurons: ACTIVITY needs one row per measure with 4 columns (name, table.column, transformation, organoids)");
end
tables = ["Firing", "Bursting", "NETbursting0", "NETbursting1", "NETbursting2", "cNEs"];
for a=1:size(C.ACTIVITY, 1)
    parts = split(C.ACTIVITY(a,2), ".");
    if numel(parts) ~= 2 || ~ismember(parts(1), tables)
        error("CompareNeurons: activity '%s': the column must be <table>.<column>, with a table among %s", C.ACTIVITY(a,1), strjoin(tables, ', '));
    end
    if ~ismember(C.ACTIVITY(a,3), ["none", "log", "yes/no"])
        error("CompareNeurons: activity '%s': the transformation must be none, log or yes/no", C.ACTIVITY(a,1));
    end
    if ~ismember(C.ACTIVITY(a,4), ["all", "bursting"])
        error("CompareNeurons: activity '%s': the organoids must be all or bursting", C.ACTIVITY(a,1));
    end
end
if numel(unique(C.ACTIVITY(:,1))) < size(C.ACTIVITY, 1)
    error("CompareNeurons: two activity measures have the same name");
end
check_list(C.COMPOSITION, ["MeanWaveformScore1", "MeanWaveformScore2", "MeanFiringPatternScore"], "COMPOSITION");
check_list(C.COMPOSITION_MODE, ["separate", "joint"], "COMPOSITION_MODE");
check_list(C.ROLES, ["onset", "participation"], "ROLES");
check_list(C.ROLE_SCORES, ["waveformScore1", "waveformScore2", "firingPatternScore"], "ROLE_SCORES");
check_list(C.ROLE_COVARIATES, ["noiseSimilarity", "logRateOutside"], "ROLE_COVARIATES");
if isempty(C.COMPOSITION) || isempty(C.ROLE_SCORES)
    error("CompareNeurons: COMPOSITION and ROLE_SCORES need at least one score");
end
if isempty(C.OUTPUT)
    C.OUTPUT = here;
end
if isempty(C.NAME) || strlength(C.NAME) == 0
    C.NAME = "comparison";
end
end

function S = run_settings(file)
% the structure C built by a settings file (a script)
C = struct();
run(char(file));
S = C;
end

function check_list(values, allowed, name)
bad = setdiff(values, allowed);
if ~isempty(bad)
    error("CompareNeurons: %s: unknown value(s) %s (allowed: %s)", name, strjoin(bad, ', '), strjoin(allowed, ', '));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   ORGANOIDS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [O, Tab] = find_organoids(C)
% every "_Analyzed.mat" file of the group folders of the tables, with the row of the table it comes from (for its folder and columns)
Tab = cell(numel(C.TABLES), 1);
O = table('Size', [0 7], 'VariableTypes', {'string', 'string', 'string', 'string', 'double', 'double', 'string'}, ...
    'VariableNames', {'Recording', 'Group', 'Project', 'File', 'TableIndex', 'TableRow', 'RecFolder'});
for t=1:numel(C.TABLES)
    if ~isfile(C.TABLES(t))
        error("CompareNeurons: the recordings table %s does not exist", C.TABLES(t));
    end
    T = read_table(C.TABLES(t));
    Tab{t} = T;
    table_dir = string(fileparts(char(C.TABLES(t))));
    if isempty(C.PROJECTS)
        project = table_dir;
    else
        project = C.PROJECTS(t);
    end
    groups = unique(T.group(strlength(T.group) > 0), 'stable');
    for g=groups'
        d = dir(fullfile(project, g, '*_Analyzed.mat'));
        for k=1:numel(d)
            rec = string(extractBefore(d(k).name, '_Analyzed'));
            [row, folder] = match_row(T, g, rec, table_dir);
            O(end+1, :) = {rec, g, project, string(fullfile(d(k).folder, d(k).name)), t, row, folder};   %#ok<AGROW>
        end
    end
end
end

function T = read_table(file)
% recordings table, every column as text (as "run_pipeline")
if endsWith(lower(file), ".csv")
    opts = detectImportOptions(file, 'TextType', 'string', 'VariableNamingRule', 'preserve', 'FileType', 'text', 'Delimiter', ',');
else
    opts = detectImportOptions(file, 'TextType', 'string', 'VariableNamingRule', 'preserve');
end
opts = setvartype(opts, opts.VariableNames, 'string');
T = readtable(file, opts);
if ~all(ismember(["file_name", "group"], string(T.Properties.VariableNames)))
    error("CompareNeurons: the recordings table %s needs the columns file_name and group", file);
end
for v=string(T.Properties.VariableNames)
    x = T.(v);  x(ismissing(x)) = "";  T.(v) = strtrim(x);
end
end

function [row, folder] = match_row(T, group, rec, table_dir)
% row of the table whose file name (and organoid label) starts the name of the analyzed file
row = NaN;  folder = "";  best = 0;
has_label = ismember("organoid_label", string(T.Properties.VariableNames));
has_folder = ismember("folder", string(T.Properties.VariableNames));
for r=find(T.group == group)'
    prefix = T.file_name(r) + "_";
    if has_label && strlength(T.organoid_label(r)) > 0 && startsWith(rec, T.file_name(r) + "_" + T.organoid_label(r) + "_")
        prefix = T.file_name(r) + "_" + T.organoid_label(r) + "_";
    end
    if strlength(T.file_name(r)) > 0 && startsWith(rec, prefix) && strlength(prefix) > best
        best = strlength(prefix);  row = r;
    end
end
if ~isnan(row) && has_folder
    folder = T.folder(row);
    if strlength(folder) > 0 && ~(startsWith(folder, filesep) || ~isempty(regexp(folder, '^[A-Za-z]:', 'once')))
        folder = string(fullfile(table_dir, folder));
    end
end
end

function d = duplicated_values(x)
[~, first] = unique(x, 'stable');
d = true(size(x));  d(first) = false;
end

function [NS, B, source, note] = load_organoid(o, min_units)
% NeuronScore of the recording (saved, or computed here from the baseline and the templates) and its baseline condition
NS = [];  source = "none";  note = "";
v = whos('-file', o.File);
names = {v.name};
L = load(o.File, 'PerCondition');
P = L.PerCondition;
if isfield(P, 'file')
    k = baselineCondition({P.file});
else
    k = 1;
end
B = P(k);
n_units = 0;
i = find(strcmp(names, 'spikeTimestamp'), 1);
if ~isempty(i)
    n_units = v(i).size(2);
end
if ismember('NeuronScore', names)
    NS = load(o.File, 'NeuronScore').NeuronScore;
    if ~isempty(NS)
        source = "saved";
        return
    end
end
if n_units < min_units || ~isfield(B, 'spikeTrain') || isempty(B.spikeTrain)
    note = "too few units";
    return
end
st = load(o.File, 'spikeTimestamp').spikeTimestamp;
templates = [];  unit_type = [];  rate = NaN;
source = "computed here";
try
    [templates, unit_type, rate] = loadKilosortTemplates(char(o.RecFolder), st);
catch err
    source = "computed here, firing pattern only";
    note = "no Kilosort templates (" + string(err.message) + ")";
end
NS = NeuronScores(B.spikeTrain, B.window, B.NetworkBurstDetected, templates, rate, unit_type);
NS.cluster = double([st.cluster])';
NS.condition = B.condition;
end

function R = fixed_reference()
% the fixed axes of "NeuronScores" (without any parameter set)
previous = parameter_set('get');
parameter_set('clear');
c = onCleanup(@() restore_set(previous));
R = NeuronScores(struct('relativeRate', 0, 'lv', 1, 'cv', 1, 'burstIndex', 0), [], [], [], NaN).reference;
end

function NS2 = rescore(NS, axes_file)
% the scores of the saved measures and features on the new axes (the other settings of the scores as they were)
previous = parameter_set('get');
P = struct('NeuronScores', NS.parameters);
P.NeuronScores.AXES_FILE = string(axes_file);
parameter_set('set', P);
c = onCleanup(@() restore_set(previous));
NS2 = NeuronScores(NS.FiringPattern, [], [], NS.WaveformFeature, NaN);
for f=["cluster", "condition"]
    if isfield(NS, f)
        NS2.(f) = NS.(f);
    end
end
end

function restore_set(previous)
if isempty(previous)
    parameter_set('clear');
else
    parameter_set('set', previous);
end
end

function c = cosine(a, b)
c = (a(:)' * b(:)) / (norm(a) * norm(b));
end

function U = units_table(O, NS, Base, rate_offset)
% one row per unit: scores on the new axes, and roles in the primary network bursts of the baseline
U = table();
for i=1:height(O)
    S = NS{i};
    if isempty(S), continue, end
    n = numel(S.waveformScore1);
    N = Base{i}.NetworkBurstDetected;
    n1 = 0;
    if isstruct(N) && isfield(N, 'nBurst') && isscalar(N.nBurst) && ~isnan(N.nBurst)
        n1 = N.nBurst;
    end
    participant = false(n, 1);  rank = nan(n, 1);  participation = nan(n, 1);
    if n1 > 0
        if isfield(N, 'Channels') && ~isempty(N.Channels)
            ch = N.Channels(:,1);
            participant(ch) = true;
            if isfield(Base{i}, 'NetworkBurstStructure') && isstruct(Base{i}.NetworkBurstStructure) && ...
                    isfield(Base{i}.NetworkBurstStructure, 'structure') && size(Base{i}.NetworkBurstStructure.structure, 1) == numel(ch)
                rank(ch) = Base{i}.NetworkBurstStructure.structure(:,1);
            end
        end
        if isfield(N, 'timeStamp') && size(N.timeStamp, 2) == n && size(N.timeStamp, 1) >= n1
            participation = mean(~cellfun(@isempty, N.timeStamp(1:n1, :)), 1)';
        end
    end
    onset = nan(n, 1);
    if max(rank) > 1
        onset = 100 * (rank - 1) / (max(rank) - 1);
    end
    cluster = nan(n, 1);
    if isfield(S, 'cluster'), cluster = S.cluster(:); end
    rate = S.FiringPattern.rateOutside(:);
    U = [U; table(repmat(O.Recording(i), n, 1), repmat(O.Group(i), n, 1), cluster, S.waveformScore1, S.waveformScore2, ...
        S.firingPatternScore, S.noiseSimilarity, S.clean, rate, log10(rate + rate_offset), participant, onset, participation, ...
        repmat(n1, n, 1), 'VariableNames', {'Organoid', 'Group', 'Cluster', 'waveformScore1', 'waveformScore2', 'firingPatternScore', ...
        'noiseSimilarity', 'clean', 'rateOutside', 'logRateOutside', 'participant', 'onset', 'participation', 'PrimaryBursts'})];  %#ok<AGROW>
end
end

function Org = organoid_table(O, NS, Base, source, note)
% one row per organoid: units, composition (mean scores on the new axes)
n = height(O);
v = nan(n, 9);
for i=1:n
    N = Base{i}.NetworkBurstDetected;
    if isstruct(N) && isfield(N, 'nBurst') && isscalar(N.nBurst)
        v(i,9) = N.nBurst;
        if isnan(N.nBurst), v(i,9) = 0; end
    end
    if isempty(NS{i}), continue, end
    c = NS{i}.composition;
    v(i,1:8) = [c.nUnits c.nWaveformScored c.nFiringScored c.meanWaveformScore1 c.meanWaveformScore2 c.meanFiringPatternScore ...
        c.meanNoiseSimilarity log10(c.nUnits)];
end
Org = [table(O.Recording, O.Group, 'VariableNames', {'Organoid', 'Group'}), array2table(v, 'VariableNames', {'Units', ...
    'WaveformScored', 'FiringScored', 'MeanWaveformScore1', 'MeanWaveformScore2', 'MeanFiringPatternScore', 'MeanNoiseSimilarity', ...
    'LogUnits', 'PrimaryBursts'}), table(source, note, 'VariableNames', {'ScoresSource', 'Note'})];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   ACTIVITY, COVARIATES   %%%%%%%%%%%%%%%%%%%%%%%%

function M = load_compiled(O, C)
% the compiled tables of each group (only those named in the settings, and NETbursting1 for the bursting organoids)
need = ["NETbursting1"; extractBefore(C.ACTIVITY(:,2), ".")];
dotted = C.COVARIATES(contains(C.COVARIATES, "."));
need = unique([need; extractBefore(dotted(:), ".")]);
M = containers.Map();
keys = unique(O.Project + filesep + O.Group);
missing = strings(0, 1);
for k=keys'
    file = k + "_Compiled.mat";
    if ~isfile(file)
        missing(end+1) = file;                                              %#ok<AGROW>
        continue
    end
    v = whos('-file', file);
    present = string({v.name});
    absent = setdiff(need, present);
    if ~isempty(absent)
        error("CompareNeurons: %s has no table(s) %s (compile the group again)", file, strjoin(absent, ', '));
    end
    M(char(k)) = load(file, need{:});
end
if ~isempty(missing)
    error("CompareNeurons: compile these groups first (Compile_MEAproject): %s", strjoin(missing, ', '));
end
end

function [value, found] = compiled_value(M, o, B, C, spec)
% value of "<table>.<column>" for the organoid, in the condition of the settings
value = NaN;  found = false;
parts = split(spec, ".");
T = M(char(o.Project + filesep + o.Group)).(parts(1));
if ~ismember(parts(2), string(T.Properties.VariableNames))
    error("CompareNeurons: the table %s has no column %s", parts(1), parts(2));
end
rows = string(T.Organoid) == o.Recording;
if lower(C.CONDITION) == "baseline"
    if isfield(B, 'file')
        rows = rows & string(T.File) == string(B.file);
    end
else
    rows = rows & lower(string(T.Condition)) == lower(C.CONDITION);
end
k = find(rows, 1);
if isempty(k), return, end
found = true;
x = T.(parts(2))(k);
if iscell(x), x = x{1}; end
value = double(x);
end

function [Org, names] = add_activity(Org, O, Base, M, C)
% the activity measures (raw values; the transformation is applied in the models) and whether the organoid has primary network bursts
n = height(O);
names = strings(size(C.ACTIVITY, 1), 1);
bursting = nan(n, 1);
missing = false(n, 1);
for i=1:n
    [b, found] = compiled_value(M, O(i,:), Base{i}, C, "NETbursting1.Number_NETburst_Events");
    missing(i) = ~found;
    if found
        bursting(i) = double(b > 0);                                        % (NaN: bursts rejected by the pipeline = no bursts)
    end
end
Org.Bursting = bursting;
for a=1:size(C.ACTIVITY, 1)
    names(a) = "Activity_" + matlab.lang.makeValidName(C.ACTIVITY(a,1));
    x = nan(n, 1);
    for i=1:n
        x(i) = compiled_value(M, O(i,:), Base{i}, C, C.ACTIVITY(a,2));
    end
    Org.(names(a)) = x;
end
Org.Note(missing) = strip(Org.Note(missing) + "; not in the compiled tables of its group (compile the group again)", "left", ";");
Org.Note = strip(strip(Org.Note, "left", ";"));
end

function [Org, names, kind] = add_covariates(Org, O, Base, M, Tab, list, C, is_random)
% covariates (built in, compiled columns, columns of the recordings tables) or random-effect columns (always categories)
if nargin < 8, is_random = false; end
built_in = ["MeanNoiseSimilarity", "LogUnits", "PrimaryBursts", "Units"];
names = strings(numel(list), 1);  kind = strings(numel(list), 1);
for j=1:numel(list)
    name = list(j);
    if ~is_random && ismember(name, built_in)
        names(j) = name;  kind(j) = "number";
        continue
    end
    if ~is_random && contains(name, ".")
        names(j) = "Covariate_" + matlab.lang.makeValidName(strrep(name, ".", "_"));
        x = nan(height(O), 1);
        for i=1:height(O)
            x(i) = compiled_value(M, O(i,:), Base{i}, C, name);
        end
        Org.(names(j)) = x;  kind(j) = "number";
        continue
    end
    text = strings(height(O), 1);                                           % a column of the recordings tables
    found = false;
    for i=1:height(O)
        T = Tab{O.TableIndex(i)};
        if ismember(name, string(T.Properties.VariableNames))
            found = true;
            if ~isnan(O.TableRow(i))
                text(i) = T.(name)(O.TableRow(i));
            end
        end
    end
    if ~found
        error("CompareNeurons: '%s' is not a built-in covariate, a compiled column (<table>.<column>) or a column of the recordings tables", name);
    end
    prefix = "Covariate_";
    if is_random, prefix = "Random_"; end
    names(j) = prefix + matlab.lang.makeValidName(name);
    number = str2double(text);
    if ~is_random && all(isnan(number) == (text == ""))
        Org.(names(j)) = number;  kind(j) = "number";
    else
        text(text == "") = missing;
        Org.(names(j)) = text;  kind(j) = "category";
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   MODELS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = compare_composition(Org, groups, label, C, random_names, opts)
% do the groups differ in the mean scores of their organoids? score ~ group (+ MeanNoiseSimilarity for waveform score 2) (+ random)
R = table();
for s=C.COMPOSITION
    cov = strings(0, 1);
    if s == "MeanWaveformScore2"
        cov = "MeanNoiseSimilarity";                                        % (user's choice: waveform score 2 always with the similarity)
    end
    [T, random, note] = model_rows(Org, groups, s, [cov; random_names], random_names, opts);
    base = struct('Comparison', label, 'Score', s, 'Covariates', joined(cov, ', '));
    if numel(categories(T.Group)) < 2
        R = [R; group_rows(base, groups, [], [], T, s, "", "fewer than 2 groups with organoids")];   %#ok<AGROW>
        continue
    end
    formula = s + " ~ 1 + Group" + terms(cov) + random_terms(random);
    [m, msg] = fit_model(T, formula, "linear", ~isempty(random));
    R = [R; group_rows(base, groups, m, opts, T, s, model_name("linear", random), join_notes([note, msg]))];   %#ok<AGROW>
end
end

function R = group_rows(base, groups, m, opts, T, s, model, note)
% one row per group against the reference group
R = table();
for g=groups(2:end)
    r = base;
    r.Reference = groups(1);  r.Group = g;
    r.Organoids = height(T);
    r.ReferenceMean = mean(T.(s)(T.Group == groups(1)));
    r.GroupMean = mean(T.(s)(T.Group == g));
    [r.Difference, r.Lower, r.Upper, r.p] = coefficient(m, "Group_" + g, opts);
    r.DifferenceInSD = r.Difference / std(T.(s));
    r.Model = model;  r.Note = note;
    R = [R; struct2table(r, 'AsArray', true)];                              %#ok<AGROW>
end
end

function R = compare_activity(Org, groups, label, spec, name, C, covariate_names, covariate_kind, random_names, opts)
% does the composition explain the activity? activity ~ group + covariates (+ composition) (+ random), for each composition score
% (or all of them, COMPOSITION_MODE = "joint")
R = table();
kind = "linear";
if spec(3) == "yes/no", kind = "logistic"; end
if C.COMPOSITION_MODE == "joint"
    sets = {C.COMPOSITION};
else
    sets = num2cell(C.COMPOSITION);
end
Org.Activity = Org.(name);
switch spec(3)
    case "log"
        Org.Activity(~(Org.Activity > 0)) = NaN;
        Org.Activity = log10(Org.Activity);
    case "yes/no"
        Org.Activity = double(Org.Activity > 0);                            % (no value: no network bursts)
        Org.Activity(isnan(Org.(name)) & isnan(Org.Bursting)) = NaN;
end
if spec(4) == "bursting"
    Org = Org(Org.Bursting == 1, :);
end
for k=1:numel(sets)
    comp = reshape(string(sets{k}), [], 1);
    cov = covariate_names(:);
    if any(comp == "MeanWaveformScore2") && ~any(cov == "MeanNoiseSimilarity")
        cov = [cov; "MeanNoiseSimilarity"];                                 %#ok<AGROW> (user's choice)
    end
    [T, random, note] = model_rows(Org, groups, "Activity", [cov; comp; random_names], random_names, opts);
    base = struct('Comparison', label, 'Activity', spec(1), 'Column', spec(2), 'Transformation', spec(3), 'Organoids', spec(4), ...
        'N', height(T), 'Composition', joined(comp, ' + '), 'Covariates', joined(cov, ', '), 'Model', model_name(kind, random));
    if numel(categories(T.Group)) < 2 || height(T) < 4
        R = [R; activity_rows(base, groups, comp, [], [], [], [], opts, T, "fewer than 2 groups or 4 organoids with all the values")];   %#ok<AGROW>
        continue
    end
    n_terms = (numel(categories(T.Group)) - 1) + numel(comp) + term_count(T, cov);
    if height(T) < opts.MIN_ORGANOIDS_PER_TERM * n_terms
        note = join_notes([note, sprintf("%d organoids for %d terms (rule of thumb: %d per term)", height(T), n_terms, opts.MIN_ORGANOIDS_PER_TERM)]);
    end
    rt = random_terms(random);
    [m0, n0] = fit_model(T, "Activity ~ 1 + Group" + terms(cov) + rt, kind, ~isempty(random));
    [m1, n1] = fit_model(T, "Activity ~ 1 + Group" + terms(cov) + terms(comp) + rt, kind, ~isempty(random));
    [g0, ng0] = fit_model(T, "Activity ~ 1" + terms(cov) + rt, kind, ~isempty(random));
    [g1, ng1] = fit_model(T, "Activity ~ 1" + terms(cov) + terms(comp) + rt, kind, ~isempty(random));
    tests = [lrt(g0, m0), lrt(g1, m1), lrt(m0, m1)];
    R = [R; activity_rows(base, groups, comp, m0, m1, tests, std_of(T, comp), opts, T, join_notes([note, n0, n1, ng0, ng1]))];   %#ok<AGROW>
end
end

function R = activity_rows(base, groups, comp, m0, m1, tests, sds, opts, T, note)
% one row per group against the reference group and per composition score
R = table();
if isempty(tests), tests = [NaN NaN NaN]; end
if isempty(sds), sds = nan(size(comp)); end
for g=groups(2:end)
    for j=1:numel(comp)
        r = base;
        r.Group = g + " vs " + groups(1);
        [r.GroupEffect, r.GroupLower, r.GroupUpper, r.GroupP] = coefficient(m0, "Group_" + g, opts);
        [r.GroupEffectAdjusted, r.GroupLowerAdjusted, r.GroupUpperAdjusted, r.GroupPAdjusted] = coefficient(m1, "Group_" + g, opts);
        r.GroupChangePercent = 100 * (r.GroupEffect - r.GroupEffectAdjusted) / r.GroupEffect;
        r.GroupTestP = tests(1);  r.GroupTestPAdjusted = tests(2);
        r.CompositionTerm = comp(j);
        [r.CompositionEffect, r.CompositionLower, r.CompositionUpper, r.CompositionP] = coefficient(m1, comp(j), opts);
        r.CompositionEffectPerSD = r.CompositionEffect * sds(j);
        r.CompositionTestP = tests(3);
        r.Note = note;
        R = [R; struct2table(r, 'AsArray', true)];                          %#ok<AGROW>
    end
end
end

function R = compare_roles(U, groups, label, role, C, opts)
% within the organoids, do the scores relate to the role? normal scores within each organoid; mixed model with a random intercept and
% random slopes of the scores for each organoid: all the organoids, each group, and the difference between the groups
R = table();
vars = [role, C.ROLE_SCORES, C.ROLE_COVARIATES];
ok = all(~isnan(U{:, vars}), 2);
U = U(ok, :);
[o, ~, oi] = unique(U.Organoid);
enough = accumarray(oi, 1) >= C.MIN_ROLE_UNITS;
U = U(enough(oi), :);
if height(U) == 0
    base = struct('Comparison', label, 'Role', role, 'Covariates', joined(C.ROLE_COVARIATES, ', '));
    R = role_rows(base, "all groups", C.ROLE_SCORES, [], U, opts, sprintf("no organoid with %d primary network bursts and %d units with all the values", ...
        C.MIN_BURSTS, C.MIN_ROLE_UNITS));
    return
end
for v=vars
    U.(v) = normal_scores(U.(v), U.Organoid);
end
U.Group = categorical(U.Group, groups);
U.Group = removecats(U.Group);
U.Organoid = categorical(U.Organoid);
base = struct('Comparison', label, 'Role', role, 'Covariates', joined(C.ROLE_COVARIATES, ', '));
fixed = role + " ~ 1" + terms(C.ROLE_SCORES) + terms(C.ROLE_COVARIATES);
% all the organoids of the comparison
[m, note] = fit_role(U, fixed, C.ROLE_SCORES);
R = [R; role_rows(base, "all groups", C.ROLE_SCORES, m, U, opts, note)];
% each group
for g=groups
    Ug = U(U.Group == g, :);
    if numel(unique(Ug.Organoid)) < opts.MIN_GROUP_ORGANOIDS
        R = [R; role_rows(base, g, C.ROLE_SCORES, [], Ug, opts, sprintf("fewer than %d organoids", opts.MIN_GROUP_ORGANOIDS))];   %#ok<AGROW>
        continue
    end
    Ug.Organoid = removecats(Ug.Organoid);
    [m, note] = fit_role(Ug, fixed, C.ROLE_SCORES);
    R = [R; role_rows(base, g, C.ROLE_SCORES, m, Ug, opts, note)];          %#ok<AGROW>
end
% difference between the groups (interaction of the group with each score)
if numel(categories(U.Group)) >= 2
    inter = role + " ~ 1 + Group" + terms(C.ROLE_SCORES) + terms(C.ROLE_COVARIATES) + strjoin(" + Group:" + C.ROLE_SCORES, "");
    [m, note] = fit_role(U, inter, C.ROLE_SCORES);
    for g=groups(2:end)
        r = role_rows(base, g + " - " + groups(1), C.ROLE_SCORES, [], U(U.Group == g | U.Group == groups(1), :), opts, note);
        for j=1:numel(C.ROLE_SCORES)
            [r.Estimate(j), r.Lower(j), r.Upper(j), r.p(j)] = coefficient(m, ["Group_" + g + ":" + C.ROLE_SCORES(j), ...
                C.ROLE_SCORES(j) + ":Group_" + g], opts);
        end
        R = [R; r];                                                         %#ok<AGROW>
    end
end
end

function [m, note] = fit_role(U, fixed, scores)
% random intercept and random slopes of the scores (independent); a random intercept only if that model cannot be fitted
[m, note] = fit_model(U, fixed + " + (1" + terms(scores) + " | Organoid)", "linear", true, "Diagonal");
if isempty(m) || contains(note, "warning")
    [m2, note2] = fit_model(U, fixed + " + (1 | Organoid)", "linear", true);
    if ~isempty(m2)
        m = m2;
        note = join_notes(["random intercept only (" + note + ")", note2]);
    end
end
end

function R = role_rows(base, group, scores, m, U, opts, note)
R = table();
for s=scores
    r = base;
    r.Groups = group;
    r.Score = s;
    r.Organoids = numel(unique(U.Organoid));
    r.Units = height(U);
    [r.Estimate, r.Lower, r.Upper, r.p] = coefficient(m, s, opts);
    r.Note = note;
    R = [R; struct2table(r, 'AsArray', true)];                              %#ok<AGROW>
end
end

function [T, random, note] = model_rows(Org, groups, y, vars, random_names, opts)
% the organoids of the comparison with all the values of the model; group as categories (the first one is the reference); random terms
% kept when they have at least 2 levels
note = strings(0, 1);
T = Org(ismember(Org.Group, groups), :);
vars = unique([y; vars(:)], 'stable');
keep = true(height(T), 1);
for v=vars'
    x = T.(v);
    keep = keep & ~ismissing(x);
end
T = T(keep, :);
T.Group = removecats(categorical(T.Group, groups));
for v=setdiff(vars', [y, random_names(:)'])
    if isstring(T.(v))
        T.(v) = categorical(T.(v));                                         % covariates given as text: categories
    end
end
random = strings(0, 1);
for v=random_names(:)'
    T.(v) = categorical(T.(v));
    levels = numel(categories(removecats(T.(v))));
    if levels >= 2
        random(end+1) = v;                                                  %#ok<AGROW>
        if levels < opts.MIN_RANDOM_LEVELS
            note(end+1) = sprintf("%s has %d levels", extractAfter(v, "Random_"), levels);   %#ok<AGROW>
        end
    else
        note(end+1) = sprintf("%s left out (one level)", extractAfter(v, "Random_"));   %#ok<AGROW>
    end
end
note = join_notes(note);
end

function [m, note] = fit_model(T, formula, kind, has_random, pattern)
% linear (mixed) or logistic (mixed) model; warnings and errors become a note
m = [];  note = "";
if nargin < 5, pattern = ""; end
w = warning;
warning('off', 'backtrace');
c = onCleanup(@() warning(w));
try
    if kind == "logistic" && has_random
        txt = evalc('m = fitglme(T, formula, ''Distribution'', ''Binomial'', ''FitMethod'', ''Laplace'');');
    elseif kind == "logistic"
        txt = evalc('m = fitglm(T, formula, ''Distribution'', ''binomial'');');
    elseif has_random && strlength(pattern) > 0
        txt = evalc('m = fitlme(T, formula, ''FitMethod'', ''ML'', ''CovariancePattern'', {char(pattern)});');
    elseif has_random
        txt = evalc('m = fitlme(T, formula, ''FitMethod'', ''ML'');');
    else
        txt = evalc('m = fitlm(T, formula);');
    end
    found = regexp(txt, 'Warning: ([^\n]*)', 'tokens');
    if ~isempty(found)
        note = "warning: " + strjoin(unique(string(cellfun(@(x) x{1}, found, 'UniformOutput', false))), " / ");
    end
catch err
    m = [];
    note = "not fitted: " + string(err.message);
end
end

function [est, lo, hi, p] = coefficient(m, names, opts)
% estimate, confidence interval and p of the coefficient (the first of the names found); NaN if the model or the coefficient is missing
est = NaN;  lo = NaN;  hi = NaN;  p = NaN;
if isempty(m), return, end
cn = string(m.CoefficientNames);
j = find(ismember(cn, string(names)), 1);
if isempty(j), return, end
if isa(m, 'LinearMixedModel')                                               % degrees of freedom of Satterthwaite (they follow the number
    [~, ~, S] = fixedEffects(m, 'DFMethod', 'satterthwaite', 'Alpha', 1 - opts.LEVEL);   % of organoids, not of units)
    est = S.Estimate(j);  p = S.pValue(j);  lo = S.Lower(j);  hi = S.Upper(j);
    if isnan(p) || isnan(lo)                                                % (not defined, e.g., a random-slope variance at 0: the
        df = max(min(grouping_levels(m)) - 1, 1);                           % levels of the grouping variable - 1, conservative)
        t = tinv(1 - (1 - opts.LEVEL) / 2, df);
        p = 2 * tcdf(-abs(est / S.SE(j)), df);
        lo = est - t * S.SE(j);  hi = est + t * S.SE(j);
    end
    return
end
est = m.Coefficients.Estimate(j);
p = m.Coefficients.pValue(j);
if isa(m, 'GeneralizedLinearMixedModel')
    ci = coefCI(m, 'Alpha', 1 - opts.LEVEL);
else
    ci = coefCI(m, 1 - opts.LEVEL);
end
lo = ci(j,1);  hi = ci(j,2);
end

function n = grouping_levels(m)
% number of levels of each grouping variable of a mixed model
gv = m.Formula.GroupingVariableNames;
names = strings(0, 1);
for q=1:numel(gv)
    names = [names; reshape(string(gv{q}), [], 1)];                         %#ok<AGROW>
end
V = m.Variables;
n = arrayfun(@(v) numel(unique(V.(v))), names);
end

function p = lrt(m0, m1)
% likelihood-ratio test of the nested models
p = NaN;
if isempty(m0) || isempty(m1), return, end
df = numel(m1.CoefficientNames) - numel(m0.CoefficientNames);
if df < 1, return, end
p = chi2cdf(max(2 * (m1.LogLikelihood - m0.LogLikelihood), 0), df, 'upper');
end

function s = terms(names)
names = names(strlength(names) > 0);
s = "";
if ~isempty(names)
    s = strjoin(" + " + reshape(names, 1, []), "");
end
end

function s = random_terms(random)
s = "";
if ~isempty(random)
    s = strjoin(" + (1 | " + reshape(random, 1, []) + ")", "");
end
end

function k = term_count(T, names)
% number of coefficients of the covariates (categories: levels - 1)
k = 0;
for v=reshape(names, 1, [])
    if iscategorical(T.(v)) || isstring(T.(v))
        k = k + numel(unique(T.(v))) - 1;
    else
        k = k + 1;
    end
end
end

function s = std_of(T, names)
s = arrayfun(@(v) std(T.(v)), names);
end

function name = model_name(kind, random)
if kind == "logistic", name = "logistic"; else, name = "linear"; end
if isempty(random)
    name = name + " regression";
else
    name = name + " mixed model (random intercepts: " + joined(extractAfter(random, "Random_"), ", ") + ")";
end
end

function s = joined(x, sep)
% strjoin, "" for an empty list
x = string(x);
s = "";
if ~isempty(x)
    s = strjoin(reshape(x, 1, []), sep);
end
end

function s = comparison_labels(C)
labels = strings(1, numel(C.COMPARISONS));
for k=1:numel(C.COMPARISONS)
    labels(k) = strjoin(C.COMPARISONS{k}, " vs ");
end
s = strjoin(labels, "; ");
end

function s = join_notes(notes)
notes = string(notes);
notes = notes(strlength(notes) > 0);
s = strjoin(notes, "; ");
if isempty(notes), s = ""; end
end

function z = normal_scores(x, organoid)
% ranks within each organoid turned into normal quantiles
z = nan(size(x));
[~, ~, g] = unique(organoid);
for k=1:max(g)
    i = find(g == k);
    z(i) = norminv((tiedrank(x(i)) - 0.5) / numel(i));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   EXCEL   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_excel(Res, file)
C = Res.settings;
A = Res.axes;
setting = ["Tables"; "Condition"; "Groups for the axes"; "Comparisons"; "Composition"; "Covariates"; "Random intercepts"; ...
    "Roles"; "Role scores"; "Role covariates"];
value = [joined(C.TABLES, '; '); C.CONDITION; joined(unique(Res.Organoids.Group, 'stable'), ', '); ...
    comparison_labels(C); C.COMPOSITION_MODE + ": " + joined(C.COMPOSITION, ', '); ...
    joined(C.COVARIATES, ', '); joined(C.RANDOM, ', '); joined(C.ROLES, ', '); joined(C.ROLE_SCORES, ', '); joined(C.ROLE_COVARIATES, ', ')];
for a=1:size(C.ACTIVITY, 1)
    setting(end+1) = "Activity " + a;                                       %#ok<AGROW>
    value(end+1) = strjoin(C.ACTIVITY(a,:), " | ");                         %#ok<AGROW>
end
writetable(table(setting, value, 'VariableNames', {'Setting', 'Value'}), file, 'Sheet', 'Settings', 'WriteMode', 'replacefile');
F = A.fit;
item = ["Version"; "Organoids (waveform / firing pattern)"; "Units (waveform / firing pattern)"; "Variance of the axes (%)"; ...
    "Stability between halves of the organoids (median |cos|)"; "Lowest |cos| between halves"; "Agreement with the fixed axes (|cos|)"];
value = [A.version; sprintf("%d / %d", F.nOrganoidsWaveform, F.nOrganoidsFiring); sprintf("%d / %d", F.nUnitsWaveform, F.nUnitsFiring); ...
    sprintf("%.0f, %.0f, %.0f", 100 * F.eigenvaluesWaveform(1:2) / numel(A.WF_MEAN), 100 * F.eigenvaluesFiring(1) / numel(A.FP_MEAN)); ...
    sprintf("%.3f, %.3f, %.3f", median(F.stabilityWaveform), median(F.stabilityFiring)); ...
    sprintf("%.3f, %.3f, %.3f", min(F.stabilityWaveform), min(F.stabilityFiring)); sprintf("%.3f, %.3f, %.3f", Res.axesAgreement)];
writetable(table(item, value, 'VariableNames', {'Axes', 'Value'}), file, 'Sheet', 'Axes');
writetable(F.organoids, file, 'Sheet', 'Axes', 'Range', 'D1');
for sheet=["Organoids", "Composition", "Activity", "Roles", "Units"]
    if height(Res.(sheet)) > 0
        writetable(Res.(sheet), file, 'Sheet', sheet);
    end
end
end
