function api = compare_tab(tab, ctx)

% TAB "4. Neuron composition" OF THE WINDOW "CellCircuit": comparison of groups of organoids by their neurons ("CompareNeurons") (new, November 2026;
% named "4. Comparison" until 22 November 2026, renamed by the user's choice)
% The tab builds the settings of "CompareNeurons" (the same as a copy of "comparison_settings.m") and runs it:
%   1. recordings tables (the table of tab 1 and/or others; the groups of a table are in the folder of its results, as for the pipeline)
%   2. groups found (organoids = "_Analyzed.mat" files of <project>/<group>; compiled = <project>/<group>_Compiled.mat), and the groups used
%      to fit the axes of the neuron scores (at least 20 organoids; groups that only help the fitting can be ticked without being compared)
%   3. comparisons (2 or more groups; the first one is the reference)
%   4. activity measures (default list of "CompareNeurons"; any column of the compiled tables can be added)
%   5. models (composition scores, covariates, random intercepts, condition) and 6. roles in the network bursts
%   7. run, save the settings as a file (to run them again with CompareNeurons("<file>")), open a settings file, open the results
% Called by "CellCircuit" with the tab and ctx, a structure of functions of the window: table_file() (file of the table of tab 1, saved
% first if needed; "" if none), results_folder() (folder of the results of tab 1), message(txt, title), guard(fn) (callback that reports
% its errors), with_dialog(fn) (file dialogs in front), open_file(f), test (true: no dialogs)
% Returns the functions and controls used by the tests

D = CompareNeurons("defaults");
STANDARD = ["include", "folder", "file_name", "group", "to_discard", "source", "organoid_label", "picture", "picture_organoid", ...
    "picture_3D_side", "picture_n_organoids", "crop_ellipse", "crop_excel", "file_name_3Brain"];   % columns of the pipeline (not covariates)
BUILT_IN = ["MeanNoiseSimilarity", "LogUnits", "PrimaryBursts"];
SCORE_TEXT = ["Waveform score 1", "Waveform score 2", "Firing-pattern score"];
SCORE_TIP = ["Mean waveform score 1 of the organoid (+ broader, slower spikes)", ...
    "Mean waveform score 2 (+ larger second peak; always with the similarity to the QC noise shape)", ...
    "Mean firing-pattern score (+ higher relative rate and CV, lower LV and burst index)"];
ID_COLUMNS = ["Organoid", "File", "Condition", "dose", "ttx_qc", "T_start", "T_end", "Note"];

tables = strings(0, 1);  projects = strings(0, 1);      % recordings tables and the folder of their groups
comparisons = {};                                       % each one: groups (the first one is the reference)
last_results = "";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Layout

g = uigridlayout(tab, [2 2], 'RowHeight', {'1x', 150}, 'ColumnWidth', {470, '1x'});     % (the run, at the bottom, always visible)
left = uigridlayout(g, [4 1], 'RowHeight', {'fit', 'fit', 'fit', 'fit'}, 'Padding', [0 0 0 0], 'Scrollable', 'on');
right = uigridlayout(g, [3 1], 'RowHeight', {'fit', 'fit', 'fit'}, 'Padding', [0 0 0 0], 'Scrollable', 'on');

uilabel(left, 'WordWrap', 'on', 'Text', ['Compares groups of organoids by the neurons they are made of. The groups must have gone through ' ...
    'the pipeline (stage 8 compiles them). The axes of the neuron scores are fitted again on the groups ticked in "Use for the axes" ' ...
    '(at least 20 organoids; groups that are not compared can help), then only the comparisons listed are run. Your files are not changed.']);

% 1. recordings tables
p1 = uipanel(left, 'Title', '1. Recordings tables (their groups)');
g1 = uigridlayout(p1, [2 1], 'RowHeight', {64, 22});
table_list = uilistbox(g1, 'Items', {}, 'Tooltip', 'Each table, and the folder of its groups and compiled files');
b1 = uigridlayout(g1, [1 4], 'ColumnWidth', {130, 95, 90, '1x'}, 'Padding', [0 0 0 0]);
uibutton(b1, 'Text', 'Add the table of tab 1', 'ButtonPushedFcn', ctx.guard(@(~,~) add_tab1_table()), ...
    'Tooltip', 'The recordings table of tab 1 (saved first), with the results folder of tab 2');
uibutton(b1, 'Text', 'Add a table...', 'ButtonPushedFcn', ctx.guard(@(~,~) on_add_table()), ...
    'Tooltip', 'Another recordings table (.xlsx or .csv); its groups are in its folder');
uibutton(b1, 'Text', 'Remove table', 'ButtonPushedFcn', ctx.guard(@(~,~) on_remove_table()));
uibutton(b1, 'Text', 'Find the groups', 'FontWeight', 'bold', 'ButtonPushedFcn', ctx.guard(@(~,~) find_groups()), ...
    'Tooltip', 'List the groups of the tables, with their organoids');

% 2. groups
p2 = uipanel(left, 'Title', '2. Groups');
g2 = uigridlayout(p2, [2 1], 'RowHeight', {130, 'fit'});
group_table = uitable(g2, 'Data', empty_groups(), 'ColumnEditable', [false false false true], 'ColumnWidth', {150, 80, 75, 125}, ...
    'CellEditCallback', ctx.guard(@(~,~) update_counts()));
count_label = uilabel(g2, 'WordWrap', 'on', 'FontWeight', 'bold', 'Text', 'Add a table and press "Find the groups".');

% 3. comparisons
p3 = uipanel(left, 'Title', '3. Comparisons (only these are run)');
g3 = uigridlayout(p3, [5 2], 'RowHeight', {22, 22, 52, 22, 70}, 'ColumnWidth', {'1x', '1x'});
uilabel(g3, 'Text', 'Reference group:');
uilabel(g3, 'Text', 'Compared with (one or more):');
reference_list = uidropdown(g3, 'Items', {}, 'Tooltip', 'The first group of the comparison (the others are compared with it)');
reference_list.Layout.Row = 2;  reference_list.Layout.Column = 1;
others_list = uilistbox(g3, 'Items', {}, 'Multiselect', 'on', 'Value', {}, 'Tooltip', 'Hold Cmd or Ctrl to select several');
others_list.Layout.Row = [2 3];  others_list.Layout.Column = 2;
add_button = uibutton(g3, 'Text', 'Add comparison', 'ButtonPushedFcn', ctx.guard(@(~,~) on_add_comparison()));
add_button.Layout.Row = 4;  add_button.Layout.Column = 1;
remove_button = uibutton(g3, 'Text', 'Remove comparison', 'ButtonPushedFcn', ctx.guard(@(~,~) on_remove_comparison()));
remove_button.Layout.Row = 4;  remove_button.Layout.Column = 2;
comparison_list = uilistbox(g3, 'Items', {}, 'Tooltip', 'The comparisons that will be run');
comparison_list.Layout.Row = 5;  comparison_list.Layout.Column = [1 2];

% 4. activity measures
p4 = uipanel(right, 'Title', '4. Activity of each organoid (does the composition explain it?)');
g4 = uigridlayout(p4, [3 1], 'RowHeight', {130, 22, 'fit'});
activity_table = uitable(g4, 'Data', activity_data(D.ACTIVITY), 'ColumnEditable', true, 'ColumnWidth', {160, 245, 118, 85}, ...
    'SelectionType', 'row', 'Multiselect', 'on');
b4 = uigridlayout(g4, [1 4], 'ColumnWidth', {'1x', 60, 110, 80}, 'Padding', [0 0 0 0]);
column_list = uidropdown(b4, 'Items', {'(press "Find the groups" to list the compiled columns)'}, ...
    'Tooltip', 'Columns of the compiled tables of the groups');
uibutton(b4, 'Text', 'Add', 'ButtonPushedFcn', ctx.guard(@(~,~) on_add_measure()), 'Tooltip', 'Add the column chosen on the left as a measure');
uibutton(b4, 'Text', 'Remove measure', 'ButtonPushedFcn', ctx.guard(@(~,~) on_remove_measures()));
uibutton(b4, 'Text', 'Default list', 'ButtonPushedFcn', ctx.guard(@(~,~) set_activity(D.ACTIVITY)));
uilabel(g4, 'WordWrap', 'on', 'FontColor', [0.4 0.4 0.4], 'Text', ['Column: <table>.<column> of the compiled tables (Firing, Bursting, ' ...
    'NETbursting0 = pooled, NETbursting1 = primary, NETbursting2 = all clusters, cNEs). Transformation: none, log (log10), yes/no ' ...
    '(value > 0, logistic model). Organoids: all, or bursting (only those with primary network bursts).']);

% 5. models
p5 = uipanel(right, 'Title', '5. Models: activity ~ group + composition + covariates (+ random intercepts)');
g5 = uigridlayout(p5, [4 3], 'RowHeight', {22, 22, 22, 66}, 'ColumnWidth', {'1x', '1x', '1x'});
composition_checks = gobjects(1, 3);
for k=1:3
    composition_checks(k) = uicheckbox(g5, 'Text', SCORE_TEXT(k), 'Value', true, 'Tooltip', SCORE_TIP(k));
end
mode_list = uidropdown(g5, 'Items', {'One model per score', 'All the scores in one model'}, 'ItemsData', {'separate', 'joint'});
condition_label = uilabel(g5, 'Text', 'Condition of the activity:', 'HorizontalAlignment', 'right');
condition_field = uieditfield(g5, 'text', 'Value', 'baseline', 'Tooltip', ['"baseline" (found by its file name, as the pipeline) or the ' ...
    'label of a condition (column Condition of the compiled tables)']);
uilabel(g5, 'Text', 'Covariates:');
uilabel(g5, 'Text', 'Random intercepts:', 'Tooltip', 'Columns you added to the recordings tables');
uilabel(g5, 'Text', '');
covariate_list = uilistbox(g5, 'Items', cellstr(BUILT_IN), 'Multiselect', 'on', 'Value', cellstr(D.COVARIATES), 'Tooltip', ...
    ['MeanNoiseSimilarity (always added with waveform score 2), LogUnits, PrimaryBursts, the extra columns of the recordings tables, ' ...
    'or a compiled column']);
random_list = uilistbox(g5, 'Items', {}, 'Multiselect', 'on', 'Value', {}, 'Tooltip', ['Columns you added to the recordings tables, e.g. ' ...
    '"batch" (none: ordinary regression)']);
uilabel(g5, 'WordWrap', 'on', 'FontColor', [0.4 0.4 0.4], 'Text', ['A note is written when a model has fewer than 10 organoids ' ...
    'per term. The same choices apply to every measure.']);

% 6. roles
p6 = uipanel(right, 'Title', '6. Roles of the units in the network bursts of the baseline');
g6 = uigridlayout(p6, [3 4], 'RowHeight', {22, 22, 22}, 'ColumnWidth', {'1x', '1x', '1x', '1x'});
role_checks = [uicheckbox(g6, 'Text', 'Onset (who starts)', 'Value', true), uicheckbox(g6, 'Text', 'Participation', 'Value', true)];
uilabel(g6, 'Text', 'Minimum bursts:', 'HorizontalAlignment', 'right');
min_bursts_field = uieditfield(g6, 'numeric', 'Value', D.MIN_BURSTS, 'Limits', [1 Inf], 'RoundFractionalValues', 'on');
role_score_checks = [uicheckbox(g6, 'Text', 'Waveform score 1', 'Value', true), uicheckbox(g6, 'Text', 'Waveform score 2', 'Value', true), ...
    uicheckbox(g6, 'Text', 'Firing-pattern score', 'Value', true)];
uilabel(g6, 'Text', '');
role_cov_checks = [uicheckbox(g6, 'Text', 'Noise similarity', 'Value', true, 'Tooltip', 'Covariate: similarity of the template to the QC noise shape'), ...
    uicheckbox(g6, 'Text', 'Rate outside bursts', 'Value', true, 'Tooltip', 'Covariate: log10 of the firing rate outside the network bursts')];
uilabel(g6, 'Text', 'Minimum units:', 'HorizontalAlignment', 'right');
min_units_field = uieditfield(g6, 'numeric', 'Value', D.MIN_ROLE_UNITS, 'Limits', [1 Inf], 'RoundFractionalValues', 'on');

% 7. run
p7 = uipanel(g, 'Title', '7. Run the comparison (its own button: Run at the bottom of the window runs the pipeline)');
p7.Layout.Row = 2;  p7.Layout.Column = [1 2];
g7 = uigridlayout(p7, [1 2], 'ColumnWidth', {520, '1x'});
c7 = uigridlayout(g7, [3 4], 'RowHeight', {22, 22, 22}, 'ColumnWidth', {130, '1x', 75, 150}, 'Padding', [0 0 0 0]);
uilabel(c7, 'Text', 'Results folder:');
output_field = uieditfield(c7, 'text', 'Placeholder', 'the folder of the first table');
uibutton(c7, 'Text', 'Browse...', 'ButtonPushedFcn', ctx.guard(@(~,~) on_browse_output()));
name_field = uieditfield(c7, 'text', 'Value', char(D.NAME), 'Tooltip', 'Name of the results: <name>_<date and time>.xlsx and .mat');
run_button = uibutton(c7, 'Text', 'Run comparison', 'FontWeight', 'bold', 'BackgroundColor', [0.20 0.45 0.75], 'FontColor', [1 1 1], ...
    'ButtonPushedFcn', ctx.guard(@(~,~) run_comparison()));
run_button.Layout.Row = [2 3];  run_button.Layout.Column = 1;
b = uibutton(c7, 'Text', 'Save settings...', 'ButtonPushedFcn', ctx.guard(@(~,~) on_save_settings()), ...
    'Tooltip', 'Save these settings as a file, to run them again: CompareNeurons("<file>")');
b.Layout.Row = 2;  b.Layout.Column = 2;
b = uibutton(c7, 'Text', 'Open settings...', 'ButtonPushedFcn', ctx.guard(@(~,~) on_open_settings()), ...
    'Tooltip', 'Open a settings file (a copy of comparison_settings.m)');
b.Layout.Row = 3;  b.Layout.Column = 2;
b = uibutton(c7, 'Text', 'Open the Excel results', 'ButtonPushedFcn', ctx.guard(@(~,~) on_open_results(true)));
b.Layout.Row = 2;  b.Layout.Column = [3 4];
b = uibutton(c7, 'Text', 'Open the results folder', 'ButtonPushedFcn', ctx.guard(@(~,~) on_open_results(false)));
b.Layout.Row = 3;  b.Layout.Column = [3 4];
log_area = uitextarea(g7, 'Editable', 'off', 'Value', {'The messages of the comparison appear here.'});

api = struct('add_table', @add_table, 'find_groups', @find_groups, 'set_axes_groups', @set_axes_groups, 'add_comparison', @add_comparison, ...
    'settings', @build_settings, 'check', @check, 'run', @run_comparison, 'save_settings', @save_settings, 'load_settings', @load_settings, ...
    'controls', struct('tables', table_list, 'groups', group_table, 'count', count_label, 'reference', reference_list, 'others', others_list, ...
    'comparisons', comparison_list, 'activity', activity_table, 'columns', column_list, 'composition', composition_checks, 'mode', mode_list, ...
    'condition', condition_field, 'covariates', covariate_list, 'random', random_list, 'roles', role_checks, 'role_scores', role_score_checks, ...
    'role_covariates', role_cov_checks, 'min_bursts', min_bursts_field, 'min_units', min_units_field, 'output', output_field, ...
    'name', name_field, 'run', run_button, 'log', log_area));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tables and groups

    function add_tab1_table()
        file = ctx.table_file();
        if strlength(file) == 0
            ctx.message('Open or build a recordings table in tab 1 first (and save it).', 'Neuron composition');
            return
        end
        add_table(file, ctx.results_folder());
    end

    function on_add_table()
        [f, d] = ctx.with_dialog(@() uigetfile({'*.xlsx;*.xls;*.csv', 'Recordings table'}, 'Recordings table'));
        if isequal(f, 0), return, end
        add_table(string(fullfile(d, f)), "");
    end

    function add_table(file, project)
        % adds a table and the folder of its groups ("" = the folder of the table)
        file = string(file);
        if nargin < 2 || strlength(string(project)) == 0
            project = string(fileparts(char(file)));
        end
        if ~isfile(file)
            error("The recordings table %s does not exist", file);
        end
        if ~any(tables == file)
            tables(end+1, 1) = file;
            projects(end+1, 1) = string(project);
        end
        table_list.Items = cellstr(tables + "   ->   groups in " + projects);
        if strlength(string(output_field.Value)) == 0
            output_field.Value = char(fileparts(char(tables(1))));
        end
    end

    function on_remove_table()
        k = find(strcmp(table_list.Items, table_list.Value), 1);
        if isempty(k), return, end
        tables(k) = [];  projects(k) = [];
        table_list.Items = cellstr(tables + "   ->   groups in " + projects);
    end

    function find_groups()
        % groups of the tables, their organoids and compiled files; extra columns of the tables (covariates, random intercepts); columns
        % of the compiled tables (activity measures)
        if isempty(tables)
            ctx.message('Add a recordings table first.', 'Neuron composition');
            return
        end
        old = group_table.Data;
        G = empty_groups();
        extra = strings(0, 1);
        first_compiled = "";
        for t=1:numel(tables)
            T = read_table(tables(t));
            extra = union(extra, setdiff(string(T.Properties.VariableNames(:)), STANDARD));
            for grp=unique(T.group(strlength(T.group) > 0), 'stable')'
                n = numel(dir(fullfile(projects(t), grp, '*_Analyzed.mat')));
                compiled = isfile(fullfile(projects(t), grp + "_Compiled.mat"));
                if compiled && strlength(first_compiled) == 0
                    first_compiled = fullfile(projects(t), grp + "_Compiled.mat");
                end
                k = find(G.Group == grp, 1);
                if isempty(k)
                    use = true;
                    j = find(old.Group == grp, 1);
                    if ~isempty(j), use = old.("Use for the axes")(j); end
                    G = [G; table(grp, n, compiled, use & n > 0, 'VariableNames', G.Properties.VariableNames)];   %#ok<AGROW>
                else
                    G.Organoids(k) = G.Organoids(k) + n;
                    G.Compiled(k) = G.Compiled(k) && compiled;
                end
            end
        end
        group_table.Data = G;
        names = cellstr(G.Group(G.Organoids > 0));
        reference_list.Items = names;
        others_list.Items = names;  others_list.Value = {};
        keep_cov = intersect(string(covariate_list.Value), [BUILT_IN(:); extra; compiled_names(first_compiled)]);
        covariate_list.Items = cellstr([BUILT_IN(:); extra; compiled_names(first_compiled)]);
        covariate_list.Value = cellstr(keep_cov);
        keep_random = intersect(string(random_list.Value), extra);
        random_list.Items = cellstr(extra);
        random_list.Value = cellstr(keep_random);
        cols = compiled_names(first_compiled);
        if isempty(cols)
            column_list.Items = {'(no compiled tables found: run stage 8 for these groups)'};
        else
            column_list.Items = cellstr(cols);
        end
        update_counts();
    end

    function names = compiled_names(file)
        % numeric columns of the compiled tables ("<table>.<column>"), without the identification columns
        names = strings(0, 1);
        if strlength(file) == 0, return, end
        v = whos('-file', file);
        present = string({v.name});
        for tname = intersect(["Firing", "Bursting", "NETbursting0", "NETbursting1", "NETbursting2", "cNEs"], present, 'stable')
            T = load(file, tname).(tname);
            for v = string(T.Properties.VariableNames)
                if ~ismember(v, ID_COLUMNS) && isnumeric(T.(v))
                    names(end+1, 1) = tname + "." + v;                      %#ok<AGROW>
                end
            end
        end
    end

    function update_counts()
        G = group_table.Data;
        if isempty(G)
            return
        end
        n = sum(G.Organoids(G.("Use for the axes")));
        txt = sprintf('%d organoids of %d groups ticked for the axes (at least 20 are needed).', n, sum(G.("Use for the axes")));
        if any(~G.Compiled & G.Organoids > 0)
            txt = [txt sprintf(' Not compiled yet: %s.', strjoin(G.Group(~G.Compiled & G.Organoids > 0), ', '))];
        end
        count_label.Text = txt;
        if n < 20
            count_label.FontColor = [0.8 0 0];
        else
            count_label.FontColor = [0 0.45 0];
        end
    end

    function set_axes_groups(groups)
        % (for the tests and the settings files) the groups used for the axes ([] = every group)
        G = group_table.Data;
        if isempty(groups)
            G.("Use for the axes") = G.Organoids > 0;
        else
            G.("Use for the axes") = ismember(G.Group, string(groups));
        end
        group_table.Data = G;
        update_counts();
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparisons and measures

    function on_add_comparison()
        if isempty(reference_list.Items)
            ctx.message('Press "Find the groups" first.', 'Neuron composition');
            return
        end
        add_comparison([string(reference_list.Value), setdiff(string(others_list.Value), string(reference_list.Value), 'stable')]);
    end

    function add_comparison(groups)
        groups = reshape(string(groups), 1, []);
        if numel(groups) < 2
            ctx.message('Choose a reference group and at least one group to compare with it.', 'Neuron composition');
            return
        end
        if ~any(cellfun(@(c) isequal(c, groups), comparisons))
            comparisons{end+1} = groups;
        end
        show_comparisons();
    end

    function on_remove_comparison()
        k = find(strcmp(comparison_list.Items, comparison_list.Value), 1);
        if isempty(k), return, end
        comparisons(k) = [];
        show_comparisons();
    end

    function show_comparisons()
        labels = cellfun(@(c) char(strjoin(c, " vs ")), comparisons, 'UniformOutput', false);
        comparison_list.Items = labels;
    end

    function T = activity_data(A)
        A = string(A);
        if isempty(A), A = strings(0, 4); end
        T = table(A(:,1), A(:,2), categorical(A(:,3), ["none", "log", "yes/no"]), categorical(A(:,4), ["all", "bursting"]), ...
            'VariableNames', {'Name', 'Column', 'Transformation', 'Organoids'});
    end

    function set_activity(A)
        activity_table.Data = activity_data(A);
    end

    function on_add_measure()
        col = string(column_list.Value);
        if ~contains(col, ".")
            ctx.message('Press "Find the groups" first: the list shows the columns of the compiled tables.', 'Activity');
            return
        end
        T = activity_table.Data;
        name = extractAfter(col, ".");
        if any(T.Name == name), name = col; end
        set_activity([string(T.Name), string(T.Column), string(T.Transformation), string(T.Organoids); name, col, "none", "all"]);
    end

    function on_remove_measures()
        sel = activity_table.Selection;
        if isempty(sel)
            ctx.message('Click a row of the measures first.', 'Activity');
            return
        end
        T = activity_table.Data;
        T(unique(sel(:,1)), :) = [];
        activity_table.Data = T;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

    function C = build_settings()
        % the settings of "CompareNeurons" given by the tab
        G = group_table.Data;
        T = activity_table.Data;
        C.TABLES = tables';
        C.PROJECTS = projects';
        C.CONDITION = string(strtrim(condition_field.Value));
        C.FIT_GROUPS = G.Group(G.("Use for the axes"))';
        C.COMPARISONS = comparisons;
        C.ACTIVITY = [strtrim(string(T.Name)), strtrim(string(T.Column)), string(T.Transformation), string(T.Organoids)];
        C.COMPOSITION = D.COMPOSITION([composition_checks.Value]);
        C.COMPOSITION_MODE = string(mode_list.Value);
        C.COVARIATES = reshape(string(covariate_list.Value), 1, []);
        C.RANDOM = reshape(string(random_list.Value), 1, []);
        C.ROLES = D.ROLES([role_checks.Value]);
        C.ROLE_SCORES = D.ROLE_SCORES([role_score_checks.Value]);
        C.ROLE_COVARIATES = D.ROLE_COVARIATES([role_cov_checks.Value]);
        C.MIN_BURSTS = min_bursts_field.Value;
        C.MIN_ROLE_UNITS = min_units_field.Value;
        C.OUTPUT = string(strtrim(output_field.Value));
        C.NAME = string(strtrim(name_field.Value));
    end

    function problems = check()
        % what would stop the comparison
        problems = strings(0, 1);
        C = build_settings();
        G = group_table.Data;
        if isempty(C.TABLES)
            problems(end+1) = "Add a recordings table (1) and press Find the groups.";
        elseif isempty(G)
            problems(end+1) = "Press Find the groups (1).";
        end
        n = sum(G.Organoids(G.("Use for the axes")));
        if ~isempty(G) && n < 20
            problems(end+1) = sprintf("Only %d organoids are ticked for the axes (2): at least 20 are needed.", n);
        end
        if isempty(C.COMPARISONS)
            problems(end+1) = "Add at least one comparison (3).";
        else
            missing = setdiff([C.COMPARISONS{:}], C.FIT_GROUPS);
            if ~isempty(missing)
                problems(end+1) = "The groups compared must be ticked for the axes (2): " + strjoin(missing, ", ") + ".";
            end
            not_compiled = intersect([C.COMPARISONS{:}], G.Group(~G.Compiled));
            if ~isempty(not_compiled)
                problems(end+1) = "Compile these groups first (stage 8): " + strjoin(not_compiled, ", ") + ".";
            end
        end
        if isempty(C.ACTIVITY) || any(C.ACTIVITY(:) == "" | ismissing(C.ACTIVITY(:)))
            problems(end+1) = "Every activity measure (4) needs a name, a column, a transformation and the organoids.";
        end
        if isempty(C.COMPOSITION)
            problems(end+1) = "Tick at least one composition score (5).";
        end
        if isempty(C.ROLE_SCORES) && ~isempty(C.ROLES)
            problems(end+1) = "Tick at least one score for the roles (6), or untick the roles.";
        end
    end

    function on_save_settings()
        [f, d] = ctx.with_dialog(@() uiputfile('*.m', 'Save the settings', 'my_comparison.m'));
        if isequal(f, 0), return, end
        save_settings(string(fullfile(d, f)));
    end

    function save_settings(file)
        % the settings as a file like "comparison_settings.m" (to run them again: CompareNeurons("<file>"))
        C = build_settings();
        L = ["% Settings of ""CompareNeurons"", saved by the window CellCircuit (tab 4) on " + string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
             "% To run them again: CompareNeurons(""" + file + """) (see comparison_settings.m for each setting)";
             "C.TABLES = " + text_list(C.TABLES) + ";";
             "C.PROJECTS = " + text_list(C.PROJECTS) + ";";
             "C.CONDITION = " + text_list(C.CONDITION) + ";";
             "C.FIT_GROUPS = " + text_list(C.FIT_GROUPS) + ";";
             "C.COMPARISONS = {" + strjoin(string(cellfun(@text_list, C.COMPARISONS, 'UniformOutput', false)), ", ") + "};";
             "C.ACTIVITY = ["];
        for a=1:size(C.ACTIVITY, 1)
            L(end+1, 1) = "    " + strjoin(quoted(C.ACTIVITY(a,:)), ", ");          %#ok<AGROW>
        end
        L = [L; "    ];";
             "C.COMPOSITION = " + text_list(C.COMPOSITION) + ";";
             "C.COMPOSITION_MODE = " + text_list(C.COMPOSITION_MODE) + ";";
             "C.COVARIATES = " + text_list(C.COVARIATES) + ";";
             "C.RANDOM = " + text_list(C.RANDOM) + ";";
             "C.ROLES = " + text_list(C.ROLES) + ";";
             "C.ROLE_SCORES = " + text_list(C.ROLE_SCORES) + ";";
             "C.ROLE_COVARIATES = " + text_list(C.ROLE_COVARIATES) + ";";
             "C.MIN_BURSTS = " + C.MIN_BURSTS + ";";
             "C.MIN_ROLE_UNITS = " + C.MIN_ROLE_UNITS + ";";
             "C.OUTPUT = " + text_list(C.OUTPUT) + ";";
             "C.NAME = " + text_list(C.NAME) + ";"];
        writelines(L, file);
        add_log("Settings saved: " + file + " (to run them again: CompareNeurons(""" + file + """))");
    end

    function on_open_settings()
        [f, d] = ctx.with_dialog(@() uigetfile('*.m', 'Settings of CompareNeurons'));
        if isequal(f, 0), return, end
        load_settings(string(fullfile(d, f)));
    end

    function load_settings(file)
        % the tab filled in with a settings file (the defaults for the settings it does not give)
        S = run_settings(file);
        C = D;
        for f=string(fieldnames(S))'
            C.(f) = S.(f);
        end
        tables = strings(0, 1);  projects = strings(0, 1);
        T = reshape(string(C.TABLES), [], 1);
        P = reshape(string(C.PROJECTS), [], 1);
        for t=1:numel(T)
            if isempty(P), add_table(T(t), ""); else, add_table(T(t), P(t)); end
        end
        find_groups();
        set_axes_groups(C.FIT_GROUPS);
        comparisons = cellfun(@(c) reshape(string(c), 1, []), C.COMPARISONS, 'UniformOutput', false);
        show_comparisons();
        set_activity(C.ACTIVITY);
        for k=1:3
            composition_checks(k).Value = ismember(D.COMPOSITION(k), string(C.COMPOSITION));
        end
        mode_list.Value = char(string(C.COMPOSITION_MODE));
        condition_field.Value = char(string(C.CONDITION));
        covariate_list.Items = cellstr(union(string(covariate_list.Items), string(C.COVARIATES), 'stable'));
        covariate_list.Value = cellstr(string(C.COVARIATES));
        random_list.Items = cellstr(union(string(random_list.Items), string(C.RANDOM), 'stable'));
        random_list.Value = cellstr(string(C.RANDOM));
        for k=1:2
            role_checks(k).Value = ismember(D.ROLES(k), string(C.ROLES));
            role_cov_checks(k).Value = ismember(D.ROLE_COVARIATES(k), string(C.ROLE_COVARIATES));
        end
        for k=1:3
            role_score_checks(k).Value = ismember(D.ROLE_SCORES(k), string(C.ROLE_SCORES));
        end
        min_bursts_field.Value = C.MIN_BURSTS;
        min_units_field.Value = C.MIN_ROLE_UNITS;
        if ~isempty(C.OUTPUT), output_field.Value = char(string(C.OUTPUT)); end
        name_field.Value = char(string(C.NAME));
        add_log("Settings opened: " + file);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run, results

    function ok = run_comparison()
        ok = false;
        problems = check();
        if ~isempty(problems)
            ctx.message(sprintf('Before running:\n- %s', strjoin(problems, '\n- ')), 'The comparison cannot start');
            return
        end
        C = build_settings();
        if isempty(C.OUTPUT) || strlength(C.OUTPUT) == 0
            C.OUTPUT = string(fileparts(char(C.TABLES(1))));
        end
        run_button.Enable = 'off';
        restore = onCleanup(@() set_enable(run_button));
        add_log("Comparison started: " + strjoin(string(cellfun(@(c) strjoin(c, " vs "), C.COMPARISONS, 'UniformOutput', false)), "; "));
        drawnow;
        try
            [txt, R] = evalc('CompareNeurons(C)');
            fprintf('%s', txt);                                             % (also in the Command Window)
            for line = splitlines(string(strtrim(txt)))'
                add_log(line);
            end
            last_results = R.files.excel;
            add_log("Done: " + R.files.excel);
            ok = true;
        catch ME
            add_log("The comparison stopped: " + ME.message);
            ctx.message(sprintf('The comparison stopped:\n%s', ME.message), 'Neuron composition');
        end
    end

    function on_open_results(excel)
        if strlength(last_results) == 0
            ctx.message('Run a comparison first (the results of this window are opened here).', 'Results');
            return
        end
        if excel
            ctx.open_file(last_results);
        else
            ctx.open_file(string(fileparts(char(last_results))));
        end
    end

    function on_browse_output()
        d = ctx.with_dialog(@() uigetdir(char(output_field.Value), 'Results folder of the comparison'));
        if isequal(d, 0), return, end
        output_field.Value = d;
    end

    function add_log(txt)
        v = log_area.Value;
        if numel(v) == 1 && startsWith(v{1}, 'The messages of the comparison')
            v = {};
        end
        v{end+1} = sprintf('%s  %s', char(datetime('now', 'Format', 'HH:mm:ss')), char(txt));
        log_area.Value = v;
        scroll(log_area, 'bottom');
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   LOCAL FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = empty_groups()
G = table(strings(0, 1), zeros(0, 1), false(0, 1), false(0, 1), 'VariableNames', {'Group', 'Organoids', 'Compiled', 'Use for the axes'});
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
if ~ismember("group", string(T.Properties.VariableNames))
    error("The recordings table %s has no column group", file);
end
T.group(ismissing(T.group)) = "";
T.group = strtrim(T.group);
end

function S = run_settings(file)
% the structure C built by a settings file (a script)
C = struct();
run(char(file));
S = C;
end

function s = text_list(x)
% MATLAB text of a list of texts: "a" or ["a", "b"] ([] if empty)
x = reshape(string(x), 1, []);
if isempty(x)
    s = "[]";
elseif isscalar(x)
    s = quoted(x);
else
    s = "[" + strjoin(quoted(x), ", ") + "]";
end
end

function q = quoted(x)
q = """" + replace(string(x), """", """""") + """";
end

function set_enable(b)
if isvalid(b)
    b.Enable = 'on';
end
end
