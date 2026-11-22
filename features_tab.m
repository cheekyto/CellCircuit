function api = features_tab(tab, ctx)

% TAB "6. Features plotting" OF THE WINDOW "CellCircuit": plots of any measure of the compiled tables, comparing the groups chosen
% ("plotFeatures") (new, November 2026; user's choices)
%   1. groups: their compiled files ("<group>_Compiled.mat", written by stage 8 next to the "<group>_Compiled.xlsx" file): the groups of
%      the results folder of tab 1, or files / folders added; which groups are plotted (ticked), their names on the plots and their order
%   2. table (a sheet of the Excel file: Firing, Bursting, Network Bursting ..., cNEs ..., Neuron scores, Neuron composition) and measures
%      (one panel each)
%   3. recording conditions (one block of groups per condition)
%   4. plot: points with a bar (mean +/- SEM or SD), a violin plot, or a box plot (median, 25th-75th percentiles, whiskers at the
%      percentiles chosen, mean); options
% A preview follows every choice (the first 12 measures); Draw opens every page full size; Save figure... (.png, .pdf -several pages-, .svg,
% .fig); Export values... writes the values drawn and their summary (n, mean, SD, SEM, median, percentiles) to an Excel file
% Called by "CellCircuit" with the tab and ctx, a structure of functions of the window: results_folder() (folder of the results of tab 1),
% message(txt, title), guard(fn) (callback that reports its errors), with_dialog(fn) (file dialogs in front), open_file(f), test (true: no
% dialogs). Returns the functions and controls used by the tests

STYLES = ["bar", "violin", "box"];
WHISKERS = {[5 95], [10 90], [2.5 97.5], [0 100]};

files = strings(0, 1);      % compiled files read (the order they were added: each group keeps its color)
D = [];                     % output "data" of "plotFeatures"
order = zeros(0, 1);        % rows of the group table: index of each group in D (the order of the plot)
hold_preview = false;       % (while the lists are filled, no preview)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Layout

g = uigridlayout(tab, [1 2], 'ColumnWidth', {450, '1x'});
left = uigridlayout(g, [5 1], 'RowHeight', {'fit', 'fit', 'fit', 'fit', 'fit'}, 'Padding', [0 0 0 0], 'Scrollable', 'on');
right = uigridlayout(g, [3 1], 'RowHeight', {'1x', 'fit', 30}, 'Padding', [0 0 0 0]);     % (preview, information, buttons)

uilabel(left, 'WordWrap', 'on', 'Text', ['Plots any measure of the compiled tables (the sheets of "<group>_Compiled.xlsx", written by stage 8) ' ...
    'for the groups you choose. Each point is one organoid in one condition (one unit in "Neuron scores"); values that are not numbers ' ...
    '(e.g., network-burst measures without enough events) are left out. Your files are not changed.']);

% 1. groups
p1 = uipanel(left, 'Title', '1. Groups (compiled files)');
g1 = uigridlayout(p1, [3 1], 'RowHeight', {150, 22, 'fit'});
group_table = uitable(g1, 'Data', empty_rows(), 'ColumnEditable', [true true false], 'ColumnWidth', {40, 130, 'auto'}, ...
    'SelectionType', 'row', 'CellEditCallback', ctx.guard(@(~,~) on_table_edit()), ...
    'Tooltip', 'Tick the groups to plot; the names can be changed (as they appear on the plots); the order of the rows is the order of the plot');
b1 = uigridlayout(g1, [1 6], 'ColumnWidth', {95, 75, 85, 60, 40, 45}, 'Padding', [0 0 0 0], 'ColumnSpacing', 4);
uibutton(b1, 'Text', 'Groups of tab 1', 'ButtonPushedFcn', ctx.guard(@(~,~) add_tab1_groups()), ...
    'Tooltip', 'The compiled files of the results folder of tab 1 (the folder of the recordings table, or the project folder)');
uibutton(b1, 'Text', 'Add files...', 'ButtonPushedFcn', ctx.guard(@(~,~) on_add_files()), 'Tooltip', 'Compiled files of any project ("<group>_Compiled.mat")');
uibutton(b1, 'Text', 'Add folder...', 'ButtonPushedFcn', ctx.guard(@(~,~) on_add_folder()), 'Tooltip', 'Every compiled file of a folder');
uibutton(b1, 'Text', 'Remove', 'ButtonPushedFcn', ctx.guard(@(~,~) remove_selected()), 'Tooltip', 'Remove the selected rows');
uibutton(b1, 'Text', 'Up', 'ButtonPushedFcn', ctx.guard(@(~,~) move_selected(-1)), 'Tooltip', 'Move the selected row up (to the left on the plot)');
uibutton(b1, 'Text', 'Down', 'ButtonPushedFcn', ctx.guard(@(~,~) move_selected(1)), 'Tooltip', 'Move the selected row down (to the right on the plot)');
group_info = uilabel(g1, 'WordWrap', 'on', 'FontColor', [0.3 0.3 0.3], 'Text', 'Add the groups of tab 1, or compiled files.');

% 2. table and measures
p2 = uipanel(left, 'Title', '2. Table and measures');
g2 = uigridlayout(p2, [3 2], 'RowHeight', {22, 170, 22}, 'ColumnWidth', {'1x', '1x'});
sheet_list = uidropdown(g2, 'Items', {'(add groups first)'}, 'ItemsData', {''}, 'ValueChangedFcn', ctx.guard(@(~,~) on_sheet()), ...
    'Tooltip', 'A sheet of the Excel files of the groups');
sheet_list.Layout.Column = [1 2];
measure_list = uilistbox(g2, 'Items', {}, 'Multiselect', 'on', 'ValueChangedFcn', ctx.guard(@(~,~) preview()), ...
    'Tooltip', 'One panel per measure (Ctrl / Cmd + click to choose several; the units come from the Excel sheet)');
measure_list.Layout.Column = [1 2];
uibutton(g2, 'Text', 'All measures', 'ButtonPushedFcn', ctx.guard(@(~,~) set_all_measures(true)));
uibutton(g2, 'Text', 'Only the first', 'ButtonPushedFcn', ctx.guard(@(~,~) set_all_measures(false)));

% 3. conditions
p3 = uipanel(left, 'Title', '3. Recording conditions');
g3 = uigridlayout(p3, [1 1], 'RowHeight', {80});
condition_list = uilistbox(g3, 'Items', {}, 'Multiselect', 'on', 'ValueChangedFcn', ctx.guard(@(~,~) preview()), ...
    'Tooltip', 'Condition and dose, as the compiled tables group the recordings; several conditions: one block of groups each');

% 4. plot
p4 = uipanel(left, 'Title', '4. Plot');
g4 = uigridlayout(p4, [5 2], 'RowHeight', {22, 22, 22, 22, 22}, 'ColumnWidth', {'1x', 170});
uilabel(g4, 'Text', 'Plot:');
style_list = uidropdown(g4, 'Items', {'Points + bar (mean +/- error bar)', 'Points + violin plot', 'Points + box plot (percentiles)'}, ...
    'ItemsData', cellstr(STYLES), 'ValueChangedFcn', ctx.guard(@(~,~) update_options()));
error_label = uilabel(g4, 'Text', 'Error bar:');
error_list = uidropdown(g4, 'Items', {'SEM (standard error)', 'SD (standard deviation)'}, 'ItemsData', {'sem', 'sd'}, ...
    'ValueChangedFcn', ctx.guard(@(~,~) preview()));
whisker_label = uilabel(g4, 'Text', 'Whiskers (box from the 25th to the 75th):');
whisker_list = uidropdown(g4, 'Items', {'5th-95th percentiles', '10th-90th percentiles', '2.5th-97.5th percentiles', 'Minimum-maximum'}, ...
    'ItemsData', {1, 2, 3, 4}, 'ValueChangedFcn', ctx.guard(@(~,~) preview()));
points_check = uicheckbox(g4, 'Text', 'Individual values (points)', 'Value', true, 'ValueChangedFcn', ctx.guard(@(~,~) preview()));
mean_check = uicheckbox(g4, 'Text', 'Mean (diamond)', 'Value', true, 'ValueChangedFcn', ctx.guard(@(~,~) preview()), ...
    'Tooltip', 'Violin and box plots: the mean as a diamond (the median is the white dot or the line)');
n_check = uicheckbox(g4, 'Text', 'Number of values (n)', 'Value', true, 'ValueChangedFcn', ctx.guard(@(~,~) preview()));
column_box = uigridlayout(g4, [1 2], 'ColumnWidth', {'1x', 45}, 'Padding', [0 0 0 0]);
uilabel(column_box, 'Text', 'Panels per row (0 = auto):');
column_field = uieditfield(column_box, 'numeric', 'Value', 0, 'Limits', [0 12], 'RoundFractionalValues', 'on', ...
    'ValueChangedFcn', ctx.guard(@(~,~) preview()));

% preview and buttons
preview_box = uipanel(right, 'Title', 'Preview', 'BackgroundColor', 'w');
info_label = uilabel(right, 'WordWrap', 'on', 'FontColor', [0.4 0.4 0.4], 'Text', 'The preview appears here once groups are added.');
gb = uigridlayout(right, [1 4], 'ColumnWidth', {120, 120, 130, '1x'}, 'Padding', [0 0 0 0]);
uibutton(gb, 'Text', 'Draw', 'FontWeight', 'bold', 'ButtonPushedFcn', ctx.guard(@(~,~) draw()), ...
    'Tooltip', 'Open the plots full size (one window per 12 measures)');
uibutton(gb, 'Text', 'Save figure...', 'ButtonPushedFcn', ctx.guard(@(~,~) on_save()), ...
    'Tooltip', '.png, .pdf (one page per 12 measures, vector), .svg or .fig');
uibutton(gb, 'Text', 'Export values...', 'ButtonPushedFcn', ctx.guard(@(~,~) on_export()), ...
    'Tooltip', 'Excel file: the values drawn (sheet Values) and their summary per group (sheet Summary)');
update_options();

api = struct('add_files', @add_files, 'inputs', @inputs, 'draw', @draw, 'save', @save_figure, 'export', @export_values, ...
    'preview', @preview, 'data', @get_data, 'order', @get_order, ...
    'controls', struct('groups', group_table, 'sheet', sheet_list, 'measures', measure_list, 'conditions', condition_list, ...
    'style', style_list, 'error', error_list, 'whiskers', whisker_list, 'points', points_check, 'mean', mean_check, 'n', n_check, ...
    'columns', column_field, 'preview', preview_box, 'info', info_label, 'group_info', group_info));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Groups

    function add_tab1_groups()
        d = string(ctx.results_folder());
        if strlength(d) == 0 || ~isfolder(d)
            ctx.message('Open or save a recordings table in tab 1 first (its results folder holds the compiled files).', 'Features plotting');
            return
        end
        found = compiled_in(d);
        if isempty(found)
            ctx.message(sprintf('No compiled file ("<group>_Compiled.mat") in %s: run stage 8 (Compile the Excel tables) first.', d), 'Features plotting');
            return
        end
        add_files(found);
    end

    function on_add_files()
        [f, d] = ctx.with_dialog(@() uigetfile('*_Compiled.mat', 'Compiled files of the groups', 'MultiSelect', 'on'));
        if isequal(f, 0), return, end
        add_files(string(fullfile(d, f)));
    end

    function on_add_folder()
        d = ctx.with_dialog(@() uigetdir('', 'Folder with compiled files ("<group>_Compiled.mat")'));
        if isequal(d, 0), return, end
        found = compiled_in(string(d));
        if isempty(found)
            ctx.message(sprintf('No compiled file ("<group>_Compiled.mat") in %s.', d), 'Features plotting');
            return
        end
        add_files(found);
    end

    function add_files(new_files)
        % adds compiled files (the ones already listed are skipped) and reads them all again
        new_files = string(new_files(:));
        new_files = new_files(~ismember(new_files, files));
        if isempty(new_files)
            return
        end
        read_files([files; new_files], [order; numel(files) + (1:numel(new_files))']);
    end

    function remove_selected()
        sel = selected_rows();
        if isempty(sel), return, end
        keep_rows = setdiff(1:numel(order), sel);
        kept = order(keep_rows);
        new_files = files(sort(kept));                                      % (the order they were added)
        [~, new_order] = ismember(files(kept), new_files);
        read_files(new_files, new_order(:));
    end

    function move_selected(step)
        sel = selected_rows();
        if ~isscalar(sel), return, end
        to = sel + step;
        if to < 1 || to > numel(order), return, end
        Dt = group_table.Data;
        swap = [to sel];
        order([sel to]) = order(swap);
        Dt([sel to], :) = Dt(swap, :);
        group_table.Data = Dt;
        group_table.Selection = to;
        preview();
    end

    function read_files(new_files, new_order)
        % reads the compiled files (plotFeatures "data") and fills the lists, keeping the choices that still apply
        old = group_table.Data;
        if isempty(new_files)
            files = strings(0, 1);  D = [];     order = zeros(0, 1);
            group_table.Data = empty_rows();
            fill_lists();
            return
        end
        [~, Dn] = plotFeatures(new_files, "data");
        files = new_files;  D = Dn;     order = new_order(:);
        rows = empty_rows();
        for r=1:numel(order)
            k = order(r);
            use = true;     name = D.groups(k);
            j = find(old.File == files(k), 1);
            if ~isempty(j)
                use = old.Plot(j);  name = old.Group(j);
            end
            rows = [rows; table(use, name, files(k), 'VariableNames', rows.Properties.VariableNames)]; %#ok<AGROW>
        end
        group_table.Data = rows;
        fill_lists();
    end

    function on_table_edit()
        Dt = group_table.Data;
        empty = strlength(strtrim(Dt.Group)) == 0;
        if any(empty)
            Dt.Group(empty) = D.groups(order(empty));                       % (an empty name: the name of the file)
            group_table.Data = Dt;
        end
        preview();
    end

    function sel = selected_rows()
        sel = group_table.Selection;                                        % (SelectionType 'row': the indices of the rows)
        if ~isempty(sel), sel = unique(sel(:))'; end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lists and options

    function fill_lists()
        hold_preview = true;
        restore = onCleanup(@() release_preview());
        if isempty(D)
            sheet_list.Items = {'(add groups first)'};   sheet_list.ItemsData = {''};
            measure_list.Items = {};    condition_list.Items = {};
            group_info.Text = 'Add the groups of tab 1, or compiled files.';
        else
            old = string(sheet_list.Value);
            names = [D.sheets.name];
            sheet_list.Items = cellstr(names);  sheet_list.ItemsData = cellstr(names);
            if ismember(old, names), sheet_list.Value = char(old); else, sheet_list.Value = char(names(1)); end
            group_info.Text = sprintf('%d groups read: %d tables found.', numel(D.groups), numel(D.sheets));
            fill_sheet();
        end
        clear restore
        preview();
    end

    function on_sheet()
        hold_preview = true;
        restore = onCleanup(@() release_preview());
        fill_sheet();
        clear restore
        preview();
    end

    function fill_sheet()
        % measures and conditions of the table chosen (the measures / conditions chosen before are kept when the table has them)
        S = current_sheet();
        if isempty(S), return, end
        old_m = string(measure_list.Value);     old_c = string(condition_list.Value);
        items = S.measures;
        shown = items;
        has_unit = strlength(S.units) > 0;
        shown(has_unit) = items(has_unit) + "  (" + S.units(has_unit) + ")";
        measure_list.Items = cellstr(shown);    measure_list.ItemsData = cellstr(items);
        keep = intersect(old_m, items, 'stable');
        if isempty(keep), keep = items(1); end
        measure_list.Value = cellstr(keep);
        condition_list.Items = cellstr(S.conditions);    condition_list.ItemsData = cellstr(S.conditions);
        keep = intersect(old_c, S.conditions, 'stable');
        if isempty(keep)
            base = S.conditions(lower(S.conditions) == "baseline");
            if isempty(base), base = S.conditions(1); end
            keep = base;
        end
        condition_list.Value = cellstr(keep);
    end

    function set_all_measures(every)
        S = current_sheet();
        if isempty(S), return, end
        if every
            measure_list.Value = cellstr(S.measures);
        else
            measure_list.Value = cellstr(S.measures(1));
        end
        preview();
    end

    function update_options()
        st = string(style_list.Value);
        error_label.Visible = st == "bar";     error_list.Visible = st == "bar";
        whisker_label.Visible = st == "box";   whisker_list.Visible = st == "box";
        mean_check.Enable = st ~= "bar";
        preview();
    end

    function S = current_sheet()
        S = [];
        if isempty(D), return, end
        k = find([D.sheets.name] == string(sheet_list.Value), 1);
        if ~isempty(k), S = D.sheets(k); end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots

    function [style, in] = inputs()
        % the choices of the tab, as the inputs of "plotFeatures"
        if isempty(D)
            error('Add the groups first (their compiled files).');
        end
        Dt = group_table.Data;
        use = order(Dt.Plot);
        if isempty(use)
            error('Tick at least one group.');
        end
        if isempty(measure_list.Value)
            error('Choose at least one measure.');
        end
        if isempty(condition_list.Value)
            error('Choose at least one recording condition.');
        end
        labels = D.groups;
        labels(order) = string(Dt.Group);
        style = string(style_list.Value);
        in = struct('groups', use(:)', 'labels', labels, 'sheet', string(sheet_list.Value), 'measures', string(measure_list.Value), ...
            'conditions', string(condition_list.Value), 'points', points_check.Value, 'error_bar', string(error_list.Value), ...
            'whiskers', WHISKERS{whisker_list.Value}, 'show_mean', mean_check.Value, 'show_n', n_check.Value, 'columns', column_field.Value);
    end

    function preview()
        if hold_preview || ~isvalid(preview_box)
            return
        end
        delete(preview_box.Children);
        info_label.FontColor = [0.4 0.4 0.4];
        if isempty(D)
            info_label.Text = 'The preview appears here once groups are added.';
            return
        end
        try
            [style, in] = inputs();
            in.parent = preview_box;
            [~, res] = plotFeatures(D, style, in);
        catch ME
            info_label.Text = ['Not drawn: ' ME.message];
            info_label.FontColor = [0.8 0 0];
            return
        end
        nm = numel(in.measures);
        txt = sprintf('%d groups, %d condition(s), %d measure(s); %d values drawn.', numel(in.groups), numel(in.conditions), nm, height(res.values));
        if nm > 12
            txt = sprintf('%s The preview shows the first 12 measures: Draw and Save show them all (%d pages).', txt, ceil(nm/12));
        end
        info_label.Text = txt;
    end

    function figs = draw()
        [style, in] = inputs();
        figs = plotFeatures(D, style, in);
    end

    function on_save()
        [f, d] = ctx.with_dialog(@() uiputfile({'*.png', 'PNG image'; '*.pdf', 'PDF (vector; several pages)'; '*.svg', 'SVG (vector)'; ...
            '*.fig', 'MATLAB figure'}, 'Save the figure', 'features.png'));
        if isequal(f, 0), return, end
        file = save_figure(string(fullfile(d, f)));
        ctx.message(sprintf('Saved: %s', file), 'Features plotting');
    end

    function file = save_figure(file)
        [style, in] = inputs();
        in.visible = false;     in.output = string(file);
        [figs, res] = plotFeatures(D, style, in);
        close(figs);
        file = res.file;
    end

    function on_export()
        [f, d] = ctx.with_dialog(@() uiputfile('*.xlsx', 'Export the values drawn', 'features_values.xlsx'));
        if isequal(f, 0), return, end
        file = export_values(string(fullfile(d, f)));
        ctx.message(sprintf('Saved: %s (sheets Values and Summary)', file), 'Features plotting');
    end

    function file = export_values(file)
        % the values drawn and their summary (sheets "Values" and "Summary")
        [style, in] = inputs();
        in.visible = false;
        [figs, res] = plotFeatures(D, style, in);
        close(figs);
        if isfile(file), delete(file); end
        writetable(res.values, file, 'Sheet', 'Values');
        writetable(res.summary, file, 'Sheet', 'Summary');
    end

    function release_preview()
        hold_preview = false;
    end

    function out = get_data()
        out = D;
    end

    function out = get_order()
        out = order;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = empty_rows()

T = table(false(0, 1), strings(0, 1), strings(0, 1), 'VariableNames', {'Plot', 'Group', 'File'});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function found = compiled_in(d)

% compiled files of a folder ("<group>_Compiled.mat")

L = dir(fullfile(d, '*_Compiled.mat'));
found = strings(0, 1);
for k=1:numel(L)
    found(end+1, 1) = string(fullfile(L(k).folder, L(k).name)); %#ok<AGROW>
end

end
