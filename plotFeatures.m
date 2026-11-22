function [fig, result] = plotFeatures(data, plot_type, inputs)

% PLOTS OF THE COMPILED MEASURES OF GROUPS OF ORGANOIDS (new, November 2026; tab "6. Features plotting" of the window "CellCircuit")
% Compares groups on any measure of the compiled tables of the pipeline ("Compile_MEAproject": "<group>_Compiled.mat", whose tables are
% the sheets of "<group>_Compiled.xlsx"): one panel per measure, the groups side by side (one block of groups per recording condition when
% several conditions are chosen). Each point is one row of the table: an organoid in a condition (a unit for the sheet "Neuron scores")
%   "bar"       the values and a bar at their mean, with an error bar (SEM or SD)
%   "violin"    the values and their distribution (Gaussian kernel density, cut at the smallest and largest value), with the median (white
%               dot), the 25th-75th percentiles (thick line) and the mean (diamond)
%   "box"       the values, a box from the 25th to the 75th percentile, the median (line), the mean (diamond) and whiskers at the percentiles
%               chosen (5th and 95th by default; [0 100] = minimum and maximum)
%   "data"      no plot: reads the compiled files once (groups, tables, measures, conditions, units of the measures), to draw several plots
%               without reading them again
% Values that are not numbers (NaN, e.g., the network-burst measures of a recording with too few events) are left out and counted (missing)
% (user's choices, 21 November 2026) individual points with a bar (mean +/- error bar), a violin plot, or a box plot (mean, median,
% percentiles); the groups chosen by the user; any measure of the compiled tables
%
% Examples:
%   plotFeatures(["A_Compiled.mat", "B_Compiled.mat"], "bar", struct('sheet', "Firing", 'measures', "Firing_Rate"))
%   D = plotFeatures(files, "data");    plotFeatures(D, "box", struct('sheet', "Bursting", 'conditions', ["Baseline", "+Gabazine 20uM"]))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   data:           compiled files ("<group>_Compiled.mat"; string array), or the structure returned by "data"
%   plot_type:      "bar" (default), "violin", "box" or "data"
%   inputs:         (optional) structure that overwrites any of the choices below and of the PARAMETERS

groups = [];                    % groups to plot, by name or number (in the order given); [] = every group of "data"
labels = [];                    % names of the groups on the plots (one per group of "data"); [] = the names of the files without "_Compiled"
sheet = "Firing";               % table: its sheet in the Excel file ("Firing", "Bursting", "Network Bursting Pooled", "Network Bursting Primary",
                                % "Network Bursting ALL clusters", "cNEs", "cNEs vs baseline", "cNEs vs previous", "Neuron scores",
                                % "Neuron composition"), or its name in the .mat file (e.g., "NETbursting0")
measures = [];                  % measures to plot (columns of the table); [] = every measure of the table
conditions = "Baseline";        % recording conditions, by their labels in the tables (condition and dose, e.g., "+Gabazine 20uM"); [] or "all" =
                                % every condition found (tables without conditions: all their rows)
points = true;                  % individual values (one point per row of the table)
error_bar = "sem";              % (bar) "sem" (standard error of the mean) or "sd" (standard deviation)
whiskers = [5 95];              % (box) percentiles of the ends of the whiskers; [0 100] = minimum and maximum
show_mean = true;               % (violin, box) the mean as a diamond
show_n = true;                  % number of values of each group, under it
columns = 0;                    % panels per row; 0 = automatic
per_page = 12;                  % most panels in one figure; more measures make more figures (pages of the PDF)
output = "";                    % file to save: .png, .pdf (several pages if needed), .svg or .fig; "" = not saved
visible = true;                 % =false, the figure is created without being shown (e.g., to save it)
parent = [];                    % container to draw in (e.g., the preview of the window "CellCircuit"), instead of a new figure; only the
                                % first page is drawn there

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   fig:            the figures (one per page), or "parent" when it is given ([] for "data")
%   result:         structure with what was drawn:
%       .values     table, one row per value drawn: Group, Condition, Organoid (or the first column of the table), File, Measure, Value
%       .summary    table, one row per measure, condition and group: N, N_missing, Mean, SD, SEM, Median, P25, P75, Whisker_low, Whisker_high
%                   (the percentiles of "whiskers"), Min, Max
%       .groups, .conditions, .measures (drawn), .sheet, .n_pages, .file (file saved, "" if none)
%       "data"      the compiled files: .files, .groups (names), .labels, .tables{g}.<table name> (the tables of each group), .sheets
%                   (one element per table found: .name -sheet-, .var -name in the .mat file-, .measures, .units -of each measure, from
%                   the first row of the Excel sheet when the .xlsx file is next to the .mat file; "" otherwise-, .conditions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TABLES = ["Firing", "Bursting", "NETbursting0", "NETbursting1", "NETbursting2", "cNEs", "cNEsVsBaseline", "cNEsVsPrevious", "Neurons", ...
    "NeuronComposition"];                                                  % tables of the compiled .mat file ("Compile_MEAproject") ...
SHEETS = ["Firing", "Bursting", "Network Bursting Pooled", "Network Bursting Primary", "Network Bursting ALL clusters", "cNEs", ...
    "cNEs vs baseline", "cNEs vs previous", "Neuron scores", "Neuron composition"];   % ... and their sheets in the Excel file
NOT_MEASURES = ["ttx_qc", "T_start", "T_end"];                              % numeric columns that are not measures (TTX flag, window)
COLORS = [0 114 178; 230 159 0; 0 158 115; 213 94 0; 204 121 167; 86 180 233; 90 90 90; 160 130 0]/255;   % colors of the groups, in
                                % this order (Okabe-Ito, with grey and a dark yellow instead of its pale yellow: distinct for color-blind readers);
                                % each group keeps its color (its place in "data") when others are left out or the order changes
MARKERS = 'osd^v';              % markers of the points: groups 1-8 circles, the next 8 squares, ...
INK = [0.15 0.15 0.15];         % error bars, whiskers, medians
BAR_ALPHA = 0.35;               % transparency of the bars, violins and boxes (the points stay visible on them)
HALF_WIDTH = 0.36;              % half width of a bar, violin or box (the groups are 1 apart)
JITTER = 0.22;                  % (bar, box) largest horizontal spread of the points from the center of their group
CONDITION_GAP = 0.9;            % space between two blocks of conditions, in group widths
POINT_SIZE = 22;                % size of the points (scatter)

if nargin < 2 || strlength(string(plot_type)) == 0, plot_type = "bar"; end
if nargin < 3, inputs = struct(); end
plot_type = lower(string(plot_type));
if ~ismember(plot_type, ["bar", "violin", "box", "data"])
    error('plotFeatures: the plot must be "bar", "violin", "box" or "data" (not "%s")', plot_type);
end
override_inputs(inputs);                    % overwrite the choices / PARAMETERS above with the ones provided as input (if any)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compiled files

if isstruct(data)
    D = data;
else
    D = read_compiled(string(data), TABLES, SHEETS, NOT_MEASURES);
end
if ~isempty(labels)
    labels = string(labels);
    if numel(labels) ~= numel(D.groups)
        error('plotFeatures: give one label per group (%d groups, %d labels)', numel(D.groups), numel(labels));
    end
    D.labels = labels(:);
end
fig = [];
if plot_type == "data"
    result = D;
    return
end

% table, groups, measures and conditions chosen
s = find([D.sheets.name] == string(sheet) | [D.sheets.var] == string(sheet), 1);
if isempty(s)
    error('plotFeatures: the table "%s" is not in these compiled files (tables found: %s)', sheet, strjoin([D.sheets.name], ', '));
end
S = D.sheets(s);
gsel = group_list(D, groups);
if isempty(measures)
    measures = S.measures;
end
measures = string(measures(:))';
unknown = setdiff(measures, S.measures, 'stable');
if ~isempty(unknown)
    error('plotFeatures: not a measure of the table "%s": %s', S.name, strjoin(unknown, ', '));
end
if isempty(measures)
    error('plotFeatures: the table "%s" has no measures', S.name);
end
all_conditions = isempty(conditions) || (isscalar(string(conditions)) && lower(string(conditions)) == "all");
if all_conditions || isequal(S.conditions, "(all)")
    conds = S.conditions;
else
    conds = string(conditions(:))';
    unknown = setdiff(conds, S.conditions, 'stable');
    if ~isempty(unknown)
        error('plotFeatures: condition not found in the table "%s": %s (conditions found: %s)', S.name, strjoin(unknown, ', '), ...
            strjoin(S.conditions, ', '));
    end
end
whiskers = double(whiskers);
if plot_type == "box" && (numel(whiskers) ~= 2 || any(whiskers < 0 | whiskers > 100) || whiskers(1) >= 25 || whiskers(2) <= 75)
    error('plotFeatures: whiskers must be two percentiles, the first below 25 and the second above 75 (e.g., [5 95]; [0 100] = min and max)');
end
error_bar = lower(string(error_bar));
if ~ismember(error_bar, ["sem", "sd"])
    error('plotFeatures: error_bar must be "sem" or "sd" (not "%s")', error_bar);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Values and summary

ng = numel(gsel);   nc = numel(conds);  nm = numel(measures);
vals = cell(nm, ng, nc);    missing = zeros(nm, ng, nc);
V = table();
for g=1:ng
    T = group_table(D, gsel(g), S.var);
    if isempty(T)
        continue
    end
    lab = condition_labels(T);
    ids = row_names(T);
    for c=1:nc
        rows = lab == conds(c);
        for m=1:nm
            col_values = double(T.(measures(m))(rows));
            ok = isfinite(col_values);
            vals{m,g,c} = col_values(ok);
            missing(m,g,c) = sum(~ok);
            k = find(rows);     k = k(ok);
            if ~isempty(k)
                V = [V; table(repmat(D.labels(gsel(g)), numel(k), 1), repmat(conds(c), numel(k), 1), ids.organoid(k), ids.file(k), ...
                    repmat(measures(m), numel(k), 1), col_values(ok), 'VariableNames', {'Group', 'Condition', 'Organoid', 'File', 'Measure', 'Value'})]; %#ok<AGROW>
            end
        end
    end
end
if isempty(V)
    V = table(strings(0,1), strings(0,1), strings(0,1), strings(0,1), strings(0,1), zeros(0,1), ...
        'VariableNames', {'Group', 'Condition', 'Organoid', 'File', 'Measure', 'Value'});
end
summary = summary_table(vals, missing, measures, D.labels(gsel), conds, whiskers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots (one panel per measure, at most "per_page" per figure)

per_page = max(1, round(per_page));
pages = ceil(nm/per_page);
if ~isempty(parent)
    pages = 1;
end
figs = gobjects(0);
for p=1:pages
    idx = (p-1)*per_page + 1 : min(p*per_page, nm);
    n_panels = numel(idx);
    nc_grid = columns;
    if nc_grid <= 0
        nc_grid = min(n_panels, max(1, ceil(sqrt(n_panels*1.3))));
    end
    nr_grid = ceil(n_panels/nc_grid);
    if isempty(parent)
        f = figure('Name', sprintf('Features - %s%s', S.name, page_text(p, pages)), 'Color', 'w', 'Visible', on_off(visible), ...
            'Position', [60 60 min(1800, 120 + nc_grid*(150 + 95*ng*nc)) min(1000, 80 + nr_grid*300)]);
        container = f;
        figs(end+1) = f; %#ok<AGROW>
    else
        container = parent;
    end
    tl = tiledlayout(container, nr_grid, nc_grid, 'TileSpacing', 'compact', 'Padding', 'compact');
    for ip=1:n_panels                                                       % (names not used by the nested functions, which share
        m = idx(ip);                                                        %  the variables of this workspace)
        ax = nexttile(tl);
        draw_measure(ax, squeeze_cells(vals(m,:,:), ng, nc), gsel, D.labels(gsel), conds, plot_type);
        unit = S.units(S.measures == measures(m));
        ylab = measures(m);
        if ~isempty(unit) && strlength(unit) > 0
            ylab = ylab + " (" + unit + ")";
        end
        ylabel(ax, ylab, 'Interpreter', 'none');
    end
    if nc > 1                                                               % (several conditions: the groups in a legend)
        h = gobjects(ng, 1);
        ax1 = nexttile(tl, 1);
        hold(ax1, 'on');
        for g=1:ng
            h(g) = patch(ax1, NaN, NaN, group_color(gsel(g)), 'FaceAlpha', 0.6, 'EdgeColor', group_color(gsel(g)), 'DisplayName', D.labels(gsel(g)));
        end
        lg = legend(ax1, h, 'Interpreter', 'none', 'Orientation', 'horizontal', 'Box', 'off');
        lg.Layout.Tile = 'north';
    end
end

% save
file = "";
if strlength(string(output)) > 0 && isempty(parent)
    file = string(output);
    [~, ~, ext] = fileparts(file);
    for p=1:numel(figs)
        if lower(ext) == ".fig"
            savefig(figs(p), page_file(file, p, pages));
        elseif lower(ext) == ".pdf"                                         % (vector, to edit the plots; the advice to use an
            warn_state = warning('off', 'MATLAB:print:ContentTypeImageSuggested');   %  image instead, for the transparency, is not shown)
            exportgraphics(figs(p), file, 'ContentType', 'vector', 'Append', p > 1);
            warning(warn_state);
        else
            exportgraphics(figs(p), page_file(file, p, pages), 'Resolution', 200);
        end
    end
end
if isempty(parent)
    fig = figs;
else
    fig = parent;
end
result = struct('values', V, 'summary', summary, 'groups', D.labels(gsel), 'conditions', conds, 'measures', measures, 'sheet', S.name, ...
    'n_pages', pages, 'file', file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nested functions (they use the choices and PARAMETERS above)

    function draw_measure(ax, vc, gidx, glabels, conds, style)
        % one panel: groups side by side, one block per condition
        hold(ax, 'on');     box(ax, 'off');
        ax.YGrid = 'on';    ax.GridAlpha = 0.12;    ax.TickDir = 'out';     ax.TickLabelInterpreter = 'none';
        n_g = numel(gidx);  n_c = numel(conds);
        X = zeros(n_g, n_c);
        lo = inf;   hi = -inf;
        for cc=1:n_c
            for gg=1:n_g
                x = (cc-1)*(n_g + CONDITION_GAP) + gg;
                X(gg,cc) = x;
                v = vc{gg,cc};
                if isempty(v)
                    continue
                end
                col = group_color(gidx(gg));
                lo = min(lo, min(v));   hi = max(hi, max(v));
                switch style
                    case "bar"
                        mu = mean(v);
                        patch(ax, x + HALF_WIDTH*[-1 1 1 -1], [0 0 mu mu], col, 'FaceAlpha', BAR_ALPHA, 'EdgeColor', col, 'LineWidth', 1);
                        lo = min(lo, 0);    hi = max(hi, 0);
                        if numel(v) > 1
                            e = std(v);
                            if error_bar == "sem", e = e/sqrt(numel(v)); end
                            line(ax, [x x], mu + [-e e], 'Color', INK, 'LineWidth', 1.2);
                            line(ax, x + 0.12*[-1 1; -1 1]', [mu-e mu-e; mu+e mu+e]', 'Color', INK, 'LineWidth', 1.2);
                            lo = min(lo, mu - e);   hi = max(hi, mu + e);
                        end
                    case "violin"
                        if numel(v) >= 2 && max(v) > min(v)
                            [yy, dd] = density(v);
                            dd = dd/max(dd)*HALF_WIDTH;
                            patch(ax, [x - dd, fliplr(x + dd)], [yy, fliplr(yy)], col, 'FaceAlpha', BAR_ALPHA, 'EdgeColor', col, 'LineWidth', 1);
                        else
                            line(ax, x + HALF_WIDTH*[-1 1], [v(1) v(1)], 'Color', col, 'LineWidth', 1.5);
                        end
                        q = prctile(v, [25 50 75]);
                        line(ax, [x x], q([1 3]), 'Color', INK, 'LineWidth', 3.5);
                        plot(ax, x, q(2), 'o', 'MarkerSize', 5.5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', INK, 'LineWidth', 1);
                    case "box"
                        q = prctile(v, [whiskers(1) 25 50 75 whiskers(2)]);
                        w = HALF_WIDTH*0.8;
                        patch(ax, x + w*[-1 1 1 -1], q([2 2 4 4]), col, 'FaceAlpha', BAR_ALPHA, 'EdgeColor', col, 'LineWidth', 1);
                        line(ax, x + w*[-1 1], q([3 3]), 'Color', INK, 'LineWidth', 2);
                        line(ax, [x x; x x]', [q(1) q(2); q(4) q(5)]', 'Color', INK, 'LineWidth', 1);
                        line(ax, x + 0.4*w*[-1 1; -1 1]', [q(1) q(1); q(5) q(5)]', 'Color', INK, 'LineWidth', 1);
                end
                if points
                    xs = x + jitter(v, style);
                    scatter(ax, xs, v, POINT_SIZE, col, 'filled', 'Marker', MARKERS(1 + mod(floor((gidx(gg)-1)/size(COLORS,1)), numel(MARKERS))), ...
                        'MarkerEdgeColor', 'w', 'LineWidth', 0.4, 'MarkerFaceAlpha', 0.9);
                end
                if show_mean && style ~= "bar"
                    plot(ax, x, mean(v), 'd', 'MarkerSize', 6, 'MarkerFaceColor', INK, 'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
                end
            end
        end
        if ~isfinite(lo)
            lo = 0;     hi = 1;
            text(ax, mean(X(:)), 0.5, 'no values', 'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
        end
        span = hi - lo;
        if ~(span > 0)
            span = max(abs(hi), 1);                                         % (a single value, or all the values equal)
        end
        pad = 0.06*span;
        y0 = lo - pad*(1 + 2*show_n);
        ylim(ax, [y0, hi + pad]);
        xlim(ax, [X(1) - 0.7, X(end) + 0.7]);
        if show_n
            for cc=1:n_c
                for gg=1:n_g
                    text(ax, X(gg,cc), y0, sprintf('n=%d', numel(vc{gg,cc})), 'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'bottom', 'FontSize', 7.5, 'Color', [0.4 0.4 0.4]);
                end
            end
        end
        if n_c == 1
            xticks(ax, X(:,1));     xticklabels(ax, cellstr(glabels));
        else
            xticks(ax, mean(X, 1));     xticklabels(ax, cellstr(conds));
        end
        if n_c == 1 && n_g > 3
            xtickangle(ax, 30);
        end
    end

    function rgb = group_color(gi)
        rgb = COLORS(1 + mod(gi-1, size(COLORS,1)), :);
    end

    function dx = jitter(v, style)
        % horizontal places of the points: spread more where the values are dense (as the violin), in a fixed order (same plot every time)
        n = numel(v);
        if n == 1
            dx = 0;     return
        end
        u = mod((1:n)'*0.6180339887, 1) - 0.5;                              % (golden-ratio sequence: even, not random)
        d = ones(n, 1);
        if n >= 3 && max(v) > min(v)
            d = kde_at(v, v);
            d = d(:)/max(d);
        end
        if style == "violin"
            dx = 2*u.*d*HALF_WIDTH*0.8;
        else
            dx = 2*u.*sqrt(d)*JITTER;
        end
        dx = dx(:)';
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = read_compiled(files, TABLES, SHEETS, NOT_MEASURES)

% the tables of each compiled file, the measures (numeric columns) and conditions of each table, and the units of the measures

if isempty(files)
    error('plotFeatures: give at least one compiled file ("<group>_Compiled.mat")');
end
files = files(:);
D = struct('files', files, 'groups', strings(numel(files), 1), 'labels', strings(numel(files), 1), 'tables', {cell(numel(files), 1)});
for g=1:numel(files)
    if ~isfile(files(g))
        error('plotFeatures: file not found: %s', files(g));
    end
    [~, name] = fileparts(files(g));
    D.groups(g) = regexprep(name, '_Compiled$', '');
    w = whos('-file', files(g));
    present = intersect(TABLES, string({w.name}), 'stable');
    if isempty(present)
        error('plotFeatures: no compiled table in %s (it should come from "Compile_MEAproject")', files(g));
    end
    D.tables{g} = load(files(g), present{:});
end
D.labels = D.groups;
sheets = struct('name', {}, 'var', {}, 'measures', {}, 'units', {}, 'conditions', {});
for t=1:numel(TABLES)
    measures = strings(1, 0);   conds = strings(1, 0);  found = false;
    for g=1:numel(files)
        if ~isfield(D.tables{g}, TABLES(t)) || ~istable(D.tables{g}.(TABLES(t)))
            continue
        end
        found = true;
        T = D.tables{g}.(TABLES(t));
        names = string(T.Properties.VariableNames);
        numeric = arrayfun(@(k) (isnumeric(T.(names(k))) || islogical(T.(names(k)))) && size(T.(names(k)), 2) == 1, 1:numel(names));
        measures = [measures, setdiff(names(numeric), [measures, NOT_MEASURES], 'stable')]; %#ok<AGROW>
        conds = [conds, setdiff(unique(condition_labels(T), 'stable')', conds, 'stable')]; %#ok<AGROW>
    end
    if ~found
        continue
    end
    base = conds(lower(conds) == "baseline");                               % (the baseline first, then the order found)
    conds = [base, conds(lower(conds) ~= "baseline")];
    units = excel_units(files, SHEETS(t), measures);
    sheets(end+1) = struct('name', SHEETS(t), 'var', TABLES(t), 'measures', measures, 'units', units, 'conditions', conds); %#ok<AGROW>
end
D.sheets = sheets;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function units = excel_units(files, sheet, measures)

% units of the measures: the first row of the sheet of the Excel file next to the first compiled file that has one (row 2 = the names of the
% columns); "" when there is no Excel file or the sheet cannot be read

units = strings(size(measures));
for g=1:numel(files)
    xls = regexprep(files(g), '\.mat$', '.xlsx');
    if ~isfile(xls)
        continue
    end
    try
        C = readcell(xls, 'Sheet', sheet, 'Range', '1:2');
    catch
        continue
    end
    if size(C, 1) < 2
        continue
    end
    for k=1:size(C, 2)
        name = C{2,k};
        u = C{1,k};
        if (ischar(name) || isstring(name)) && (ischar(u) || isstring(u))
            units(measures == string(name)) = string(u);
        end
    end
    return
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lab = condition_labels(T)

% label of the condition of each row, as the tables group them: condition and dose ("+Gabazine 20uM"); "(all)" without a column Condition

names = string(T.Properties.VariableNames);
if ~any(names == "Condition")
    lab = repmat("(all)", height(T), 1);
    return
end
lab = to_text(T.Condition);
if any(names == "dose")
    d = to_text(T.dose);
    lab = strtrim(lab + " " + d);
end
lab(strlength(lab) == 0) = "(no condition)";

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ids = row_names(T)

% organoid (or the first column of the table) and file of each row, as text

names = string(T.Properties.VariableNames);
first = names(1);
if any(names == "Organoid"), first = "Organoid"; end
ids.organoid = to_text(T.(first));
if any(names == "File")
    ids.file = to_text(T.File);
else
    ids.file = strings(height(T), 1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = to_text(x)

% a column of a table as a string column (cells of text, numbers or empty values)

if iscell(x)
    s = strings(numel(x), 1);
    for k=1:numel(x)
        v = x{k};
        if isempty(v)
            s(k) = "";
        elseif ischar(v) || isstring(v)
            s(k) = string(v);
        elseif isnumeric(v) || islogical(v)
            s(k) = strjoin(string(v(:)'), " ");
        end
    end
else
    s = string(x(:));
    s(ismissing(s)) = "";
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gsel = group_list(D, groups)

% groups chosen, as indices of D.groups (by name -or label- or number, in the order given)

if isempty(groups)
    gsel = 1:numel(D.groups);
    return
end
if isnumeric(groups)
    gsel = groups(:)';
    if any(gsel < 1 | gsel > numel(D.groups) | gsel ~= round(gsel))
        error('plotFeatures: groups must be numbers from 1 to %d', numel(D.groups));
    end
    return
end
groups = string(groups(:))';
gsel = zeros(1, numel(groups));
for k=1:numel(groups)
    j = find(D.groups == groups(k) | D.labels == groups(k), 1);
    if isempty(j)
        error('plotFeatures: group not found: %s (groups: %s)', groups(k), strjoin(D.groups, ', '));
    end
    gsel(k) = j;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = group_table(D, g, var)

T = [];
if isfield(D.tables{g}, var) && istable(D.tables{g}.(var))
    T = D.tables{g}.(var);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = summary_table(vals, missing, measures, glabels, conds, whiskers)

% one row per measure, condition and group

[nm, ng, nc] = size(vals);
n = nm*ng*nc;
Measure = strings(n, 1);    Condition = strings(n, 1);  Group = strings(n, 1);
st = nan(n, 12);
r = 0;
for m=1:nm
    for c=1:nc
        for g=1:ng
            r = r + 1;
            Measure(r) = measures(m);   Condition(r) = conds(c);    Group(r) = glabels(g);
            v = vals{m,g,c};
            k = numel(v);
            st(r, 1:2) = [k, missing(m,g,c)];
            if k > 0
                q = prctile(v, [50 25 75 whiskers(1) whiskers(2)]);
                sd = NaN;
                if k > 1, sd = std(v); end
                st(r, 3:12) = [mean(v), sd, sd/sqrt(k), q, min(v), max(v)];
            end
        end
    end
end
S = [table(Measure, Condition, Group), array2table(st, 'VariableNames', {'N', 'N_missing', 'Mean', 'SD', 'SEM', 'Median', 'P25', 'P75', ...
    'Whisker_low', 'Whisker_high', 'Min', 'Max'})];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [yy, dd] = density(v)

% Gaussian kernel density of the values, on 120 points from the smallest to the largest value

yy = linspace(min(v), max(v), 120);
dd = kde_at(v, yy);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = kde_at(v, y)

% Gaussian kernel density of the values v at the points y (bandwidth: Silverman's rule of thumb)

v = v(:);   n = numel(v);
s = std(v);
iq = diff(prctile(v, [25 75]))/1.34;
h = 0.9*min(s, iq)*n^(-1/5);
if ~(h > 0)
    h = 0.9*s*n^(-1/5);
end
if ~(h > 0)
    h = max(abs(v(1)), 1)*1e-3;
end
y = y(:)';
d = sum(exp(-0.5*((y - v)/h).^2), 1)/(n*h*sqrt(2*pi));
d = reshape(d, size(y));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vc = squeeze_cells(c, ng, nc)

vc = reshape(c, ng, nc);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = page_text(p, pages)

t = "";
if pages > 1
    t = sprintf(' (%d of %d)', p, pages);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = page_file(file, p, pages)

% file of page p: "<name>_2.png" for the second page, ... (the first page keeps the name)

f = file;
if pages > 1 && p > 1
    [d, n, e] = fileparts(file);
    f = fullfile(d, n + "_" + p + e);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = on_off(tf)

s = 'off';
if tf, s = 'on'; end

end
