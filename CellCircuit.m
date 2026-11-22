function out = CellCircuit(varargin)

% WINDOW TO RUN THE PIPELINE ("run_pipeline") WITHOUT WRITING CODE (new, November 2026)
%
%   CellCircuit                         opens the window
%   CellCircuit('table', "rec.xlsx")    opens the window with a recordings table already loaded
%
% Six tabs, and a bar with the buttons Run and Stop that is always visible:
%   1. Recordings:  the recordings table of "run_pipeline" (one row per recording): open it, add the recordings found in a folder, edit it
%                   (tick boxes "include", group, units to discard, ...) and save it (.xlsx or .csv)
%   2. Settings:    what to run (the whole pipeline, from a stage to the end, or some stages), which parameters (the default file
%                   "pipeline_parameters.m", the set saved for each recording by its last run, or another file), single values changed for
%                   this run, and options (npy-matlab folder, results folder). The values are in two lists: the common adjustments (at most
%                   15; default: "common_parameters.m"; each user can move parameters from and to the full list, and the list is kept in
%                   the MATLAB preferences of that user) and all the parameters (advanced), with a search
%   3. Run:         summary of the run and the equivalent command, progress of each recording and stage (filled in during the run), the
%                   messages of the run, and buttons to open the logs and the Excel results
%   4. Neuron composition: (new, November 2026; "Comparison" until 22 November 2026) comparison of groups of organoids by their neurons ("CompareNeurons"; see "compare_tab"): the
%                   groups of the recordings tables, the groups used to fit the axes of the neuron scores, the comparisons, the activity
%                   measures and the models; it has its own button "Run comparison" (the bar at the bottom runs the pipeline)
%   5. Visualization: (new, November 2026) plots of a recording of the table ("plotMEArecording"): raster plot (+ population firing rate),
%                   bursts of each unit, electrode map, spike waveforms, and a video of the electrode map (MP4); the recording, its
%                   dataset, the time window (whole recording, a condition, or typed), the units and the options are chosen in the tab
%   6. Features plotting: (new, November 2026) plots of any measure of the compiled tables ("<group>_Compiled.xlsx"), comparing the groups
%                   chosen ("plotFeatures"; see "features_tab"): points with a bar (mean +/- SEM or SD), a violin plot or a box plot
%                   (median, percentiles, mean), one panel per measure, one block of groups per recording condition; saved as an image or
%                   a PDF, and the values drawn exported to Excel
% The run is the same as calling "run_pipeline" with the equivalent command shown in the tab Run. Stop: the run stops after the stage that
% is running (a stage cannot be interrupted). The folder of npy-matlab and the last table are remembered (MATLAB preferences "CellCircuit").
%
% (option for the tests) CellCircuit(..., 'test', true): no dialog windows, and fig.UserData.test gives functions to drive the window

opt = struct('table', "", 'test', false);
for iarg=1:2:numel(varargin)                                                % (not "k": the nested functions below use their own "k")
    opt.(char(varargin{iarg})) = varargin{iarg+1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants

code_dir = fileparts(mfilename('fullpath'));
default_param_file = string(fullfile(code_dir, 'pipeline_parameters.m'));
STAGES = ["ShapeSpikeTrainData", "ShapeSpikeTrainData_3Brain", "QualityControl_General", "QualityControl_HyperSynch", "QualityControl_Manual", ...
          "SelectElectrodes", "IntegrateDataSets", "Analyze_spikeTrain", "Compile_MEAproject"];      % order of the pipeline (as "run_pipeline")
END_TO_END = STAGES(STAGES ~= "ShapeSpikeTrainData_3Brain");
STAGE_TEXT = ["1. Read the spike-sorted data (Kilosort)", "1b. Read 3Brain data (3Brain recordings only)", ...
              "2. Quality control: activity, noise and TTX", "3. Quality control: artifactual synchrony between units", ...
              "4. Quality control: units discarded by hand (column to_discard)", "5. Keep the electrodes under the organoid (cropping columns)", ...
              "6. Combine Kilosort and 3Brain units", "7. Analyze firing, bursts and network bursts", "8. Compile the Excel tables of each group"];
STAGE_SHORT = ["Read data", "Read 3Brain", "QC activity", "QC synchrony", "QC by hand", "Cropping", "Combine", "Analyze", "Compile"];
COLUMNS = ["include", "folder", "file_name", "group", "to_discard", "source", "organoid_label", "picture", "picture_organoid", ...
    "picture_3D_side", "picture_n_organoids", "crop_ellipse", "crop_excel", "file_name_3Brain"];
COLOR_DONE = [0.80 0.94 0.80];      COLOR_RUN = [1.00 0.93 0.60];      COLOR_ERROR = [1.00 0.75 0.75];      COLOR_SKIP = [0.90 0.90 0.90];
COLOR_CHANGED = [1.00 0.93 0.60];
MAX_COMMON = 15;                    % (user's choice) at most this number of parameters in the list of common adjustments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% State of the window

common_default = common_parameters();   % the few parameters shown first ("Common adjustments"; see common_parameters.m)
common = load_common_list();        % the list of this user: the default one, or the one changed in the window (MATLAB preferences)
selected_rows = [];                 % rows of the parameter table selected by the user
table_file = "";                    % file of the recordings table ("" = not saved yet)
added_folder = "";                  % (new, November 2026) last folder whose recordings were added: until the table is saved, its folder
                                    % is suggested for the table (and so for the Excel results, saved next to it)
table_dirty = false;                % the table was edited since it was opened or saved
params = [];                        % list of every parameter (see "read_parameter_list")
changes = containers.Map('KeyType', 'char', 'ValueType', 'any');   % values changed for this run: 'Function.PARAMETER' -> text typed
running = false;                    % a run is in progress
stop_requested = false;             % the user asked to stop
close_after_run = false;            % the window closes at the end of the run
run_info = struct('recordings', strings(0), 'groups', strings(0), 'stages', strings(0), 'project', "", 'excel', strings(0));
viz_data = [];                      % (new, November 2026) dataset read in the tab Visualization (output "data" of "plotMEArecording")
viz_photo = struct('photo', "", 'side', "", 'n', 0);  % photo of the chip of the recording chosen (columns picture, picture_3D_side,
                                    % picture_n_organoids of its row), for the outline of the organoids
viz_last = "";                      % last file saved by the tab Visualization (figure or video)
viz_hold = false;                   % (November 2026) no preview while a dataset is read (one preview at the end)
viz_seen_row = [];                  % row selected in tab 1 at the last visit of the tab Visualization (a new selection is proposed there)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Window

fig = uifigure('Name', 'CellCircuit - run the pipeline', 'Position', [80 80 1180 760], 'CloseRequestFcn', @on_close);
main = uigridlayout(fig, [2 1], 'RowHeight', {'1x', 44}, 'Padding', [8 8 8 8]);
tabs = uitabgroup(main);
tab_rec = uitab(tabs, 'Title', '1. Recordings');
tab_set = uitab(tabs, 'Title', '2. Settings');
tab_run = uitab(tabs, 'Title', '3. Run');
tab_cmp = uitab(tabs, 'Title', '4. Neuron composition');                    % (new, November 2026; renamed from "4. Comparison", user's choice)
tab_viz = uitab(tabs, 'Title', '5. Visualization');
tab_feat = uitab(tabs, 'Title', '6. Features plotting');                    % (new, November 2026)

% bar at the bottom (always visible)
bar = uigridlayout(main, [1 5], 'ColumnWidth', {'1x', 90, 120, 90, 70}, 'Padding', [0 0 0 0]);
status_label = uilabel(bar, 'Text', 'Open or build a recordings table (tab 1), check the settings (tab 2), then press Run.', 'FontWeight', 'bold');
uibutton(bar, 'Text', 'Help', 'ButtonPushedFcn', guard(@on_help), 'Tooltip', 'How to use this window');
check_button = uibutton(bar, 'Text', 'Check settings', 'ButtonPushedFcn', guard(@on_check), 'Tooltip', 'Check the table and the settings without running');
run_button = uibutton(bar, 'Text', 'Run', 'FontWeight', 'bold', 'BackgroundColor', [0.30 0.65 0.35], 'FontColor', [1 1 1], ...
    'ButtonPushedFcn', guard(@on_run), 'Tooltip', 'Run the pipeline with the recordings of tab 1 and the settings of tab 2');
stop_button = uibutton(bar, 'Text', 'Stop', 'Enable', 'off', 'ButtonPushedFcn', guard(@on_stop), ...
    'Tooltip', 'Stop after the stage that is running (a stage cannot be interrupted)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tab 1: recordings

g1 = uigridlayout(tab_rec, [5 1], 'RowHeight', {'fit', 'fit', 'fit', '1x', 'fit'});
% (revised, user's request, November 2026: shorter; the details of the cropping columns are in the Help)
uilabel(g1, 'WordWrap', 'on', 'Text', ['One row per recording; tick "include" to run it. "group": recordings compiled together in one Excel ' ...
    'file. "file_name": common name of the files of the recording. Optional: "to_discard" (units to remove by hand, e.g., 3 17 25), ' ...
    '"source" (Kilosort or 3Brain), "file_name_3Brain" (a 3Brain dataset to combine), and the cropping columns ("picture", ' ...
    '"organoid_label"...: keep only the electrodes under each organoid; see Help).']);
b1 = uigridlayout(g1, [1 9], 'ColumnWidth', {105, 170, 105, 115, 90, 160, 60, 75, '1x'}, 'Padding', [0 0 0 0]);
uibutton(b1, 'Text', 'Open table...', 'ButtonPushedFcn', guard(@on_open_table), 'Tooltip', 'Open a recordings table (.xlsx or .csv)');
uibutton(b1, 'Text', 'Add recordings from folder...', 'ButtonPushedFcn', guard(@on_add_folder), ...
    'Tooltip', 'Add one row for each recording found in a folder and its subfolders (Kilosort output or files of the pipeline)');
uibutton(b1, 'Text', 'Add empty row', 'ButtonPushedFcn', guard(@(~,~) add_rows(empty_rows(1))));
uibutton(b1, 'Text', 'Remove selected', 'ButtonPushedFcn', guard(@on_remove_rows), 'Tooltip', 'Remove the selected rows (click a row first)');
group_field = uieditfield(b1, 'text', 'Placeholder', 'group', 'Tooltip', 'Type a group, select rows, and press "Set group of selected rows"');
uibutton(b1, 'Text', 'Set group of selected rows', 'ButtonPushedFcn', guard(@on_group_rows), ...
    'Tooltip', 'The group typed on the left goes to the selected rows (hold Cmd or Ctrl to select several)');
uibutton(b1, 'Text', 'Save', 'ButtonPushedFcn', guard(@(~,~) save_table(false)));
uibutton(b1, 'Text', 'Save as...', 'ButtonPushedFcn', guard(@(~,~) save_table(true)));
table_label = uilabel(b1, 'Text', 'No table yet', 'FontAngle', 'italic');
% (new, user's request) tick or untick "include" for every recording in one click, e.g., untick all and tick the few to run
b_inc = uigridlayout(g1, [1 4], 'ColumnWidth', {80, 90, 210, '1x'}, 'Padding', [0 0 0 0]);
tick_all_button = uibutton(b_inc, 'Text', 'Tick all', 'ButtonPushedFcn', guard(@(~,~) set_include(true)), ...
    'Tooltip', 'Tick "include" for every recording');
untick_all_button = uibutton(b_inc, 'Text', 'Untick all', 'ButtonPushedFcn', guard(@(~,~) set_include(false)), ...
    'Tooltip', 'Untick "include" for every recording, then tick the few recordings to run');
uibutton(b_inc, 'Text', 'Organoids in the picture...', 'ButtonPushedFcn', guard(@on_show_organoids), ...
    'Tooltip', 'Show the organoids found in the picture of the selected row, with their numbers (for the column picture_organoid)');
include_label = uilabel(b_inc, 'Text', '', 'FontWeight', 'bold', 'FontColor', [0 0.45 0]);
rec_table =uitable(g1, 'Data', empty_rows(0), 'ColumnEditable', true, 'SelectionType', 'row', 'Multiselect', 'on', ...
    'CellEditCallback', guard(@on_table_edit), 'ColumnWidth', {55, 330, 170, 90, 90, 80, 100, 150, 70, 70, 70, 110, 90, 120});
uilabel(g1, 'Text', 'Tip: you can also edit the table in Excel and open it again here. Save the table before closing; Run saves it for you.', ...
    'FontColor', [0.4 0.4 0.4]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tab 2: settings

g2 = uigridlayout(tab_set, [1 2], 'ColumnWidth', {390, '1x'});
left = uigridlayout(g2, [3 1], 'RowHeight', {'fit', 'fit', 'fit'}, 'Padding', [0 0 0 0], 'Scrollable', 'on');   % (the rows of each panel
                                                                            %  that the choices do not need are hidden: see "update_settings_view")

% what to run
p_what = uipanel(left, 'Title', 'What to run');
gw = uigridlayout(p_what, [5 1], 'RowHeight', {22, 18, 22, 0, 'fit'});
what_mode = uidropdown(gw, 'Items', {'The whole pipeline (from the spike-sorted data to the Excel tables)', 'From a stage to the end', 'Only some stages'}, ...
    'ItemsData', {'all', 'from', 'only'}, 'ValueChangedFcn', guard(@(~,~) update_settings_view()));
from_label = uilabel(gw, 'Text', 'First stage:');
from_stage = uidropdown(gw, 'Items', cellstr(STAGE_TEXT(STAGES ~= "ShapeSpikeTrainData_3Brain")), 'ItemsData', cellstr(END_TO_END), ...
    'Value', 'Analyze_spikeTrain', 'ValueChangedFcn', guard(@(~,~) update_summary()));
% (revised, user's report: with a list, the selection was not visible enough) one tick box per stage
only_grid = uigridlayout(gw, [numel(STAGES) 1], 'RowHeight', repmat({18}, 1, numel(STAGES)), 'RowSpacing', 3, 'Padding', [4 0 0 0]);
only_checks = gobjects(1, numel(STAGES));
for istage = 1:numel(STAGES)
    only_checks(istage) = uicheckbox(only_grid, 'Text', STAGE_TEXT(istage), 'Value', ismember(STAGES(istage), ["Analyze_spikeTrain", "Compile_MEAproject"]), ...
        'ValueChangedFcn', guard(@(~,~) update_summary()));
end
what_summary = uilabel(gw, 'WordWrap', 'on', 'FontWeight', 'bold', 'Text', '');   % the stages of this run, for every choice above

% parameters
p_par = uipanel(left, 'Title', 'Parameters');
gp = uigridlayout(p_par, [5 2], 'RowHeight', {22, 22, 'fit', 22, 90}, 'ColumnWidth', {'1x', 70});
param_mode = uidropdown(gp, 'Items', {'Default values (pipeline_parameters.m)', 'The values saved for each recording by its last run', 'Another parameter file'}, ...
    'ItemsData', {'default', 'saved', 'file'}, 'ValueChangedFcn', guard(@on_param_mode));
param_mode.Layout.Column = [1 2];
param_file_field = uieditfield(gp, 'text', 'Editable', 'off', 'Placeholder', 'parameter file (.m)');
param_file_button = uibutton(gp, 'Text', 'Browse...', 'ButtonPushedFcn', guard(@on_browse_param_file));
param_help = uilabel(gp, 'WordWrap', 'on', 'FontColor', [0.4 0.4 0.4], 'Text', '');
param_help.Layout.Column = [1 2];
renew_label = uilabel(gp, 'Text', 'Take these functions from the default file (after an update):', 'WordWrap', 'on');
renew_label.Layout.Column = [1 2];
renew_list = uilistbox(gp, 'Items', {}, 'Multiselect', 'on', 'Value', {}, 'ValueChangedFcn', guard(@(~,~) update_summary()), ...
    'Tooltip', ['Only needed when the pipeline was updated after the last run and a function changed its parameters: those ' ...
    'functions take their values from the default file; every other value (e.g., the settings of a recording) is kept']);
renew_list.Layout.Column = [1 2];

% options
p_opt = uipanel(left, 'Title', 'Options');
go = uigridlayout(p_opt, [6 2], 'RowHeight', {32, 22, 18, 22, 22, 22}, 'ColumnWidth', {'1x', 70});
npy_label = uilabel(go, 'Text', 'npy-matlab folder (needed to read Kilosort data, stage 1; it contains readNPY.m):', 'WordWrap', 'on');
npy_label.Layout.Column = [1 2];
npy_field = uieditfield(go, 'text', 'Value', char(get_pref('npy_folder', "")), 'Placeholder', 'folder with readNPY.m', ...
    'ValueChangedFcn', guard(@(src,~) remember('npy_folder', string(src.Value))));
uibutton(go, 'Text', 'Browse...', 'ButtonPushedFcn', guard(@on_browse_npy));
results_label = uilabel(go, 'Text', 'Results folder (Excel files; empty: the folder of the table):');
results_label.Layout.Column = [1 2];
project_field = uieditfield(go, 'text', 'Placeholder', 'the folder of the recordings table', 'ValueChangedFcn', guard(@(~,~) update_summary()));
uibutton(go, 'Text', 'Browse...', 'ButtonPushedFcn', guard(@on_browse_project));
stop_error_box = uicheckbox(go, 'Text', 'Stop at the first error', 'Value', false, 'ValueChangedFcn', guard(@(~,~) update_summary()), ...
    'Tooltip', 'Otherwise the run continues with the next recording, and the errors are listed at the end');
stop_error_box.Layout.Column = [1 2];
verbose_box = uicheckbox(go, 'Text', 'Show the messages of each stage in the MATLAB Command Window', 'Value', false, ...
    'ValueChangedFcn', guard(@(~,~) update_summary()), 'Tooltip', 'They are always saved in the log file of each recording');
verbose_box.Layout.Column = [1 2];

% changes of single values
p_chg = uipanel(g2, 'Title', 'Change single values for this run (on top of the parameters chosen on the left)');
gc = uigridlayout(p_chg, [5 1], 'RowHeight', {22, 22, '1x', 'fit', 'fit'});
gs = uigridlayout(gc, [1 3], 'ColumnWidth', {190, '1x', 200}, 'Padding', [0 0 0 0]);
view_mode = uidropdown(gs, 'Items', {'Common adjustments', 'All parameters (advanced)'}, 'ItemsData', {'common', 'all'}, ...
    'ValueChangedFcn', guard(@(~,~) fill_param_table()), 'Tooltip', sprintf(['Common adjustments: your list of up to %d values that may need ' ...
    'changing for a recording (default: common_parameters.m). All parameters: every value of the parameter file'], MAX_COMMON));
search_field = uieditfield(gs, 'text', 'Placeholder', 'Search (name, function or words of its description)', ...
    'ValueChangingFcn', guard(@(~,ev) fill_param_table(ev.Value)));
changed_only_box = uicheckbox(gs, 'Text', 'Show only the changed values', 'ValueChangedFcn', guard(@(~,~) fill_param_table()));
gs2 = uigridlayout(gc, [1 5], 'ColumnWidth', {125, 105, 240, 100, '1x'}, 'Padding', [0 0 0 0]);
uibutton(gs2, 'Text', 'Undo all changes', 'ButtonPushedFcn', guard(@on_clear_changes));
uibutton(gs2, 'Text', 'Save as file...', 'ButtonPushedFcn', guard(@on_save_param_file), ...
    'Tooltip', 'Save the parameters with your changes as a new parameter file (to reuse them, choose "Another parameter file")');
list_button = uibutton(gs2, 'Text', 'Remove from this list', 'ButtonPushedFcn', guard(@on_list_button), ...
    'Tooltip', sprintf('Select rows first. The list of common adjustments holds up to %d parameters; it is kept for you (MATLAB preferences)', MAX_COMMON));
uibutton(gs2, 'Text', 'Default list', 'ButtonPushedFcn', guard(@on_default_list), ...
    'Tooltip', 'Go back to the default list of common adjustments (common_parameters.m)');
list_label = uilabel(gs2, 'Text', '', 'FontColor', [0.4 0.4 0.4]);
param_table = uitable(gc, 'ColumnEditable', [false false true false false], 'CellEditCallback', guard(@on_param_edit), ...   % (column 3: new value)
    'CellSelectionCallback', guard(@on_param_select));
details_label = uilabel(gc, 'WordWrap', 'on', 'FontColor', [0.25 0.25 0.25], 'Text', 'Click a row to read what it does.');
changes_label = uilabel(gc, 'WordWrap', 'on', 'Text', 'No value changed. Type a new value in the column "new value"; delete it to go back to the value of the file.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tab 3: run

g3 = uigridlayout(tab_run, [4 1], 'RowHeight', {110, '1x', 130, 26});
p_sum = uipanel(g3, 'Title', 'This run');
gsum = uigridlayout(p_sum, [1 2], 'ColumnWidth', {'1x', 120});
summary_area = uitextarea(gsum, 'Editable', 'off', 'Value', {''});
uibutton(gsum, 'Text', 'Copy command', 'ButtonPushedFcn', guard(@on_copy_command), ...
    'Tooltip', 'Copy the equivalent run_pipeline command (to paste in the MATLAB Command Window or in a script)');
progress_table = uitable(g3, 'Data', cell(0,3), 'ColumnName', {'Recording', 'Group', ''}, 'SelectionType', 'row', 'Multiselect', 'off');
log_area = uitextarea(g3, 'Editable', 'off', 'Value', {'The messages of the run appear here.'});
gr = uigridlayout(g3, [1 4], 'ColumnWidth', {200, 160, 160, '1x'}, 'Padding', [0 0 0 0]);
uibutton(gr, 'Text', 'Open the log of the selected row', 'ButtonPushedFcn', guard(@on_open_log), ...
    'Tooltip', 'Messages of every stage of the selected recording (or group), saved at each run');
uibutton(gr, 'Text', 'Open the results folder', 'ButtonPushedFcn', guard(@on_open_results));
uibutton(gr, 'Text', 'Open the Excel results', 'ButtonPushedFcn', guard(@on_open_excel));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tab 4: comparison of groups by their neurons (new, November 2026; built by "compare_tab")

compare = compare_tab(tab_cmp, struct('table_file', @tab1_table_file, 'results_folder', @results_folder, 'message', @message, ...
    'guard', @guard, 'with_dialog', @with_dialog, 'open_file', @open_file, 'test', opt.test));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tab 6: plots of the compiled measures, comparing groups (new, November 2026; built by "features_tab")

features = features_tab(tab_feat, struct('results_folder', @results_folder, 'message', @message, 'guard', @guard, ...
    'with_dialog', @with_dialog, 'open_file', @open_file, 'test', opt.test));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tab 5: visualization (new, November 2026; see "plotMEArecording")
% (revised, user's choices, 14 November 2026) a preview of the plot instead of the equivalent command; the options under the units

g4 = uigridlayout(tab_viz, [1 2], 'ColumnWidth', {440, '1x'});
left4 = uigridlayout(g4, [5 1], 'RowHeight', {'fit', 'fit', 'fit', 'fit', 'fit'}, 'Padding', [0 0 0 0], 'Scrollable', 'on');
right4 = uigridlayout(g4, [3 1], 'RowHeight', {'1x', 'fit', 40}, 'Padding', [0 0 0 0]);   % (preview, result, buttons)

% recording and dataset
p_vrec = uipanel(left4, 'Title', 'Recording');
gv1 = uigridlayout(p_vrec, [4 2], 'RowHeight', {22, 18, 22, 'fit'}, 'ColumnWidth', {'1x', 95});
viz_rec = uidropdown(gv1, 'Items', {'(open or build a recordings table in tab 1)'}, 'ItemsData', {0}, ...
    'ValueChangedFcn', guard(@(~,~) viz_choose_recording()), 'Tooltip', 'A recording of the table of tab 1');
uibutton(gv1, 'Text', 'Refresh', 'ButtonPushedFcn', guard(@(~,~) viz_fill_recordings()), 'Tooltip', 'Read the recordings table of tab 1 again');
viz_dataset_label = uilabel(gv1, 'Text', 'Dataset (file of the pipeline, after the quality control):');
viz_dataset_label.Layout.Column = [1 2];
viz_dataset = uidropdown(gv1, 'Items', {'(choose a recording)'}, 'ItemsData', {''}, 'ValueChangedFcn', guard(@(src,~) viz_load(src.Value)), ...
    'Tooltip', 'The analyzed dataset (_curated or _integrated) is chosen first; the cropped datasets have the label of their organoid');
uibutton(gv1, 'Text', 'Other file...', 'ButtonPushedFcn', guard(@on_viz_file), 'Tooltip', 'Choose the "_QCmanual.mat" file of any recording');
viz_info = uilabel(gv1, 'WordWrap', 'on', 'Text', 'No dataset read yet.', 'FontColor', [0.3 0.3 0.3]);
viz_info.Layout.Column = [1 2];

% what to draw
p_vwhat = uipanel(left4, 'Title', 'What to draw');
gv2 = uigridlayout(p_vwhat, [2 1], 'RowHeight', {22, 'fit'});
viz_type = uidropdown(gv2, 'Items', {'Raster plot', 'Bursts of each unit (raster plot)', 'Electrode map (firing rate)', ...
    'Spike waveforms (Kilosort templates)', 'Video of the electrode map (MP4)'}, 'ItemsData', {'raster', 'bursts', 'map', 'waveforms', 'video'}, ...
    'ValueChangedFcn', guard(@(~,~) viz_update_view()));
viz_help = uilabel(gv2, 'WordWrap', 'on', 'FontColor', [0.4 0.4 0.4], 'Text', '');

% time window
p_vtime = uipanel(left4, 'Title', 'Time window');
gv3 = uigridlayout(p_vtime, [2 4], 'RowHeight', {22, 22}, 'ColumnWidth', {65, '1x', 50, '1x'});
viz_time = uidropdown(gv3, 'Items', {'Whole recording'}, 'ItemsData', {0}, 'ValueChangedFcn', guard(@(~,~) viz_time_choice()), ...
    'Tooltip', 'The whole recording, one recording condition, or a window typed below (Custom)');
viz_time.Layout.Column = [1 4];
uilabel(gv3, 'Text', 'From (s):');
viz_from = uieditfield(gv3, 'numeric', 'ValueChangedFcn', guard(@(~,~) viz_custom_time()));
uilabel(gv3, 'Text', 'to (s):', 'HorizontalAlignment', 'right');
viz_to = uieditfield(gv3, 'numeric', 'ValueChangedFcn', guard(@(~,~) viz_custom_time()));

% units
p_vunits = uipanel(left4, 'Title', 'Units');
gv4 = uigridlayout(p_vunits, [1 2], 'RowHeight', {22}, 'ColumnWidth', {130, '1x'});
viz_units_mode = uidropdown(gv4, 'Items', {'All units', 'Only these units:'}, 'ItemsData', {'all', 'list'}, ...
    'ValueChangedFcn', guard(@(~,~) viz_update_view()));
viz_units = uieditfield(gv4, 'text', 'Placeholder', 'unit numbers, e.g., 1:20 35 40', 'ValueChangedFcn', guard(@(~,~) viz_units_typed()), ...
    'Tooltip', 'The numbers of the units, as in the raster plot (rows): single numbers and ranges such as 1:20, separated by spaces or commas');

% options (only the ones of the plot chosen are shown, see "viz_update_view")
p_vopt = uipanel(left4, 'Title', 'Options');
gv5 = uigridlayout(p_vopt, [14 2], 'ColumnWidth', {'1x', 150}, 'RowSpacing', 6);
viz_pop = uicheckbox(gv5, 'Text', 'Population firing rate', 'Value', true, 'ValueChangedFcn', guard(@(~,~) viz_update_view()), ...
    'Tooltip', 'Spikes of all the neurons drawn in each bin / number of neurons / bin (Hz per neuron), in bins of the length given below');
viz_pop.Layout.Column = [1 2];
viz_bin_label = uilabel(gv5, 'Text', 'Bin of the population firing rate (ms):');
viz_bin = uieditfield(gv5, 'numeric', 'Value', 200, 'Limits', [1 600000], 'ValueChangedFcn', guard(@(~,~) viz_update_view()), ...
    'Tooltip', 'Length of the bins of the population firing rate, in ms (default: 200 ms)');
viz_raster = uicheckbox(gv5, 'Text', 'Raster plot above the map (with a cursor that moves with the video)', 'Value', true, ...
    'ValueChangedFcn', guard(@(~,~) viz_update_view()));
viz_raster.Layout.Column = [1 2];
viz_cond = uicheckbox(gv5, 'Text', 'Start and name of each recording condition', 'Value', true, 'ValueChangedFcn', guard(@(~,~) viz_update_view()));
viz_cond.Layout.Column = [1 2];
viz_names_label = uilabel(gv5, 'Text', 'Names of the conditions:');
viz_names = uidropdown(gv5, 'Items', {'Condition and dose', 'File name'}, 'ItemsData', {'treatment', 'file'}, ...
    'ValueChangedFcn', guard(@(~,~) viz_names_changed()), 'Tooltip', ['Condition and dose: as the compiled tables group the recordings ' ...
    '(e.g., +Gabazine 20uM). File name: the part of the file name of the condition (e.g., gabazine_20uM)']);
viz_tunit_label = uilabel(gv5, 'Text', 'Time axis in:');
viz_tunit = uidropdown(gv5, 'Items', {'seconds (s)', 'milliseconds (ms)', 'minutes (min)'}, 'ItemsData', {'s', 'ms', 'min'}, ...
    'ValueChangedFcn', guard(@(~,~) viz_update_view()), 'Tooltip', 'Unit of the time axis (the window above is always in seconds)');
viz_measure_label = uilabel(gv5, 'Text', 'Rate of each unit:');
viz_measure = uidropdown(gv5, 'Items', {'Mean rate', '1 / median interval', '1 / 10% interval (bursts)'}, 'ItemsData', {'mean', 'median', 'fast'}, ...
    'ValueChangedFcn', guard(@(~,~) viz_update_view()), 'Tooltip', ['Mean rate: spikes / length of the window. 1 / median interval: rate ' ...
    'of the usual firing. 1 / 10% interval: rate inside the bursts (10th percentile of the intervals between spikes)']);
viz_limit_label = uilabel(gv5, 'Text', 'Top of the color scale (Hz; empty = automatic):');
viz_limit = uieditfield(gv5, 'text', 'Placeholder', 'automatic', 'ValueChangedFcn', guard(@(~,~) viz_update_view()), ...
    'Tooltip', 'Automatic: the 99th percentile of the rates of the electrodes with units');
viz_cmap_label = uilabel(gv5, 'Text', 'Colors of the map:');            % (user's request, November 2026)
viz_cmap = uidropdown(gv5, 'Items', {'Jet (blue to red)', 'Turbo (blue to red, smoother)', 'Parula (blue to yellow)', 'Hot (black to white)', ...
    'Gray (black to white)'}, 'ItemsData', {'jet', 'turbo', 'parula', 'hot', 'gray'}, 'ValueChangedFcn', guard(@(~,~) viz_update_view()), ...
    'Tooltip', 'Color palette of the firing rate on the electrode map and in the video');
viz_empty_label = uilabel(gv5, 'Text', 'Electrodes without neurons:');   % (user's request, November 2026)
viz_empty = uidropdown(gv5, 'Items', {'Black', 'White'}, 'ItemsData', {'black', 'white'}, 'ValueChangedFcn', guard(@(~,~) viz_update_view()), ...
    'Tooltip', 'Color of the electrodes without neurons on the map and in the video (the outline of the organoid is drawn in the other color)');
viz_outline = uicheckbox(gv5, 'Text', 'Outline of the organoid (from the cropping, or the photo of the chip)', 'Value', true, ...
    'ValueChangedFcn', guard(@(~,~) viz_update_view()));
viz_outline.Layout.Column = [1 2];
viz_step_label = uilabel(gv5, 'Text', 'Step of the video (s; one frame per step):');
viz_step = uieditfield(gv5, 'numeric', 'Value', 0.4, 'Limits', [0.01 60], 'ValueChangedFcn', guard(@(~,~) viz_update_view()));
viz_speed_label = uilabel(gv5, 'Text', 'Speed (times the real time):');
viz_speed = uieditfield(gv5, 'numeric', 'Value', 4, 'Limits', [0.01 1000], 'ValueChangedFcn', guard(@(~,~) viz_update_view()));
viz_save = uicheckbox(gv5, 'Text', 'Also save the figure (PNG) next to the dataset', 'Value', false, 'ValueChangedFcn', guard(@(~,~) viz_update_view()));
viz_save.Layout.Column = [1 2];

% preview, result and (user's request) the buttons Draw / Open below them
p_vprev = uipanel(right4, 'Title', 'Preview (lower resolution; Draw opens the full plot)');
p_vprev.Layout.Row = 1;
gv8 = uigridlayout(p_vprev, [2 1], 'RowHeight', {'fit', '1x'}, 'Padding', [4 4 4 4]);
viz_preview_label = uilabel(gv8, 'WordWrap', 'on', 'FontColor', [0.4 0.4 0.4], 'Text', 'The preview appears here once a dataset is read.');
viz_preview_box = uipanel(gv8, 'BorderType', 'none', 'BackgroundColor', [1 1 1]);
viz_result = uilabel(right4, 'WordWrap', 'on', 'Text', 'Choose a recording, what to draw, the time window and the units, then press Draw.');
viz_result.Layout.Row = 2;
gv6 = uigridlayout(right4, [1 3], 'ColumnWidth', {'1x', 140, 115}, 'Padding', [0 0 0 0]);
gv6.Layout.Row = 3;
viz_draw_button = uibutton(gv6, 'Text', 'Draw', 'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', [0.30 0.65 0.35], ...
    'FontColor', [1 1 1], 'ButtonPushedFcn', guard(@on_viz_draw), 'Tooltip', 'Draw the plot in its own window, full size (or make the video)');
viz_open_button = uibutton(gv6, 'Text', 'Open the saved file', 'Enable', 'off', 'ButtonPushedFcn', guard(@(~,~) viz_open(false)), ...
    'Tooltip', ['Opens the last file saved: a video, or a figure drawn with "Also save the figure (PNG)" ticked (available once one ' ...
    'is saved)']);
viz_folder_button = uibutton(gv6, 'Text', 'Open its folder', 'Enable', 'off', 'ButtonPushedFcn', guard(@(~,~) viz_open(true)), ...
    'Tooltip', 'Opens the folder of the last file saved (the folder of the dataset; available once a file is saved)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start

load_param_list(default_param_file);
update_settings_view();
if strlength(string(opt.table)) > 0
    open_table(string(opt.table));
else
    last = get_pref('last_table', "");
    if strlength(last) > 0 && isfile(last) && ~opt.test
        open_table(last);
    end
end
update_summary();
tabs.SelectionChangedFcn = guard(@on_tab_change);
viz_update_view();
fig.UserData.test = struct('open_table', @open_table, 'add_folder', @add_folder, 'start_folder', @start_folder, 'set_change', @set_change, 'run', @on_run, ...
    'add_to_common', @add_to_common, 'set_change_common', @set_change_common, 'remove_from_common', @remove_from_common, 'default_list', @on_default_list, 'common', @get_common, ...
    'stop', @on_stop, 'check', @check_settings, 'command', @command_text, 'save_param_file', @save_param_file, 'compare', compare, 'features', features, ...
    'viz_fill', @viz_fill_recordings, 'viz_load', @viz_load, 'viz_draw', @on_viz_draw, 'viz_inputs', @viz_inputs, 'viz_data', @get_viz_data, ...
    'controls', struct('group', group_field, 'what_mode', what_mode, 'from_stage', from_stage, 'only_checks', only_checks, 'param_mode', param_mode, ...
    'param_file', param_file_field, 'renew', renew_list, 'npy', npy_field, 'project', project_field, 'rec_table', rec_table, ...
    'progress', progress_table, 'log', log_area, 'status', status_label, 'param_table', param_table, 'stop_error', stop_error_box, ...
    'tabs', tabs, 'tab_viz', tab_viz, 'tab_feat', tab_feat, 'viz_rec', viz_rec, 'viz_dataset', viz_dataset, 'viz_type', viz_type, 'viz_time', viz_time, ...
    'viz_from', viz_from, 'viz_to', viz_to, 'viz_units_mode', viz_units_mode, 'viz_units', viz_units, 'viz_pop', viz_pop, ...
    'viz_raster', viz_raster, 'viz_cond', viz_cond, 'viz_measure', viz_measure, 'viz_limit', viz_limit, 'viz_outline', viz_outline, ...
    'viz_step', viz_step, 'viz_speed', viz_speed, 'viz_save', viz_save, 'viz_info', viz_info, 'viz_result', viz_result, ...
    'viz_draw', viz_draw_button, 'viz_options', gv5, 'viz_bin', viz_bin, 'viz_tunit', viz_tunit, 'viz_names', viz_names, ...
    'viz_preview', viz_preview_box, 'viz_preview_label', viz_preview_label, 'viz_cmap', viz_cmap, 'viz_empty', viz_empty));
if nargout > 0                                                              % (the window is returned only when asked: "fig" itself must stay,
    out = fig;                                                              %  the functions of the window use it)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   RECORDINGS TABLE   %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function T = empty_rows(n)
        % rows of the recordings table, with every column of "run_pipeline" (include ticked, source Kilosort)
        T = table(true(n,1), 'VariableNames', {'include'});
        for c = COLUMNS(2:end)
            if c == "source"
                T.source = categorical(repmat("Kilosort", n, 1), ["Kilosort", "3Brain"]);
            else
                T.(c) = strings(n, 1);
            end
        end
    end

    function open_table(file)
        % reads a recordings table (.xlsx or .csv); the missing optional columns are added
        file = string(file);
        try
            R = readtable(file, 'TextType', 'string', 'Delimiter', ',', 'FileType', file_type(file));
        catch
            R = readtable(file, 'TextType', 'string');
        end
        names = string(R.Properties.VariableNames);
        missing = setdiff(["folder", "file_name", "group"], names);
        if ~isempty(missing)
            message(sprintf('The table %s has no column %s (needed: folder, file_name, group).', file, strjoin(missing, ', ')), 'Recordings table');
            return
        end
        T = empty_rows(height(R));
        for c = COLUMNS
            if ~ismember(c, names)
                continue
            end
            v = R.(c);
            switch c
                case "include"
                    if isnumeric(v) || islogical(v)
                        T.include = ~(v == 0);              % (as "run_pipeline": empty = included)
                    else
                        T.include = ~ismember(strtrim(lower(string(v))), ["0", "false", "no"]);
                    end
                case "source"
                    s = string(v);   s(ismissing(s) | strlength(strtrim(s)) == 0) = "Kilosort";
                    s(strcmpi(s, "3Brain")) = "3Brain";     s(s ~= "3Brain") = "Kilosort";
                    T.source = categorical(s, ["Kilosort", "3Brain"]);
                otherwise
                    s = string(v);   s(ismissing(s)) = "";
                    T.(c) = s;
            end
        end
        for c = setdiff(names, COLUMNS, 'stable')                           % other columns of the user (e.g., notes): kept, at the right,
            v = R.(c);                                                      % and saved with the table
            if ~(isnumeric(v) || islogical(v))
                v = string(v);   v(ismissing(v)) = "";
            end
            T.(c) = v;
        end
        rec_table.Data = T;
        rec_table.ColumnWidth = [{55, 330, 170, 90, 90, 80, 100, 150, 70, 70, 70, 110, 90, 120}, repmat({'auto'}, 1, width(T) - numel(COLUMNS))];
        table_file = file;
        table_dirty = false;
        remember('last_table', file);
        update_table_label();
        update_summary();
        set_status(sprintf('Table opened: %d recordings.', height(T)));
    end

    function on_open_table(~, ~)
        [f, p] = with_dialog(@() uigetfile({'*.xlsx;*.csv', 'Recordings table (*.xlsx, *.csv)'}, 'Open a recordings table', char(start_folder())));
        if isequal(f, 0)
            return
        end
        if ~confirm_discard()
            return
        end
        open_table(string(fullfile(p, f)));
    end

    function on_add_folder(~, ~)
        d = with_dialog(@() uigetdir(char(start_folder()), 'Folder with the recordings (it and its subfolders are searched)'));
        if isequal(d, 0)
            return
        end
        add_folder(string(d));
    end

    function n = add_folder(d)
        % adds one row for each recording found in the folder d and its subfolders: a folder with files of the pipeline
        % ("<name>_raw.mat", "<name>_QCmanual.mat", ...) gives "<name>"; a folder with Kilosort output (spike_times.npy) and no file of the
        % pipeline gives the name of the folder (to be checked)
        found = find_recordings(d);
        T = rec_table.Data;
        known = T.folder + "|" + T.file_name;
        found = found(~ismember(found.folder + "|" + found.file_name, known), :);
        n = height(found);
        if n == 0
            message(sprintf('No new recording was found in %s (folders with Kilosort output, or with files of the pipeline such as <name>_raw.mat).', d), 'Add recordings');
            return
        end
        add_rows(found);
        added_folder = d;
        set_status(sprintf('%d recordings added from %s. Fill in their group, then save the table.', n, d));
    end

    function add_rows(R)
        D = rec_table.Data;
        for c = string(D.Properties.VariableNames(~ismember(D.Properties.VariableNames, R.Properties.VariableNames)))
            v = D.(c);                                                      % (the other columns of the user: empty for the new rows)
            if isnumeric(v)
                R.(c) = nan(height(R), 1);
            elseif islogical(v)
                R.(c) = false(height(R), 1);
            else
                R.(c) = strings(height(R), 1);
            end
        end
        rec_table.Data = [D; R(:, D.Properties.VariableNames)];
        table_dirty = true;
        update_table_label();
        update_summary();
    end

    function on_remove_rows(~, ~)
        sel = rec_table.Selection;
        if isempty(sel)
            message('Click the rows to remove first (hold Cmd or Ctrl to select several).', 'Remove rows');
            return
        end
        rows = selected_table_rows();
        T = rec_table.Data;
        T(rows, :) = [];
        rec_table.Data = T;
        table_dirty = true;
        update_table_label();
        update_summary();
    end

    function rows = selected_table_rows()
        % rows selected in the recordings table (with SelectionType 'row', "Selection" is the list of their numbers; before, only the
        % first one was read)
        rows = unique(rec_table.Selection(:));
        rows = rows(rows >= 1 & rows <= height(rec_table.Data));
    end

    function on_show_organoids(~, ~)
        % (November 2026) figure of the organoids found in the picture of the selected row ("OrganoidsFromPhoto"), with their numbers
        rows = selected_table_rows();
        if isempty(rows)
            message('Click the row of a recording with a picture first.', 'Organoids in the picture');
            return
        end
        T = rec_table.Data;     r = rows(1);
        pic = strtrim(T.picture(r));
        folder = strtrim(T.folder(r));
        if isempty(regexp(char(folder), '^(/|~|[A-Za-z]:[\\/]|\\\\)', 'once')) && strlength(table_file) > 0
            folder = string(fullfile(fileparts(table_file), folder));       % (relative to the folder of the table)
        end
        if strlength(pic) == 0                                              % (as the run: the photo in the folder of the recording)
            [pic, candidates] = recording_picture(folder, strtrim(T.file_name(r)));
            if strlength(pic) == 0 && numel(candidates) > 1
                message(sprintf('Several photos in the folder of this recording (%s): write the one to use in the column picture.', ...
                    strjoin(candidates, ', ')), 'Organoids in the picture');
                return
            end
        end
        if strlength(pic) == 0 || lower(pic) == "none"
            message('This recording has no photo of the chip (column picture, or the folder of the recording).', 'Organoids in the picture');
            return
        end
        if isempty(regexp(char(pic), '^(/|~|[A-Za-z]:[\\/]|\\\\)', 'once'))
            pic = string(fullfile(folder, pic));                            % (a file in the folder of the recording)
        end
        n = str2double(T.picture_n_organoids(r));
        if isnan(n), n = 0; end
        set_status(sprintf('Finding the organoids in %s...', pic));     drawnow;
        try
            res = OrganoidsFromPhoto(pic, struct('TEXT_SIDE', lower(strtrim(T.picture_3D_side(r))), 'N_ORGANOIDS', n));
            set_status(sprintf(['%d organoid(s) in the picture of row %d (the labels of organoid_label, or org1, org2..., are given in ' ...
                'this order at the first run); %s'], res.n_organoids, r, res.orientation.message));
        catch ME
            message(sprintf('The organoids could not be found in %s: %s', pic, ME.message), 'Organoids in the picture');
        end
    end

    function on_group_rows(~, ~)
        % the group typed in the field goes to the selected rows (a field of the window: the separate dialog used before did not
        % always pass the group to the table)
        group = strtrim(string(group_field.Value));
        rows = selected_table_rows();
        if strlength(group) == 0
            message('Type the group in the field on the left of the button first.', 'Group');
            return
        end
        if isempty(rows)
            message('Click the rows first (hold Cmd or Ctrl to select several).', 'Group');
            return
        end
        T = rec_table.Data;
        T.group(rows) = group;
        rec_table.Data = T;
        set_status(sprintf('Group "%s" set for %d recording(s). Save the table to keep it.', group, numel(rows)));
        table_dirty = true;
        update_table_label();
        update_summary();
    end

    function set_include(value)
        % ticks (true) or unticks (false) "include" for every recording of the table
        T = rec_table.Data;
        if height(T) == 0
            set_status('The table has no recording yet.');
            return
        end
        T.include(:) = value;
        rec_table.Data = T;
        table_dirty = true;
        update_table_label();
        update_summary();
        if value
            set_status(sprintf('All %d recordings are ticked.', height(T)));
        else
            set_status('All recordings are unticked: tick the ones to run (column include).');
        end
    end

    function update_include_label()
        % number of recordings ticked, next to the buttons Tick all / Untick all
        T = rec_table.Data;
        n = nnz(T.include);
        include_label.Text = sprintf('%d of %d recordings ticked (column include)', n, height(T));
        include_label.FontColor = ifelse(n == 0 && height(T) > 0, [0.8 0 0], [0 0.45 0]);
    end

    function on_table_edit(~, ~)
        table_dirty = true;
        update_table_label();
        update_summary();
    end

    function ok = save_table(ask)
        % saves the table (.xlsx or .csv); asks for a file name if it has none yet, or with "Save as..."
        ok = false;
        file = table_file;
        if ask || strlength(file) == 0
            [f, p] = with_dialog(@() uiputfile({'*.xlsx', 'Excel table (*.xlsx)'; '*.csv', 'Text table (*.csv)'}, 'Save the recordings table', char(fullfile(start_folder(), 'recordings.xlsx'))));
            if isequal(f, 0)
                return
            end
            file = string(fullfile(p, f));
        end
        T = rec_table.Data;
        W = T;
        W.include = double(T.include);
        W.source = string(T.source);
        try
            if isfile(file) && endsWith(lower(file), [".xlsx", ".xls"])
                delete(file);                                               % (writetable keeps the cells of an older, longer table)
            end
            writetable(W, file);
        catch ME
            message(sprintf('The table could not be saved: %s', ME.message), 'Save');
            return
        end
        table_file = file;
        table_dirty = false;
        remember('last_table', file);
        update_table_label();
        update_summary();
        set_status(sprintf('Table saved: %s', file));
        ok = true;
    end

    function update_table_label()
        if strlength(table_file) == 0
            txt = 'Not saved yet';
        else
            txt = char(table_file);
        end
        if table_dirty
            txt = [txt '  (modified)'];
        end
        table_label.Text = txt;
        table_label.Tooltip = txt;
    end

    function ok = confirm_discard()
        ok = true;
        if table_dirty && height(rec_table.Data) > 0 && ~opt.test
            choice = uiconfirm(fig, 'The table was modified and not saved. Discard the changes?', 'Recordings table', ...
                'Options', {'Discard', 'Cancel'}, 'DefaultOption', 2, 'CancelOption', 2);
            ok = strcmp(choice, 'Discard');
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   SETTINGS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function update_settings_view()
        % shows only the controls that the choices need (the other rows of the panels have no height)
        mode = what_mode.Value;
        from_label.Visible = ~strcmp(mode, 'all');
        from_stage.Visible = strcmp(mode, 'from');
        only_grid.Visible = strcmp(mode, 'only');
        if strcmp(mode, 'only')
            from_label.Text = 'Tick the stages to run (they run in the order of the pipeline):';
        else
            from_label.Text = 'First stage:';
        end
        gw.RowHeight = {22, 18*~strcmp(mode, 'all'), 22*strcmp(mode, 'from'), 21*numel(STAGES)*strcmp(mode, 'only'), 'fit'};
        pm = param_mode.Value;
        param_file_field.Visible = strcmp(pm, 'file');
        param_file_button.Visible = strcmp(pm, 'file');
        renew_label.Visible = strcmp(pm, 'saved');
        renew_list.Visible = strcmp(pm, 'saved');
        gp.RowHeight = {22, 22*strcmp(pm, 'file'), 'fit', 32*strcmp(pm, 'saved'), 90*strcmp(pm, 'saved')};
        param_file_field.Enable = strcmp(pm, 'file');
        param_file_button.Enable = strcmp(pm, 'file');
        renew_label.Enable = strcmp(pm, 'saved');
        renew_list.Enable = strcmp(pm, 'saved');
        switch pm
            case 'default'
                param_help.Text = 'Every value comes from the file pipeline_parameters.m of the pipeline.';
            case 'saved'
                param_help.Text = ['Each recording reruns with its own values of its last run (e.g., the settings of a recording ' ...
                    'that you changed). For a recording never run before, choose the default values.'];
            case 'file'
                param_help.Text = 'Every value comes from the file you choose (a copy of pipeline_parameters.m, edited).';
        end
        update_summary();
    end

    function on_param_mode(~, ~)
        if strcmp(param_mode.Value, 'file') && strlength(string(param_file_field.Value)) > 0
            load_param_list(string(param_file_field.Value));
        elseif ~strcmp(param_mode.Value, 'file')
            load_param_list(default_param_file);
        end
        update_settings_view();
    end

    function on_browse_param_file(~, ~)
        [f, p] = with_dialog(@() uigetfile({'*.m', 'Parameter file (*.m)'}, 'Choose a parameter file', char(start_folder())));
        if isequal(f, 0)
            return
        end
        param_file_field.Value = fullfile(p, f);
        load_param_list(string(param_file_field.Value));
        update_summary();
    end

    function on_browse_npy(~, ~)
        d = with_dialog(@() uigetdir(char(start_folder()), 'Folder of npy-matlab (it contains readNPY.m)'));
        if isequal(d, 0)
            return
        end
        npy_field.Value = d;
        remember('npy_folder', string(d));
        update_summary();
    end

    function on_browse_project(~, ~)
        d = with_dialog(@() uigetdir(char(start_folder()), 'Folder for the results of each group (Excel files)'));
        if isequal(d, 0)
            return
        end
        project_field.Value = d;
        update_summary();
    end

    function load_param_list(file)
        % reads the parameters of a parameter file (values and descriptions) for the list of the tab Settings
        try
            params = read_parameter_list(file);
        catch ME
            message(sprintf('The parameter file %s could not be read: %s', file, ME.message), 'Parameters');
            return
        end
        functions = unique(params.Function, 'stable');
        renew_list.Items = cellstr(functions);
        renew_list.ItemsData = cellstr(functions);
        renew_list.Value = {};
        fill_param_table();
    end

    function fill_param_table(search)
        % shows the parameters that match the search (and only the changed ones, if asked)
        if isempty(params)
            return
        end
        if nargin < 1
            search = search_field.Value;
        end
        pkeys = params.Function + "." + params.Parameter;
        % (two levels) the common adjustments, in the order of "common_parameters", with their short name and helptext; or every parameter
        common_view = strcmp(view_mode.Value, 'common');
        if common_view
            [found, where] = ismember(common.key, pkeys);
            order = where(found);
            label = common.setting(found);      helptext = common.help(found);      topic = common.topic(found);
            added = strlength(label) == 0;                                  % (parameters added by the user: no short name, the description
            label(added) = params.Parameter(order(added)) + " (" + params.Function(order(added)) + ")";   %  of the file)
            helptext(added) = params.Description(order(added));
        else
            order = (1:height(params))';
            label = params.Parameter;           helptext = params.Description;      topic = params.Stage;
        end
        show = true(numel(order), 1);
        words = split(lower(strtrim(string(search))));
        words = words(strlength(words) > 0);
        text_all = lower(topic + " " + label + " " + pkeys(order) + " " + helptext);
        for w = words'
            show = show & contains(text_all, w);
        end
        if changed_only_box.Value
            show = show & isKey_list(pkeys(order));
        end
        idx = order(show);
        new_value = strings(numel(idx), 1);
        for k=1:numel(idx)
            if isKey(changes, char(pkeys(idx(k))))
                new_value(k) = changes(char(pkeys(idx(k))));
            end
        end
        file_value = params.Value(idx);
        if common_view                                                      % (parameters shown as a period: a rate in events per second,
            period = common.display(found);  period = period(show) == "period";   %  shown and typed in seconds, 1/rate)
            for k = find(period)'
                file_value(k) = string(num2str(1/params.Default{idx(k)}, 8));
                if strlength(new_value(k)) > 0
                    new_value(k) = string(num2str(1/parse_value(new_value(k), params.Default{idx(k)}), 8));
                end
            end
        end
        if common_view
            D = table(label(show), file_value, new_value, helptext(show), pkeys(idx), ...
                'VariableNames', {'setting', 'in the file', 'new value', 'what it does', 'parameter'});
            param_table.ColumnWidth = {250, 95, 80, 420, 330};
        else
            D = table(pkeys(idx), params.Value(idx), new_value, params.Description(idx), params.Stage(idx), ...
                'VariableNames', {'parameter', 'in the file', 'new value', 'description', 'stage'});
            param_table.ColumnWidth = {300, 95, 80, 600, 220};
        end
        param_table.Data = D;
        param_table.UserData = idx;                                         % rows of "params" shown
        removeStyle(param_table);
        changed = find(strlength(new_value) > 0);
        if ~isempty(changed)
            addStyle(param_table, uistyle('BackgroundColor', COLOR_CHANGED), 'row', changed);
        end
        in_list = ismember(pkeys(idx), common.key);
        if ~common_view && any(in_list)
            addStyle(param_table, uistyle('FontWeight', 'bold'), 'row', find(in_list));
        end
        selected_rows = [];
        if common_view
            list_button.Text = 'Remove from this list';
        else
            list_button.Text = 'Add to the common adjustments';
        end
        list_label.Text = sprintf('List%s: %d of %d', ifelse(common_view, '', ' (bold)'), height(common), MAX_COMMON);
        update_changes_label();
    end

    function tf = isKey_list(keys)
        tf = false(numel(keys), 1);
        for k=1:numel(keys)
            tf(k) = isKey(changes, char(keys(k)));
        end
    end

    function on_param_select(~, ev)
        % the whole description of the row clicked (the column is too narrow for most of them)
        if isempty(ev.Indices) || isempty(param_table.UserData)
            return
        end
        selected_rows = unique(ev.Indices(:,1));
        idx = param_table.UserData(ev.Indices(1));
        key = params.Function(idx) + "." + params.Parameter(idx);
        k = find(common.key == key, 1);
        if strcmp(view_mode.Value, 'common') && ~isempty(k) && strlength(common.help(k)) > 0
            details_label.Text = char(common.setting(k) + ": " + common.help(k) + "  (" + key + "; in the file: " + params.Description(idx) + ")");
        else
            details_label.Text = char(key + ": " + params.Description(idx));
        end
    end

    function C = load_common_list()
        % the list of common adjustments of this user (MATLAB preferences), or the default one; parameters added by the user get their
        % short name and description from the parameter file (see "fill_param_table")
        C = common_default(1:min(height(common_default), MAX_COMMON), :);
        if opt.test || ~ispref('CellCircuit', 'common_keys')
            return
        end
        list = string(getpref('CellCircuit', 'common_keys'));
        list = list(1:min(numel(list), MAX_COMMON));
        C = common_rows(list(:));
    end

    function C = common_rows(list)
        % rows of the list of common adjustments for the parameters "list" (Function.PARAMETER), in this order
        n = numel(list);
        C = table(list, strings(n,1), strings(n,1), strings(n,1), strings(n,1), 'VariableNames', {'key', 'topic', 'setting', 'help', 'display'});
        [found, where] = ismember(list, common_default.key);
        C(found, :) = common_default(where(found), :);
    end

    function save_common_list()
        if ~opt.test
            setpref('CellCircuit', 'common_keys', cellstr(common.key));
        end
    end

    function ok = add_to_common(list)
        % adds parameters (Function.PARAMETER) to the list of common adjustments, at most MAX_COMMON in all
        list = string(list(:));
        list = unique(list(~ismember(list, common.key)), 'stable');
        ok = height(common) + numel(list) <= MAX_COMMON;
        if ~ok
            message(sprintf(['The common adjustments hold at most %d parameters (now %d): remove one first (view "Common adjustments", ' ...
                'select it, "Remove from this list").'], MAX_COMMON, height(common)), 'Common adjustments');
            return
        end
        common = [common; common_rows(list)];
        save_common_list();
        fill_param_table();
        set_status(sprintf('%d parameter(s) added to the common adjustments (%d of %d).', numel(list), height(common), MAX_COMMON));
    end

    function remove_from_common(list)
        common(ismember(common.key, string(list)), :) = [];
        save_common_list();
        fill_param_table();
        set_status(sprintf('Removed from the common adjustments (%d of %d left).', height(common), MAX_COMMON));
    end

    function on_list_button(~, ~)
        if isempty(selected_rows) || isempty(param_table.UserData)
            message('Click the rows first (hold Cmd or Ctrl to select several).', 'Common adjustments');
            return
        end
        idx = param_table.UserData(selected_rows(selected_rows <= numel(param_table.UserData)));
        list = params.Function(idx) + "." + params.Parameter(idx);
        if strcmp(view_mode.Value, 'common')
            remove_from_common(list);
        else
            add_to_common(list);
        end
    end

    function on_default_list(~, ~)
        if ~opt.test
            choice = uiconfirm(fig, sprintf('Go back to the default list of common adjustments (%d parameters, common_parameters.m)?', ...
                min(height(common_default), MAX_COMMON)), 'Common adjustments', 'Options', {'Default list', 'Cancel'}, 'DefaultOption', 2, 'CancelOption', 2);
            if ~strcmp(choice, 'Default list')
                return
            end
            if ispref('CellCircuit', 'common_keys')
                rmpref('CellCircuit', 'common_keys');
            end
        end
        common = common_default(1:min(height(common_default), MAX_COMMON), :);
        fill_param_table();
        set_status('The default list of common adjustments is back.');
    end

    function C = get_common()
        C = common;                                                         % (for the tests: the list as it is now)
    end

    function varargout = with_dialog(fn)
        % (new, user's report) opens a file or folder dialog in front: on macOS the dialogs opened behind the window, so the window is
        % hidden while the dialog is open, and shown again (in front) when it closes
        fig.Visible = 'off';
        drawnow;
        restore = onCleanup(@() show_window());
        [varargout{1:nargout}] = fn();
        clear restore
    end

    function show_window()
        if isvalid(fig)
            fig.Visible = 'on';
            figure(fig);
        end
    end

    function cb = guard(fn)
        % (new) callback that reports its errors in the window (status bar and a message), with the details in the Command Window: before,
        % an error of a button or a list was only printed in the Command Window, and the window just seemed not to react
        cb = @(src, ev) run_guarded(fn, src, ev);
    end

    function run_guarded(fn, src, ev)
        try
            if nargin(fn) == 0
                fn();
            else
                fn(src, ev);
            end
        catch ME
            fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
            if isvalid(fig)
                set_status(['Error: ' ME.message ' (details in the MATLAB Command Window)']);
                status_label.FontColor = [0.8 0 0];
                message(sprintf('Something went wrong:\n%s\n\nThe details are in the MATLAB Command Window. If the window was opened before an update of CellCircuit, close it and open it again.', ME.message), 'Error');
            end
        end
    end

    function remember(name, value)
        % MATLAB preferences of the window (not in the tests)
        if ~opt.test
            set_pref(name, value);
        end
    end

    function on_param_edit(~, ev)
        row = ev.Indices(1);
        idx = param_table.UserData(row);
        key = params.Function(idx) + "." + params.Parameter(idx);
        if strcmp(view_mode.Value, 'common')
            set_change_common(key, string(ev.EditData));
        else
            set_change(key, string(ev.EditData));
        end
    end

    function ok = set_change_common(key, txt)
        % a value typed in the list of common adjustments: a parameter shown as a period (seconds) is changed to the rate, 1/(seconds), written
        % as an expression (e.g., "1/(300)"), so that the full list and the run get the value of the parameter file
        txt = strtrim(string(txt));
        k = find(common.key == string(key), 1);
        if ~isempty(k) && common.display(k) == "period" && strlength(txt) > 0
            seconds = str2double(txt);
            if ~(isfinite(seconds) && seconds > 0)
                message(sprintf('The value "%s" is not valid for %s: a time in seconds, above 0, is expected (e.g., 210).', txt, key), 'Change a value');
                ok = false;
                fill_param_table();
                return
            end
            txt = "1/(" + txt + ")";
        end
        ok = set_change(key, txt);
    end

    function ok = set_change(key, txt)
        % records a value typed by the user (empty: no change); the text is checked against the type of the value in the file
        key = char(key);
        txt = strtrim(string(txt));
        ok = true;
        k = find(params.Function + "." + params.Parameter == string(key), 1);
        if isempty(k)
            message(sprintf('Unknown parameter %s.', key), 'Change a value');
            ok = false;
            return
        end
        if strlength(txt) == 0
            if isKey(changes, key)
                remove(changes, key);
            end
        else
            [~, valid, why] = parse_value(txt, params.Default{k});
            if ~valid
                message(sprintf('The value "%s" is not valid for %s: %s', txt, key, why), 'Change a value');
                ok = false;
            else
                changes(key) = txt;
            end
        end
        fill_param_table();
        update_summary();
    end

    function on_clear_changes(~, ~)
        remove(changes, keys(changes));
        fill_param_table();
        update_summary();
    end

    function update_changes_label()
        n = changes.Count;
        if n == 0
            changes_label.Text = 'No value changed. Type a new value in the column "new value"; delete it to go back to the value of the file.';
        else
            changes_label.Text = sprintf('%d value(s) changed for this run (highlighted): %s', n, strjoin(string(keys(changes)), ', '));
        end
    end

    function list = changes_cell()
        % the changes as {'Function.PARAMETER', value, ...} (option 'changes' of "run_pipeline")
        list = {};
        names = keys(changes);
        for k=1:numel(names)
            row = find(params.Function + "." + params.Parameter == string(names{k}), 1);
            value = parse_value(changes(names{k}), params.Default{row});
            list = [list, {names{k}, value}]; %#ok<AGROW>
        end
    end

    function on_save_param_file(~, ~)
        [f, p] = with_dialog(@() uiputfile({'*.m', 'Parameter file (*.m)'}, 'Save the parameters with your changes', char(fullfile(start_folder(), 'my_parameters.m'))));
        if isequal(f, 0)
            return
        end
        save_param_file(string(fullfile(p, f)));
    end

    function save_param_file(file)
        % copy of the parameter file in use, with the changed values added at its end (they replace the values above them)
        base = current_param_file();
        if strcmp(param_mode.Value, 'saved')
            base = default_param_file;
        end
        try
            txt = fileread(base);
            names = keys(changes);
            extra = sprintf('\n\n%%%% Values changed in the window CellCircuit (%s), on top of %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), base);
            for k=1:numel(names)
                extra = [extra sprintf('P.%s = %s;\n', names{k}, changes(names{k}))]; %#ok<AGROW>
            end
            fid = fopen(file, 'w');
            fwrite(fid, [txt extra]);
            fclose(fid);
            load_parameters(file);                                          % check that the new file can be read
        catch ME
            message(sprintf('The parameter file could not be saved: %s', ME.message), 'Save the parameters');
            return
        end
        set_status(sprintf('Parameters saved: %s (to use them, choose "Another parameter file").', file));
    end

    function file = current_param_file()
        file = default_param_file;
        if strcmp(param_mode.Value, 'file') && strlength(string(param_file_field.Value)) > 0
            file = string(param_file_field.Value);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   RUN   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function args = run_options()
        % options of "run_pipeline" for the settings of the window (without 'progress')
        args = {};
        switch what_mode.Value
            case 'from'
                if ~strcmp(from_stage.Value, 'ShapeSpikeTrainData')
                    args = [args, {'from', string(from_stage.Value)}];
                end
            case 'only'
                args = [args, {'stages', selected_only()}];
        end
        switch param_mode.Value
            case 'saved'
                args = [args, {'parameters', "saved"}];
                if ~isempty(renew_list.Value)
                    args = [args, {'renew', string(renew_list.Value)}];
                end
            case 'file'
                args = [args, {'parameters', string(param_file_field.Value)}];
        end
        list = changes_cell();
        if ~isempty(list)
            args = [args, {'changes', list}];
        end
        if strlength(string(project_field.Value)) > 0
            args = [args, {'project', string(project_field.Value)}];
        end
        if stop_error_box.Value
            args = [args, {'stop_on_error', 1}];
        end
        if verbose_box.Value
            args = [args, {'verbose', 1}];
        end
        if strlength(string(npy_field.Value)) > 0 && needs_npy()
            args = [args, {'extra_path', string(npy_field.Value)}];
        end
    end

    function st = run_stages()
        switch what_mode.Value
            case 'all'
                st = END_TO_END;
            case 'from'
                st = END_TO_END(find(END_TO_END == string(from_stage.Value)):end);
            case 'only'
                st = STAGES(ismember(STAGES, selected_only()));
        end
    end

    function st = selected_only()
        % the stages ticked in "Only some stages"
        st = STAGES([only_checks.Value]);
    end

    function tf = needs_npy()
        tf = any(run_stages() == "ShapeSpikeTrainData");
    end

    function txt = command_text()
        % the equivalent command (to run without the window, or in a script)
        file = table_file;
        if strlength(file) == 0
            file = "recordings.xlsx";
        end
        args = run_options();
        parts = "run_pipeline(" + value_code(file);
        for k=1:2:numel(args)
            parts = parts + ", '" + args{k} + "', " + value_code(args{k+1});
        end
        txt = parts + ")";
    end

    function update_summary()
        update_include_label();
        if ~exist('summary_area', 'var') || ~isvalid(summary_area)
            return
        end
        T = rec_table.Data;
        inc = T(T.include, :);
        groups = unique(inc.group(strlength(inc.group) > 0));
        st = run_stages();
        lines = strings(0,1);
        lines(end+1) = sprintf('Recordings: %d of %d included, in %d group(s): %s', height(inc), height(T), numel(groups), strjoin(groups, ', '));
        lines(end+1) = sprintf('Stages: %s', strjoin(STAGE_SHORT(ismember(STAGES, st)), ' > '));
        if isempty(st)
            what_summary.Text = 'No stage selected: tick at least one.';
            what_summary.FontColor = [0.8 0 0];
        else
            what_summary.Text = char("This run: " + strjoin(STAGE_SHORT(ismember(STAGES, st)), " > "));
            what_summary.FontColor = [0 0.45 0];
        end
        switch param_mode.Value
            case 'default'
                p = "default values (pipeline_parameters.m)";
            case 'saved'
                p = "the values saved for each recording by its last run";
                if ~isempty(renew_list.Value)
                    p = p + ", with " + strjoin(string(renew_list.Value), ", ") + " from the default file";
                end
            case 'file'
                p = "the file " + string(param_file_field.Value);
        end
        lines(end+1) = "Parameters: " + p + sprintf('; %d value(s) changed', changes.Count);
        if any(st == "Compile_MEAproject")                                  % (new) where the Excel results go (by default, the folder of the
            lines(end+1) = "Excel results: " + fullfile(results_folder(), "<group>_Compiled.xlsx") + " (one file per group)";   %  table)
        end
        lines(end+1) = "Equivalent command:  " + command_text();
        summary_area.Value = cellstr(lines);
    end

    function on_copy_command(~, ~)
        clipboard('copy', char(command_text()));
        set_status('The command was copied: paste it in the MATLAB Command Window or in a script.');
    end

    function problems = check_settings()
        % list of the problems that would stop the run (empty: ready)
        problems = strings(0,1);
        T = rec_table.Data;
        inc = find(T.include);
        if isempty(inc)
            problems(end+1) = "No recording is included (tab 1, column include).";
        end
        base = start_folder();
        for r = inc'
            name = T.file_name(r);
            if strlength(name) == 0
                problems(end+1) = sprintf("Row %d: the column file_name is empty.", r); %#ok<AGROW>
                name = "row " + r;
            end
            folder = T.folder(r);
            if strlength(folder) == 0
                problems(end+1) = sprintf("%s: the column folder is empty.", name); %#ok<AGROW>
            elseif ~isfolder(folder) && ~isfolder(fullfile(base, folder))
                problems(end+1) = sprintf("%s: the folder %s does not exist.", name, folder); %#ok<AGROW>
            end
            if strlength(T.group(r)) == 0
                problems(end+1) = sprintf("%s: the column group is empty.", name); %#ok<AGROW>
            end
        end
        folders = T.folder(inc);
        [u, ~, j] = unique(folders);
        dup = u(accumarray(j, 1) > 1 & strlength(u) > 0);
        if ~isempty(dup) && any(ismember(["QualityControl_General", "ShapeSpikeTrainData"], run_stages()))
            problems(end+1) = "Several recordings share the folder " + strjoin(dup, ", ") + ": each recording must be in its own folder.";
        end
        if isempty(run_stages())
            problems(end+1) = "No stage is selected (tab 2, What to run).";
        end
        if strcmp(param_mode.Value, 'file') && ~isfile(string(param_file_field.Value))
            problems(end+1) = "Choose the parameter file (tab 2, Parameters).";
        end
        if needs_npy() && exist('readNPY', 'file') ~= 2
            npy = string(npy_field.Value);
            if strlength(npy) == 0 || ~isfile(fullfile(npy, 'readNPY.m'))
                problems(end+1) = "Stage 1 reads the Kilosort data with npy-matlab: choose its folder (tab 2, Options), the one that contains readNPY.m.";
            end
        end
    end

    function on_check(~, ~)
        problems = check_settings();
        if isempty(problems)
            message('Everything is ready: press Run.', 'Check settings', 'success');
        else
            message(sprintf('Before running:\n- %s', strjoin(problems, '\n- ')), 'Check settings');
        end
    end

    function on_run(~, ~)
        if running
            return
        end
        problems = check_settings();
        if ~isempty(problems)
            message(sprintf('Before running:\n- %s', strjoin(problems, '\n- ')), 'The run cannot start');
            return
        end
        if table_dirty || strlength(table_file) == 0
            if ~save_table(false)
                return
            end
        end
        args = run_options();
        running = true;
        stop_requested = false;
        set_running(true);
        tabs.SelectedTab = tab_run;
        log_area.Value = {char(sprintf('%s  Run started: %s', now_text(), command_text()))};
        progress_table.Data = cell(0,3);
        set_status('Running... (the window is updated after each stage)');
        drawnow;
        try
            run_pipeline(table_file, args{:}, 'progress', @on_progress);
        catch ME
            if isvalid(fig)
                figure(fig);                                                % (to the front, as at the end of a run)
            end
            add_log(sprintf('The run stopped with an error: %s', ME.message));
            message(sprintf('The run stopped with an error:\n%s', ME.message), 'Run');
        end
        running = false;
        if isvalid(fig)
            set_running(false);
        end
        if close_after_run && isvalid(fig)
            delete(fig);
        end
    end

    function stop = on_progress(ev)
        % called by "run_pipeline" at each step of the run (see "report_progress" in run_pipeline.m); returns true to stop
        if ~isvalid(fig)
            stop = true;
            return
        end
        switch ev.type
            case "start"
                run_info.recordings = ev.recordings;
                run_info.groups = unique(ev.groups, 'stable');
                run_info.rec_groups = ev.groups;
                run_info.stages = ev.stages;
                run_info.excel = strings(0);
                init_progress();
                add_log(sprintf('%d recording(s); stages: %s; parameters: %s', numel(ev.recordings), strjoin(STAGE_SHORT(ismember(STAGES, ev.stages)), ' > '), ev.parameters));
            case "stage"
                col = 2 + find(run_info.stages == ev.stage, 1);
                switch ev.state
                    case "running"
                        set_cell(ev.index, col, "running...", COLOR_RUN);
                        set_status(sprintf('Running %s on %s (%d of %d)...', STAGE_SHORT(STAGES == ev.stage), run_info.recordings(ev.index), ev.index, numel(run_info.recordings)));
                    case "done"
                        set_cell(ev.index, col, "done (" + time_short(ev.elapsed) + ")", COLOR_DONE);
                    case "skipped"
                        set_cell(ev.index, col, "skipped", COLOR_SKIP);
                        add_log(sprintf('%s, %s: skipped (%s)', run_info.recordings(ev.index), STAGE_SHORT(STAGES == ev.stage), ev.message));
                    case "error"
                        set_cell(ev.index, col, "ERROR", COLOR_ERROR);
                        add_log(sprintf('%s, %s: ERROR: %s', run_info.recordings(ev.index), STAGE_SHORT(STAGES == ev.stage), ev.message));
                end
            case "recording"
                for col = 3:size(progress_table.Data, 2)
                    set_cell(ev.index, col, upper(ev.state), ifelse(ev.state == "error", COLOR_ERROR, COLOR_SKIP));
                end
                add_log(sprintf('%s: %s (%s)', run_info.recordings(ev.index), ev.state, ev.message));
            case "group"
                row = numel(run_info.recordings) + find(run_info.groups == ev.group, 1);
                col = size(progress_table.Data, 2);
                switch ev.state
                    case "running"
                        set_cell(row, col, "running...", COLOR_RUN);
                        set_status(sprintf('Compiling the group %s...', ev.group));
                    case "done"
                        set_cell(row, col, "done (" + time_short(ev.elapsed) + ")", COLOR_DONE);
                        run_info.excel(end+1) = ev.excel;
                    case "error"
                        set_cell(row, col, "ERROR", COLOR_ERROR);
                        add_log(sprintf('group %s: ERROR: %s', ev.group, ev.message));
                end
            case "end"
                run_info.project = ev.project;
                if ev.stopped
                    txt = sprintf('Stopped by you after %s: the remaining stages did not run.', time_short(ev.elapsed));
                elseif isempty(ev.errors)
                    txt = sprintf('Finished in %s, without errors.', time_short(ev.elapsed));
                else
                    txt = sprintf('Finished in %s, with %d error(s) (in red; the log of each recording has the details).', time_short(ev.elapsed), numel(ev.errors));
                end
                add_log(txt);
                for k=1:numel(ev.errors)
                    add_log("   " + ev.errors(k));
                end
                for k=1:numel(run_info.excel)
                    add_log("   Excel results: " + run_info.excel(k));
                end
                set_status(txt);
                figure(fig);                                                % (new) the window comes to the front at the end of the run: the
                                                                            % plots of the stages and the MATLAB desktop had covered it
                % (revised, user's choice) no message window at the end: the run just ends (status bar, log, and the button
                % "Open the Excel results" of the tab Run)
        end
        drawnow;                                                            % (the button Stop is read here)
        stop = stop_requested;
    end

    function init_progress()
        % one row per recording (its stages) and one row per group (compilation)
        st = run_info.stages;
        rec_st = st(st ~= "Compile_MEAproject");
        n_rec = numel(run_info.recordings);
        has_compile = any(st == "Compile_MEAproject");
        n_col = 2 + numel(rec_st) + has_compile;
        D = repmat({''}, n_rec + has_compile*numel(run_info.groups), n_col);
        D(1:n_rec, 1) = cellstr(run_info.recordings);
        D(1:n_rec, 2) = cellstr(run_info.rec_groups);
        if has_compile
            for k=1:numel(run_info.groups)
                D{n_rec+k, 1} = 'Excel tables of the group';
                D{n_rec+k, 2} = char(run_info.groups(k));
            end
        end
        names = [{'Recording', 'Group'}, cellstr(STAGE_SHORT(ismember(STAGES, rec_st)))];
        if has_compile
            names{end+1} = 'Compile';
            run_info.stages = [rec_st, "Compile_MEAproject"];
        else
            run_info.stages = rec_st;
        end
        removeStyle(progress_table);
        progress_table.Data = D;
        progress_table.ColumnName = names;
        progress_table.ColumnWidth = [{220, 90}, repmat({'auto'}, 1, n_col-2)];
    end

    function set_cell(row, col, txt, color)
        if row > size(progress_table.Data, 1) || col > size(progress_table.Data, 2)
            return
        end
        progress_table.Data{row, col} = char(txt);
        addStyle(progress_table, uistyle('BackgroundColor', color), 'cell', [row col]);
    end

    function add_log(txt)
        v = log_area.Value;
        v{end+1} = sprintf('%s  %s', now_text(), txt);
        log_area.Value = v;
        scroll(log_area, 'bottom');
    end

    function set_running(on)
        run_button.Enable = ~on;
        check_button.Enable = ~on;
        stop_button.Enable = on;
        stop_button.Text = 'Stop';
        rec_table.ColumnEditable = ~on;
        tick_all_button.Enable = ~on;
        untick_all_button.Enable = ~on;
        param_table.ColumnEditable = [false false ~on false false];
        viz_draw_button.Enable = ~on;                                       % (each stage closes the figures: no plot during a run)
        for c = [what_mode, from_stage, only_checks, param_mode, renew_list]
            c.Enable = ~on;
        end
        if ~on
            update_settings_view();
        end
    end

    function on_stop(~, ~)
        if running
            stop_requested = true;
            stop_button.Text = 'Stopping...';
            stop_button.Enable = 'off';
            set_status('Stopping after the stage that is running...');
        end
    end

    function on_open_log(~, ~)
        sel = progress_table.Selection;
        T = rec_table.Data;
        if isempty(sel) || isempty(progress_table.Data)
            message('Click a row of the progress table first.', 'Log');
            return
        end
        row = sel(1);
        n_rec = numel(run_info.recordings);
        if row <= n_rec
            name = run_info.recordings(row);
            k = find(T.file_name == name, 1);
            folder = T.folder(k);
            if ~isfolder(folder)
                folder = fullfile(start_folder(), folder);
            end
            file = fullfile(folder, name + "_pipeline_log.txt");
        else
            file = fullfile(results_folder(), run_info.groups(row - n_rec) + "_pipeline_log.txt");
        end
        if isfile(file)
            open_file(file);
        else
            message(sprintf('There is no log file yet: %s', file), 'Log');
        end
    end

    function on_open_results(~, ~)
        open_file(results_folder());
    end

    function on_open_excel(~, ~)
        files = run_info.excel;
        if isempty(files)
            T = rec_table.Data;
            groups = unique(T.group(T.include & strlength(T.group) > 0));
            files = fullfile(results_folder(), groups + "_Compiled.xlsx");
            files = files(isfile(files));
        end
        if isempty(files)
            message('There are no Excel results yet (stage 8 writes them, one file per group).', 'Excel results');
            return
        end
        for f = files(:)'
            open_file(f);
        end
    end

    function f = tab1_table_file()
        % (new, November 2026) file of the recordings table of tab 1, saved first if needed ("" if there is no table), for tab 4
        f = "";
        if height(rec_table.Data) == 0
            return
        end
        if (table_dirty || strlength(table_file) == 0) && ~save_table(false)
            return
        end
        f = table_file;
    end

    function d = results_folder()
        d = string(project_field.Value);
        if strlength(d) == 0
            d = start_folder();
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   VISUALIZATION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (new, November 2026) tab 5: plots of a recording with "plotMEArecording"

    function on_tab_change(~, ev)
        % the list of recordings of the tab Visualization follows the table of tab 1
        if ev.NewValue == tab_viz
            viz_fill_recordings();
        end
    end

    function viz_fill_recordings()
        % the recordings of the table (rows with a folder and a file name); the one chosen is kept, or the row selected in tab 1 is proposed
        T = rec_table.Data;
        rows = find(strlength(strtrim(T.folder)) > 0 & strlength(strtrim(T.file_name)) > 0)';
        old = viz_rec.Value;
        old_key = "";
        if ~isequal(old, 0) && numel(viz_rec.ItemsData) >= 1 && isfield(viz_rec.UserData, 'keys')
            old_key = viz_rec.UserData.keys(viz_rec.ItemsData == old);
        end
        if isempty(rows)
            viz_rec.Items = {'(open or build a recordings table in tab 1)'};
            viz_rec.ItemsData = 0;
            viz_rec.UserData = struct('keys', "");
            viz_choose_recording();
            return
        end
        keys = row_key(T, rows);
        items = string(rows(:)) + ". " + strtrim(T.file_name(rows(:)));
        label = strtrim(T.organoid_label(rows(:)));
        items(strlength(label) > 0) = items(strlength(label) > 0) + "  (organoid " + label(strlength(label) > 0) + ")";
        viz_rec.Items = cellstr(items);
        viz_rec.ItemsData = rows;
        viz_rec.UserData = struct('keys', keys);
        sel = selected_table_rows();
        if ~isempty(sel) && ismember(sel(1), rows) && ~isequal(sel(1), viz_seen_row)   % (a row selected in tab 1 since the last visit)
            choice = sel(1);
        elseif any(keys == old_key)
            choice = rows(find(keys == old_key, 1));
        else
            choice = rows(1);
        end
        viz_seen_row = [];
        if ~isempty(sel), viz_seen_row = sel(1); end
        viz_rec.Value = choice;
        if ~(any(keys == old_key) && keys(rows == choice) == old_key) || isempty(viz_data)
            viz_choose_recording();
        end
    end

    function viz_choose_recording()
        % datasets of the recording chosen (the analyzed one first), and the photo of its chip (for the outline of the organoids)
        r = viz_rec.Value;
        T = rec_table.Data;
        viz_data = [];
        viz_photo = struct('photo', "", 'side', "", 'n', 0);
        if isequal(r, 0) || r > height(T)
            viz_dataset.Items = {'(choose a recording)'};     viz_dataset.ItemsData = {''};
            viz_info.Text = 'No dataset read yet.';
            viz_fill_times();       viz_update_view();
            return
        end
        folder = row_folder(r);
        name = strtrim(T.file_name(r));
        [files, labels] = dataset_files(folder, name, strtrim(T.organoid_label(r)));
        n = str2double(T.picture_n_organoids(r));
        if isnan(n), n = 0; end
        viz_photo = struct('photo', strtrim(T.picture(r)), 'side', lower(strtrim(T.picture_3D_side(r))), 'n', n);
        if isempty(files)
            viz_dataset.Items = {'(no dataset yet: run the pipeline up to the quality control)'};     viz_dataset.ItemsData = {''};
            viz_info.Text = sprintf('No file %s..._QCmanual.mat in %s.', name, folder);
            viz_fill_times();       viz_update_view();
            return
        end
        viz_dataset.Items = cellstr(labels);
        viz_dataset.ItemsData = cellstr(files);
        viz_dataset.Value = char(files(1));
        viz_load(files(1));
    end

    function viz_load(file)
        % reads a dataset once ("data" of "plotMEArecording"): units, conditions, outline of the cropping; then the preview
        file = string(file);
        if strlength(file) == 0
            return
        end
        set_status(sprintf('Reading %s...', file));     drawnow;
        viz_hold = true;                                                    % (one preview, at the end)
        restore = onCleanup(@() viz_release());
        [~, viz_data] = plotMEArecording(file, "data");
        viz_fill_times();
        viz_info_text();
        clear restore
        viz_update_view();
        set_status(sprintf('Dataset read: %d units, %d recording conditions. Check the preview, then press Draw.', viz_data.n_units, ...
            numel(viz_data.conditions)));
    end

    function viz_release()
        viz_hold = false;
    end

    function D = get_viz_data()
        D = viz_data;                                                       % (for the tests)
    end

    function on_viz_file(~, ~)
        [f, p] = with_dialog(@() uigetfile({'*_QCmanual.mat', 'Dataset of the pipeline (*_QCmanual.mat)'}, 'Choose a dataset', ...
            char(start_folder())));
        if isequal(f, 0)
            return
        end
        file = string(fullfile(p, f));
        viz_photo = struct('photo', "", 'side', "", 'n', 0);               % (the photo is searched in the folder of the file)
        if any(strcmp(viz_dataset.ItemsData, ''))
            viz_dataset.Items = {f};    viz_dataset.ItemsData = {char(file)};
        elseif ~any(strcmp(viz_dataset.ItemsData, char(file)))
            viz_dataset.Items = [viz_dataset.Items, {f}];   viz_dataset.ItemsData = [viz_dataset.ItemsData, {char(file)}];
        end
        viz_dataset.Value = char(file);
        viz_load(file);
    end

    function names = viz_cond_names()
        % names of the recording conditions, as chosen (condition and dose, or file name)
        names = strings(1, 0);
        if ~isempty(viz_data)
            [~, Dn] = plotMEArecording(viz_data, "data", struct('condition_names', string(viz_names.Value)));
            names = [Dn.conditions.name];
        end
    end

    function viz_info_text()
        % (revised, user's choice: without the outline of the organoid) units and conditions of the dataset
        D = viz_data;
        names = viz_cond_names();
        cond = strings(1, numel(D.conditions));
        for c=1:numel(D.conditions)
            cond(c) = sprintf('%s (%g-%g s)', names(c), D.conditions(c).start, D.conditions(c).stop);
        end
        viz_info.Text = sprintf('%d units. Conditions: %s.', D.n_units, strjoin(cond, ', '));
    end

    function viz_names_changed()
        viz_hold = true;
        restore = onCleanup(@() viz_release());
        viz_fill_times();
        if ~isempty(viz_data)
            viz_info_text();
        end
        clear restore
        viz_update_view();
    end

    function viz_fill_times()
        % the whole recording, each recording condition, or a window typed (Custom)
        if isempty(viz_data)
            viz_time.Items = {'Whole recording'};      viz_time.ItemsData = 0;      viz_time.Value = 0;
            viz_from.Value = 0;     viz_to.Value = 0;
            return
        end
        D = viz_data;
        names = viz_cond_names();
        old = viz_time.Value;
        items = sprintf("Whole recording (%g-%g s)", D.t_range(1), D.t_range(2));
        for c=1:numel(D.conditions)
            items(end+1) = sprintf("%s (%g-%g s)", names(c), D.conditions(c).start, D.conditions(c).stop); %#ok<AGROW>
        end
        items(end+1) = "Custom window (type it below)";
        viz_time.Items = cellstr(items);
        viz_time.ItemsData = [0, 1:numel(D.conditions), -1];
        if old == -1 && viz_to.Value > viz_from.Value && viz_from.Value < D.t_range(2) && viz_to.Value > D.t_range(1)
            viz_time.Value = -1;                                            % (the window typed is kept)
        elseif old > 0 && old <= numel(D.conditions)
            viz_time.Value = old;
        else
            viz_time.Value = 0;
        end
        viz_time_choice();
    end

    function viz_time_choice()
        % the fields From / to show the window chosen
        if isempty(viz_data)
            return
        end
        v = viz_time.Value;
        if v == 0
            viz_from.Value = viz_data.t_range(1);     viz_to.Value = viz_data.t_range(2);
        elseif v > 0
            viz_from.Value = viz_data.conditions(v).start;     viz_to.Value = viz_data.conditions(v).stop;
        end
        viz_update_view();
    end

    function viz_custom_time()
        if ~isempty(viz_data)
            viz_time.Value = -1;
        end
        viz_update_view();
    end

    function viz_units_typed()
        if strlength(strtrim(string(viz_units.Value))) > 0
            viz_units_mode.Value = 'list';
        end
        viz_update_view();
    end

    function in = viz_inputs()
        % choices of the tab, as the inputs of "plotMEArecording" (only the ones of the plot chosen)
        type = viz_type.Value;
        in = struct();
        if strcmp(viz_units_mode.Value, 'list')
            in.units = parse_units(viz_units.Value);
        end
        v = viz_time.Value;
        if ~strcmp(type, 'waveforms') && ~isempty(viz_data)
            if v > 0
                in.condition = v;
            elseif v < 0
                if viz_to.Value <= viz_from.Value
                    error('The end of the time window (%g s) must be after its start (%g s)', viz_to.Value, viz_from.Value);
                end
                in.window = [viz_from.Value viz_to.Value];
            end
        end
        lim = strtrim(string(viz_limit.Value));
        if strlength(lim) > 0 && ismember(type, {'map', 'video'})
            x = str2double(lim);
            if ~(x > 0)
                error('The top of the color scale must be a rate above 0 (Hz), or empty for automatic (not "%s")', lim);
            end
        end
        if ismember(type, {'raster', 'bursts', 'video'})
            in.condition_lines = viz_cond.Value;
            in.condition_names = string(viz_names.Value);
            in.time_unit = string(viz_tunit.Value);
        end
        switch type
            case 'raster'
                in.population_rate = viz_pop.Value;
                in.rate_bin = viz_bin.Value/1000;                           % (typed in ms)
            case 'map'
                in.rate_measure = string(viz_measure.Value);
            case 'video'
                in.raster = viz_raster.Value;
                in.population_rate = viz_pop.Value;
                in.rate_bin = viz_bin.Value/1000;
                in.VIDEO_STEP = viz_step.Value;
                in.VIDEO_SPEED = viz_speed.Value;
        end
        if ismember(type, {'map', 'video'})
            if strlength(lim) > 0
                in.color_limit = str2double(lim);
            end
            if ~strcmp(viz_cmap.Value, 'jet')
                in.MAP_COLORMAP = string(viz_cmap.Value);
            end
            if ~strcmp(viz_empty.Value, 'black')
                in.empty_color = string(viz_empty.Value);
            end
            in.outline = viz_outline.Value;
            if viz_outline.Value                                            % (photo of the row of the table, if given)
                if strlength(viz_photo.photo) > 0, in.photo = viz_photo.photo; end
                if strlength(viz_photo.side) > 0, in.photo_text_side = viz_photo.side; end
                if viz_photo.n > 0, in.photo_n_organoids = viz_photo.n; end
            end
        end
        if viz_save.Value && ~strcmp(type, 'video')
            in.output = "auto";
        end
    end

    function viz_update_view()
        % shows the options of the plot chosen (the other rows have no height) and its help, then updates the preview
        type = viz_type.Value;
        is = @(varargin) any(strcmp(type, varargin));
        show = [is('raster', 'video'), is('raster', 'video') && viz_pop.Value, is('video'), is('raster', 'bursts', 'video'), ...
            is('raster', 'bursts', 'video'), is('raster', 'bursts', 'video'), is('map'), is('map', 'video'), is('map', 'video'), ...
            is('map', 'video'), is('map', 'video'), is('video'), is('video'), ~is('video')];
        gv5.RowHeight = num2cell(22*show);
        controls = {viz_pop, [viz_bin_label viz_bin], viz_raster, viz_cond, [viz_names_label viz_names], [viz_tunit_label viz_tunit], ...
            [viz_measure_label viz_measure], [viz_limit_label viz_limit], [viz_cmap_label viz_cmap], [viz_empty_label viz_empty], viz_outline, [viz_step_label viz_step], ...
            [viz_speed_label viz_speed], viz_save};
        for k=1:numel(controls)
            set(controls{k}, 'Visible', show(k));
        end
        viz_pop.Text = ifelse(is('video'), 'Population firing rate above the map', 'Population firing rate (below the raster plot)');
        left4.RowHeight = {'fit', 'fit', ifelse(is('waveforms'), 0, 'fit'), 'fit', 'fit'};
        p_vtime.Visible = ~is('waveforms');
        viz_units.Enable = strcmp(viz_units_mode.Value, 'list');
        viz_draw_button.Text = ifelse(is('video'), 'Make the video', 'Draw');
        switch type
            case 'raster'
                viz_help.Text = 'One row per neuron, one tick per spike; dashed lines: start of each condition. Below: the population firing rate.';
            case 'bursts'
                viz_help.Text = 'The bursts of each neuron in color, with their first and last spike, peak and middle marked (at most 100 neurons).';
            case 'map'
                viz_help.Text = 'Each electrode colored by the firing rate of its neurons in the window; black: no neuron; dashed: outline of the organoid.';
            case 'waveforms'
                viz_help.Text = 'The Kilosort template of each neuron on its peak electrode (at most 36 neurons; the time window is not used).';
            case 'video'
                viz_help.Text = ['MP4 of the electrode map, one frame per step, with the raster plot and/or the population rate above it ' ...
                    'and a cursor that moves with the video; saved next to the dataset (it plays in PowerPoint).'];
        end
        viz_preview();
    end

    function viz_preview()
        % (user's request, November 2026) quick preview of the plot chosen, in the tab, to check the time window and the units:
        % "plotMEArecording" draws it in the panel of the preview (raster plots as images, the first frame of a video); Draw opens it full size
        if viz_hold || ~isvalid(fig)
            return
        end
        delete(viz_preview_box.Children);
        viz_preview_label.FontColor = [0.4 0.4 0.4];
        if isempty(viz_data)
            viz_preview_label.Text = 'The preview appears here once a dataset is read.';
            return
        end
        try
            in = viz_inputs();
        catch ME
            viz_preview_label.Text = ['Not drawn: ' ME.message];
            viz_preview_label.FontColor = [0.8 0 0];
            return
        end
        type = string(viz_type.Value);
        viz_find_outline(in);
        in.parent = viz_preview_box;
        in.preview = true;
        if isfield(in, 'output'), in = rmfield(in, 'output'); end
        try
            [~, res] = plotMEArecording(viz_data, type, in);
        catch ME
            delete(viz_preview_box.Children);
            viz_preview_label.Text = ['Not drawn: ' regexprep(ME.message, '^plotMEArecording: ', '')];
            viz_preview_label.FontColor = [0.8 0 0];
            return
        end
        switch type
            case "waveforms"
                viz_preview_label.Text = sprintf('%d neurons', numel(res.units));
            case "video"
                viz_preview_label.Text = sprintf('First frame of the video: %d frames, %g-%g s (%.1f s of video)', res.n_frames, ...
                    res.window(1), res.window(2), res.n_frames/res.frame_rate);
            otherwise
                viz_preview_label.Text = sprintf('%d neurons, %g-%g s', numel(res.units), res.window(1), res.window(2));
        end
    end

    function viz_find_outline(in)
        % the organoids of the photo of the chip, once per dataset (the map and the video draw their outline)
        if ismember(string(viz_type.Value), ["map", "video"]) && viz_outline.Value && isempty(viz_data.outlines) && ~viz_data.photo_searched
            set_status('Looking for the organoids in the photo of the chip...');     drawnow;
            photo_in = struct('outline', true);
            for fld = intersect(["photo", "photo_text_side", "photo_n_organoids"], string(fieldnames(in)))'
                photo_in.(fld) = in.(fld);
            end
            [~, viz_data] = plotMEArecording(viz_data, "data", photo_in);
            set_status(sprintf('Outline of the organoid: %s.', viz_data.outline_source));
        end
    end

    function [f, res] = on_viz_draw(~, ~)
        % draws the plot chosen in its own window, full size (or makes the video), with "plotMEArecording"
        f = [];     res = [];
        if isempty(viz_data)
            message('Choose a recording with a dataset first (or open a "_QCmanual.mat" file with "Other file...").', 'Visualization');
            return
        end
        type = string(viz_type.Value);
        try
            in = viz_inputs();
        catch ME                                                            % (a choice to correct: a plain message, not an error)
            viz_problem(ME.message);
            return
        end
        viz_find_outline(in);
        if opt.test
            in.visible = false;
        end
        if type == "video"
            in.progress = @viz_progress;
            viz_draw_button.Enable = 'off';
            restore = onCleanup(@() viz_enable_draw());
            set_status('Making the video...');
        else
            set_status('Drawing...');
        end
        drawnow;
        try
            [f, res] = plotMEArecording(viz_data, type, in);
        catch ME                                                            % (e.g., too many units, no templates: the message says
            fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));     %  what to change; the details stay in the
            viz_problem(regexprep(ME.message, '^plotMEArecording: ', ''));           %  Command Window)
            f = [];     res = [];
            return
        end
        if strlength(res.file) > 0
            viz_last = res.file;
            viz_open_button.Enable = 'on';      viz_folder_button.Enable = 'on';
        end
        if type == "video"
            viz_result.Text = sprintf('Video saved: %s (%d frames, %.1f s long)', res.file, res.n_frames, res.n_frames/res.frame_rate);
            set_status('The video is ready (button "Open the saved file").');
        elseif strlength(res.file) > 0
            viz_result.Text = sprintf('Drawn and saved: %s', res.file);
            set_status('The plot is drawn in its own window, and saved.');
        else
            viz_result.Text = sprintf('Drawn: %s of %s, %d neurons, %g-%g s (not saved).', viz_type.Items{strcmp(viz_type.ItemsData, ...
                viz_type.Value)}, viz_data.name, numel(res.units), res.window(1), res.window(2));
            set_status('The plot is drawn in its own window.');
        end
        if ~opt.test && ~isempty(f) && isvalid(f)
            figure(f);                                                      % (in front of the window)
        end
    end

    function viz_open(folder_only)
        % opens the last file saved (or its folder) with the program of the system
        if strlength(viz_last) == 0
            return
        end
        if folder_only
            open_file(fileparts(viz_last));
        else
            open_file(viz_last);
        end
    end

    function viz_problem(txt)
        % a choice of the tab to correct (or a plot that cannot be drawn): status bar and message
        set_status(['Not drawn: ' char(txt)]);
        status_label.FontColor = [0.8 0 0];
        viz_result.Text = ['Not drawn: ' char(txt)];
        message(char(txt), 'Visualization');
    end

    function viz_progress(x)
        if isvalid(fig)
            set_status(sprintf('Making the video: %d%%', round(100*x)));
            drawnow;
        end
    end

    function viz_enable_draw()
        if isvalid(fig)
            viz_draw_button.Enable = ~running;
        end
    end

    function folder = row_folder(r)
        % folder of the row r of the recordings table (relative to the folder of the table when it is not a full path)
        folder = strtrim(rec_table.Data.folder(r));
        if isempty(regexp(char(folder), '^(/|~|[A-Za-z]:[\\/]|\\\\)', 'once')) && strlength(table_file) > 0
            folder = string(fullfile(fileparts(table_file), folder));
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   SMALL HELPERS   %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function d = start_folder()
        % folder of the recordings table; (revised November 2026) before the table is saved, the last folder whose recordings were added
        % (before, the current folder of MATLAB: a table built from a folder was saved there by default, and so were the Excel results);
        % otherwise, the current folder
        d = string(pwd);
        if strlength(table_file) > 0
            d = string(fileparts(table_file));
        elseif strlength(added_folder) > 0
            d = added_folder;
        end
    end

    function set_status(txt)
        status_label.Text = char(txt);
        status_label.FontColor = [0 0 0];
    end

    function message(txt, title, icon)
        if nargin < 3
            icon = 'warning';
        end
        if opt.test
            fprintf('[CellCircuit] %s: %s\n', title, txt);
        else
            uialert(fig, txt, title, 'Icon', icon);
        end
    end

    function on_help(~, ~)
        txt = sprintf(['1. Recordings: open a recordings table, or press "Add recordings from folder..." and fill in the group of each ' ...
            'recording. Tick "include" for the recordings to run (to run only a few: "Untick all", then tick them).\n\n2. Settings: choose what to run (the whole pipeline the first time; ' ...
            'later, e.g., from "7. Analyze" to rerun the analysis), and the parameters (the default values the first time; "the values ' ...
            'saved for each recording" to rerun with the same settings). To change a value only for this run, search it and type the new ' ...
            'value.\n\n3. Run: press Run. The progress table fills in after each stage; the results of each recording are saved next to its ' ...
            'files, and the Excel tables of each group in the results folder. Stop ends the run after the stage that is running.\n\n' ...
            'Every run of a recording saves its parameters and a log (<file_name>_parameters.mat, <file_name>_pipeline_log.txt).\n\n' ...
            'Cropping (tab 1): "picture" is the photo of the chip (empty: the photo in the folder of the recording; "none": no cropping). ' ...
            'Every organoid of the photo becomes its own dataset, named with the labels of "organoid_label" (e.g., 56; 57; 58; empty: ' ...
            'org1, org2...); "picture_organoid" keeps only one of them (its number in the photo). Optional: "picture_3D_side" (side of ' ...
            'the photo with the "3D" text) and "picture_n_organoids"; "crop_ellipse" or "crop_excel" are used instead of the photo.\n\n' ...
            '4. Neuron composition (once the groups are compiled): add the recordings tables, press "Find the groups", tick the groups used to fit ' ...
            'the axes of the neuron scores (at least 20 organoids), add the comparisons, check the activity measures, and press "Run ' ...
            'comparison". "Save settings..." keeps the settings in a file, to run them again.\n\n' ...
            '5. Visualization: choose a recording of the table and its dataset, what to draw, the time window (the whole recording, a ' ...
            'condition, or a window you type) and the units, then press Draw. The plots open in their own windows; the video is saved ' ...
            'next to the dataset.\n\n' ...
            '6. Features plotting (once the groups are compiled): add the groups of tab 1 (or compiled files), tick the groups to plot, ' ...
            'choose the table, the measures and the conditions, and the plot (bar with SEM or SD, violin, or box plot). Draw opens the ' ...
            'plots; Save figure... saves them; Export values... writes the values drawn and their summary to Excel.']);
        if opt.test
            fprintf('%s\n', txt);
        else
            uialert(fig, txt, 'How to use CellCircuit', 'Icon', 'info');
        end
    end

    function on_close(~, ~)
        if running
            if opt.test
                choice = 'Stop and close';
            else
                choice = uiconfirm(fig, 'A run is in progress. Stop it after the stage that is running, and close the window?', 'Close', ...
                    'Options', {'Stop and close', 'Keep running'}, 'DefaultOption', 2, 'CancelOption', 2);
            end
            if strcmp(choice, 'Stop and close')
                stop_requested = true;
                close_after_run = true;
                set_status('Stopping after the stage that is running; the window closes then.');
            end
            return
        end
        if table_dirty && height(rec_table.Data) > 0 && ~opt.test
            choice = uiconfirm(fig, 'The recordings table was modified and not saved.', 'Close', ...
                'Options', {'Save and close', 'Close without saving', 'Cancel'}, 'DefaultOption', 1, 'CancelOption', 3);
            switch choice
                case 'Save and close'
                    if ~save_table(false)
                        return
                    end
                case 'Cancel'
                    return
            end
        end
        delete(fig);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%   LOCAL FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = read_parameter_list(file)

% function to list every parameter of a parameter file: stage (section of the file), function, parameter, value (text), description (the
% comment of its line and of the comment lines that continue it), and the value itself (Default, used to check the values typed)

P = load_parameters(file);
lines = splitlines(string(fileread(file)));
stage = strings(0,1);   fun = strings(0,1);   par = strings(0,1);   desc = strings(0,1);
section = "";
for i=1:numel(lines)
    s = lines(i);
    h = regexp(s, '^%%\s+(\d+[a-z]?\)\s*.*|3Brain.*)$', 'tokens', 'once');
    if ~isempty(h)
        section = strtrim(regexprep(h{1}, '\s*\(.*$', ''));
        continue
    end
    t = regexp(s, '^P\.(\w+)\.([\w.]+)\s*=\s*.*?;\s*(%.*)?$', 'tokens', 'once');
    if isempty(t) || t{1} == "links"
        continue
    end
    d = strtrim(regexprep(string(t{3}), '^%\s*', ''));
    j = i + 1;
    while j <= numel(lines) && ~isempty(regexp(lines(j), '^\s{20,}%', 'once'))   % comment lines that continue the description
        d = d + " " + strtrim(regexprep(lines(j), '^\s*%\s*', ''));
        j = j + 1;
    end
    stage(end+1,1) = section; fun(end+1,1) = t{1}; par(end+1,1) = t{2}; desc(end+1,1) = d; %#ok<AGROW>
end
[~, first] = unique(fun + "." + par, 'last');                               % (a parameter written twice: its last value counts)
first = sort(first);
stage = stage(first);   fun = fun(first);   par = par(first);   desc = desc(first);
n = numel(fun);
value = strings(n,1);   default = cell(n,1);
for k=1:n
    parts = cellstr(split(fun(k) + "." + par(k), "."));
    default{k} = getfield(P, parts{:});
    value(k) = value_text(default{k});
end
L = table(stage, fun, par, value, desc, default, 'VariableNames', {'Stage', 'Function', 'Parameter', 'Value', 'Description', 'Default'});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function txt = value_text(v)

% function to show a value of the parameter file as MATLAB text

if ischar(v)
    txt = "'" + string(v) + "'";
elseif isstring(v) && isscalar(v)
    txt = """" + v + """";
elseif islogical(v) && isscalar(v)
    txt = string(mat2str(v));
elseif isnumeric(v) || islogical(v)
    if isscalar(v)
        txt = string(num2str(v, 8));
    else
        txt = string(mat2str(v, 8));
    end
elseif iscell(v)
    parts = strings(1, numel(v));
    for k=1:numel(v)
        parts(k) = value_text(v{k});
    end
    txt = "{" + strjoin(parts, ", ") + "}";
else
    txt = "(" + class(v) + ")";
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value, ok, why] = parse_value(txt, default)

% function to read a value typed by the user, with the type of the value of the file (number(s), true/false, text, or a list {...})

txt = strtrim(string(txt));
value = [];   ok = true;   why = "";
try
    if islogical(default) && any(lower(txt) == ["true", "false"])
        value = lower(txt) == "true";
    elseif isnumeric(default) || islogical(default)
        value = str2num(char(txt)); %#ok<ST2NM>                         % (numbers, vectors, and expressions such as 1/210)
        if isempty(value) || ~(isnumeric(value) || islogical(value))
            ok = false;   why = "a number is expected (e.g., 5, 0.25, 1/210 or [1 2 3])";
        elseif islogical(default)
            value = logical(value);
        end
    elseif ischar(default) || isstring(default)
        value = regexprep(txt, '^[''"](.*)[''"]$', '$1');               % (the quotes are optional)
        if ischar(default)
            value = char(value);
        end
    elseif iscell(default)
        value = eval(char(txt));
        if ~iscell(value)
            ok = false;   why = "a list is expected, e.g., {'CGP', 'TTX'}";
        end
    else
        ok = false;   why = "this parameter can only be changed in a parameter file";
    end
catch ME
    ok = false;   why = string(ME.message);
    if iscell(default)
        why = "a list is expected, e.g., {'CGP', 'TTX'}";
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function txt = value_code(v)

% function to write a value as MATLAB code, for the equivalent command

if isstring(v) && isscalar(v)
    txt = """" + v + """";
elseif isstring(v)
    txt = "[" + strjoin("""" + v + """", ", ") + "]";
elseif iscell(v)
    parts = strings(1, numel(v));
    for k=1:numel(v)
        if ischar(v{k})
            parts(k) = "'" + string(v{k}) + "'";
        else
            parts(k) = value_code(v{k});                                    % (revised: the values of the changes with all their digits)
        end
    end
    txt = "{" + strjoin(parts, ", ") + "}";
elseif isnumeric(v) && isscalar(v) && v ~= round(v)
    txt = string(num2str(v, 15));                                           % (enough digits to give the same value, e.g., 1/210)
else
    txt = value_text(v);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = find_recordings(d)

% function to find the recordings in a folder and its subfolders: the folders with files of the pipeline ("<name>_raw.mat",
% "<name>_cleaned.mat", "<name>_QCmanual.mat"; "<name>_curated_..." and "<name>_integrated_..." give "<name>"), or with Kilosort output
% (spike_times.npy; the name of the folder is proposed). The folders of the groups (only "_Analyzed.mat" files) are not recordings, nor
% the datasets written by "SelectElectrodes" (November 2026, "cropped_datasets")

files = [dir(fullfile(d, '**', '*_raw.mat')); dir(fullfile(d, '**', '*_cleaned.mat')); dir(fullfile(d, '**', '*_QCmanual.mat'))];
folder = strings(0,1);   name = strings(0,1);
for k=1:numel(files)
    n = string(files(k).name);
    n = regexprep(n, '(_raw|_cleaned|_QCmanual)\.mat$', '');
    n = regexprep(n, '_(curated|integrated)$', '');
    folder(end+1,1) = string(files(k).folder); name(end+1,1) = n; %#ok<AGROW>
end
npy = dir(fullfile(d, '**', 'spike_times.npy'));
for k=1:numel(npy)
    f = string(npy(k).folder);
    if ~ismember(f, folder)
        [~, base, ext] = fileparts(f);
        folder(end+1,1) = f; name(end+1,1) = string(base) + string(ext); %#ok<AGROW>
    end
end
[~, keep] = unique(folder + "|" + name, 'stable');
folder = folder(keep);   name = name(keep);
[~, order] = sort(folder + "|" + name);
folder = folder(order);   name = name(order);
cropped = cropped_datasets(folder, name);                                   % (November 2026) datasets written by "SelectElectrodes" (the
folder = folder(~cropped);   name = name(~cropped);                         % cropping is run again from the row of the recording)
n = numel(folder);
picture = strings(n,1);
for k = 1:n
    picture(k) = recording_picture(folder(k), name(k));                    % (November 2026) photo of the chip in the folder of the recording
end
R = table(true(n,1), folder, name, strings(n,1), strings(n,1), categorical(repmat("Kilosort", n, 1), ["Kilosort", "3Brain"]), ...
    strings(n,1), picture, strings(n,1), strings(n,1), strings(n,1), strings(n,1), strings(n,1), strings(n,1), 'VariableNames', ...
    {'include', 'folder', 'file_name', 'group', 'to_discard', 'source', 'organoid_label', 'picture', 'picture_organoid', 'picture_3D_side', ...
    'picture_n_organoids', 'crop_ellipse', 'crop_excel', 'file_name_3Brain'});
R.source(endsWith(R.file_name, "_3Brain")) = "3Brain";

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = file_type(file)
if endsWith(lower(string(file)), ".csv")
    t = 'text';
else
    t = 'spreadsheet';
end
end

function open_file(f)
% opens a file or a folder with the program of the system (Excel, text editor, Finder / Explorer)
f = char(f);
if ispc
    winopen(f);
elseif ismac
    system(sprintf('open "%s" &', f));
else
    system(sprintf('xdg-open "%s" &', f));
end
end

function v = get_pref(name, default)
v = default;
if ispref('CellCircuit', name)
    v = string(getpref('CellCircuit', name));
end
end

function set_pref(name, value)
setpref('CellCircuit', name, char(value));
end

function txt = now_text()
txt = char(datetime('now', 'Format', 'HH:mm:ss'));
end

function txt = time_short(sec)
if sec < 60
    txt = sprintf('%.0f s', sec);
elseif sec < 3600
    txt = sprintf('%.1f min', sec/60);
else
    txt = sprintf('%.1f h', sec/3600);
end
end

function out = ifelse(c, a, b)
if c
    out = a;
else
    out = b;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [files, labels] = dataset_files(folder, name, organoid_label)

% (new, November 2026) datasets of a recording for the tab Visualization: the files "<name>..._QCmanual.mat" of its folder (after the
% quality control: "<name>", the cropped "<name>_<label>", the analyzed "..._curated" or "..._integrated"); first the ones of the row (its
% organoid, when it has one), the analyzed one first

files = strings(0,1);   labels = strings(0,1);
name = string(name);    organoid_label = string(organoid_label);
if strlength(name) == 0 || ~isfolder(folder)
    return
end
d = dir(fullfile(folder, name + "*_QCmanual.mat"));
stems = regexprep(string({d.name})', '_QCmanual\.mat$', '');
stems = stems(stems == name | startsWith(stems, name + "_"));
if isempty(stems)
    return
end
analyzed = endsWith(stems, ["_curated", "_integrated"]);
if strlength(organoid_label) > 0
    own = stems == name + "_" + organoid_label | startsWith(stems, name + "_" + organoid_label + "_");
else
    own = ismember(stems, name + ["", "_curated", "_integrated"]);
end
[~, order] = sort(2*own + analyzed, 'descend');
stems = stems(order);   analyzed = analyzed(order);
files = string(folder) + filesep + stems + "_QCmanual.mat";
labels = stems;
labels(analyzed) = labels(analyzed) + "   (analyzed)";

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function units = parse_units(txt)

% (new, November 2026) unit numbers typed in the tab Visualization: numbers and ranges ("1:20", or "1-20"), separated by spaces or commas

tokens = unit_tokens(txt);
if isempty(tokens)
    error('Type the numbers of the units (e.g., 1:20 35 40), or choose "All units"');
end
units = [];
for t = tokens
    m = regexp(t, '^(\d+):(\d+)$', 'tokens', 'once');
    if ~isempty(m)
        a = str2double(m{1});   b = str2double(m{2});
        if b < a
            error('The range %s goes down: write it from the smaller number (e.g., %d:%d)', t, b, a);
        end
        units = [units, a:b]; %#ok<AGROW>
    elseif ~isempty(regexp(t, '^\d+$', 'once'))
        units(end+1) = str2double(t); %#ok<AGROW>
    else
        error('"%s" is not a unit number or a range of units (e.g., 1:20)', t);
    end
end

end

function tokens = unit_tokens(txt)
txt = strtrim(string(txt));
txt = regexprep(txt, '\s*[:\-]\s*', ':');                                 % ("1 - 20" or "1 : 20" -> "1:20")
tokens = split(regexprep(txt, '[,;\s]+', ' '), ' ')';
tokens = tokens(strlength(tokens) > 0);
end

function keys = row_key(T, rows)
% (new, November 2026) identity of rows of the recordings table (folder, file name, organoid), to follow a recording when rows move
keys = strtrim(T.folder(rows(:))) + "|" + strtrim(T.file_name(rows(:))) + "|" + strtrim(T.organoid_label(rows(:)));
end
