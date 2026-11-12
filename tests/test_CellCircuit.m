function test_CellCircuit()

% TEST OF THE WINDOW "CellCircuit" (new, November 2026), without data: run it with the main folder on the MATLAB path.
% It checks that (1) every parameter of "common_parameters" exists in "pipeline_parameters.m", (2) the window opens, with the common
% adjustments first and every parameter in the advanced view, (3) the values typed are checked against the type of the value of the file,
% (4) the equivalent command of "run_pipeline" follows the settings, (5) the window reports the problems that would stop a run, and
% (6, November 2026) a table built from a folder is saved in that folder by default.
% A run itself needs recordings: see the "progress" option of "run_pipeline".

root = fileparts(fileparts(mfilename('fullpath')));
P = load_parameters(fullfile(root, 'pipeline_parameters.m'));
C = common_parameters();
for k=1:height(C)
    parts = cellstr(split(C.key(k), "."));
    assert(isfield_path(P, parts), "common_parameters: %s is not in pipeline_parameters.m", C.key(k));
end
fprintf('common_parameters: %d parameters, all in pipeline_parameters.m\n', height(C));

folder = tempname;   mkdir(folder);   cleanup = onCleanup(@() rmdir(folder, 's'));
mkdir(fullfile(folder, 'rec1'));
file = fullfile(folder, 'recordings.csv');
writetable(table(1, "rec1", "rec1_name", "groupA", 'VariableNames', {'include', 'folder', 'file_name', 'group'}), file);

prefs_before = current_prefs();                                             % (the user's own list and settings: the test must not change them)
fig = CellCircuit('test', true, 'table', file);
closer = onCleanup(@() delete(fig));
t = fig.UserData.test;   c = t.controls;
assert(height(c.rec_table.Data) == 1 && c.rec_table.Data.group(1) == "groupA", "the table was not loaded");
assert(height(c.param_table.Data) == height(C), "the common view does not show the common parameters");

% recordings table: the group typed goes to every selected row; Remove removes every selected row
writetable(table([1;1;1], ["rec1";"rec2";"rec3"], ["a";"b";"c"], ["groupA";"groupA";"groupA"], 'VariableNames', {'include', 'folder', 'file_name', 'group'}), ...
    fullfile(folder, 'three.csv'));
t.open_table(fullfile(folder, 'three.csv'));
b_group = findall(fig, 'Type', 'uibutton', 'Text', 'Set group of selected rows');
c.group.Value = 'groupB';   c.rec_table.Selection = [1 3];   b_group.ButtonPushedFcn(b_group, []);
assert(isequal(c.rec_table.Data.group, ["groupB"; "groupA"; "groupB"]), "the group was not set for every selected row");
b_remove = findall(fig, 'Type', 'uibutton', 'Text', 'Remove selected');
c.rec_table.Selection = [2 3];   b_remove.ButtonPushedFcn(b_remove, []);
assert(isequal(c.rec_table.Data.file_name, "a"), "Remove did not remove every selected row");
% "Untick all" / "Tick all" set the column include of every row, and the count next to them follows
t.open_table(fullfile(folder, 'three.csv'));
b_untick = findall(fig, 'Type', 'uibutton', 'Text', 'Untick all');   b_tick = findall(fig, 'Type', 'uibutton', 'Text', 'Tick all');
count = @() ticked_text(fig);
assert(count() == "3 of 3 recordings ticked (column include)", "wrong count: %s", count());
b_untick.ButtonPushedFcn(b_untick, []);
assert(~any(c.rec_table.Data.include) && count() == "0 of 3 recordings ticked (column include)", "Untick all did not untick every row");
assert(any(contains(t.check(), "No recording is included")), "no recording ticked was not reported");
c.rec_table.Data.include(2) = true;   c.rec_table.CellEditCallback(c.rec_table, []);
assert(count() == "1 of 3 recordings ticked (column include)", "the count did not follow a tick: %s", count());
b_tick.ButtonPushedFcn(b_tick, []);
assert(all(c.rec_table.Data.include) && count() == "3 of 3 recordings ticked (column include)", "Tick all did not tick every row");
% other columns of the user are kept when rows are added and when the table is saved
writetable(table([1;1], ["rec1";"rec2"], ["a";"b"], ["g";"g"], ["note 1";"note 2"], 'VariableNames', {'include', 'folder', 'file_name', 'group', 'notes'}), ...
    fullfile(folder, 'notes.xlsx'));
t.open_table(fullfile(folder, 'notes.xlsx'));
b_empty = findall(fig, 'Type', 'uibutton', 'Text', 'Add empty row');   b_empty.ButtonPushedFcn(b_empty, []);
b_save = findall(fig, 'Type', 'uibutton', 'Text', 'Save');   b_save.ButtonPushedFcn(b_save, []);
S = readtable(fullfile(folder, 'notes.xlsx'), 'TextType', 'string');
assert(ismember("notes", string(S.Properties.VariableNames)) && height(S) == 3 && S.notes(2) == "note 2", "the column notes was not kept");
t.open_table(file);

% two levels: the list of common adjustments holds at most 15 parameters (moved from and to the full list)
n0 = height(t.common());
assert(n0 == min(height(C), 15), "the list does not start with the default common adjustments");
t.remove_from_common(C.key(1:2));
assert(height(t.common()) == n0 - 2 && ~any(ismember(C.key(1:2), t.common().key)), "parameters were not removed from the list");
assert(t.add_to_common("NetworkBurstDetection.method3.MIN_ROBUST"), "a parameter could not be added");
assert(height(t.common()) == n0 - 1 && t.common().key(end) == "NetworkBurstDetection.method3.MIN_ROBUST", "the parameter was not added at the end");
if n0 == 15
    assert(~t.add_to_common(["QualityControl_HyperSynch.N_ITERATION", "Analyze_spikeTrain.MIN_N_BURST"]), "the list went over 15");
    assert(height(t.common()) == 14, "a refused addition changed the list");
end
assert(t.set_change('NetworkBurstDetection.method3.MIN_ROBUST', "5"), "a value of an added parameter was refused");
t.set_change('NetworkBurstDetection.method3.MIN_ROBUST', "");
t.default_list();
assert(isequal(t.common().key, C.key(1:n0)), "the default list did not come back");
assert(isequal(prefs_before, current_prefs()), "the test changed the preferences of the user");

% a rate shown as a period: Compile_MEAproject.MIN_NETBURST (1/210 events per s) is shown and typed in seconds in the common list
row = find(c.param_table.Data.parameter == "Compile_MEAproject.MIN_NETBURST");
assert(isscalar(row) && c.param_table.Data.("in the file")(row) == "210", "MIN_NETBURST is not shown as 210 s");
assert(t.set_change_common('Compile_MEAproject.MIN_NETBURST', "300"), "a period was refused");
row = find(c.param_table.Data.parameter == "Compile_MEAproject.MIN_NETBURST");
assert(c.param_table.Data.("new value")(row) == "300", "the new period is not shown in seconds");
assert(contains(t.command(), "'Compile_MEAproject.MIN_NETBURST', 0.00333333333333333"), "the period was not passed as a rate: %s", t.command());
assert(~t.set_change_common('Compile_MEAproject.MIN_NETBURST', "-5"), "a negative period was accepted");
assert(t.set_change('Compile_MEAproject.MIN_NETBURST', "1/200"), "a rate typed in the full list was refused");
row = find(c.param_table.Data.parameter == "Compile_MEAproject.MIN_NETBURST");
assert(c.param_table.Data.("new value")(row) == "200", "a rate typed in the full list is not shown as a period");
t.set_change('Compile_MEAproject.MIN_NETBURST', "");

% values typed by the user
assert(t.set_change('NetworkBurstDetection.method3.THRESHOLD', "6"), "a number was refused");
assert(~t.set_change('NetworkBurstDetection.method3.THRESHOLD', "six"), "a text was accepted for a number");
assert(t.set_change('Compile_MEAproject.MIN_NETBURST', "1/200"), "an expression was refused");
assert(t.set_change('Compile_MEAproject.TREATMENT', "{'BASELINE', 'TTX'}"), "a list was refused");
assert(~t.set_change('Compile_MEAproject.TREATMENT', "BASELINE"), "a text was accepted for a list");

% equivalent command
c.what_mode.Value = 'from';   c.from_stage.Value = 'Analyze_spikeTrain';
cmd = t.command();
assert(contains(cmd, "'from', ""Analyze_spikeTrain""") && contains(cmd, "'NetworkBurstDetection.method3.THRESHOLD', 6") && ...
    contains(cmd, "'Compile_MEAproject.MIN_NETBURST', 0.005"), "wrong command: %s", cmd);
c.param_mode.Value = 'saved';   c.param_mode.ValueChangedFcn(c.param_mode, []);   c.renew.Value = {'cluster_neuronal_events'};
cmd = t.command();
assert(contains(cmd, "'parameters', ""saved""") && contains(cmd, "'renew', ""cluster_neuronal_events"""), "wrong command: %s", cmd);
t.set_change('NetworkBurstDetection.method3.THRESHOLD', "");
assert(~contains(t.command(), "method3.THRESHOLD"), "an undone change is still in the command");
fprintf('equivalent command: %s\n', t.command());

% "Only some stages": one tick box per stage, and the line "This run" shows the stages
c.what_mode.Value = 'only';   c.what_mode.ValueChangedFcn(c.what_mode, []);
for k=1:numel(c.only_checks)
    c.only_checks(k).Value = contains(c.only_checks(k).Text, "7. Analyze");
end
c.only_checks(1).ValueChangedFcn(c.only_checks(1), []);
assert(contains(t.command(), "'stages', ""Analyze_spikeTrain"""), "the ticked stage is not in the command: %s", t.command());
summary = findall(fig, 'Type', 'uilabel', 'Text', 'This run: Analyze');
assert(isscalar(summary), "the line This run does not show the ticked stage");
for k=1:numel(c.only_checks)
    c.only_checks(k).Value = false;
end
c.only_checks(1).ValueChangedFcn(c.only_checks(1), []);
assert(any(contains(t.check(), "No stage")), "no stage ticked was not reported");
c.what_mode.Value = 'from';   c.what_mode.ValueChangedFcn(c.what_mode, []);

% problems that would stop a run
assert(isempty(t.check()), "unexpected problems: %s", strjoin(t.check(), " | "));
c.rec_table.Data.group(1) = "";
assert(any(contains(t.check(), "group")), "a missing group was not reported");
c.what_mode.Value = 'all';
c.npy.Value = '';
if exist('readNPY', 'file') ~= 2
    assert(any(contains(t.check(), "npy-matlab")), "the missing npy-matlab folder was not reported");
end
% opened as the users open it (no output): the buttons must still work (the window variable was cleared, November 2026)
delete(fig);
CellCircuit('test', true, 'table', file);
fig2 = findall(groot, 'Type', 'figure', 'Name', 'CellCircuit - run the pipeline');
assert(isscalar(fig2), "the window opened without output was not found");
fig2.UserData.test.check();
fig2.CloseRequestFcn(fig2, []);                                             % (the button Close of the window)
assert(~isvalid(fig2), "the window opened without output could not be closed");
% (November 2026) a table built from a folder and not saved yet: the folder whose recordings were added is suggested for the table (and so
% for the Excel results, saved next to it), instead of the current folder of MATLAB
project = fullfile(folder, 'project');   mkdir(fullfile(project, 'recX'));
x = 1;   save(fullfile(project, 'recX', 'recX_QCmanual.mat'), 'x');
fig3 = CellCircuit('test', true);
closer3 = onCleanup(@() delete(fig3));
t3 = fig3.UserData.test;
assert(t3.start_folder() == string(pwd), "without a table and without recordings added, the current folder is not suggested");
assert(t3.add_folder(string(project)) == 1, "the recording of the folder was not added");
assert(t3.start_folder() == string(project), "the folder whose recordings were added is not suggested for the table");
t3.open_table(file);
assert(t3.start_folder() == string(fileparts(file)), "the folder of the opened table is not suggested");
fprintf('a table built from a folder is saved there by default\n');
fprintf('test_CellCircuit: all checks passed\n');

end

function txt = ticked_text(fig)
% text of the count of recordings ticked, next to the buttons Tick all / Untick all
labels = findall(fig, 'Type', 'uilabel');
txt = string({labels.Text});
txt = txt(contains(txt, "recordings ticked"));
assert(isscalar(txt), "the count of recordings ticked was not found");
end

function P = current_prefs()
P = [];
if ispref('CellCircuit')
    P = getpref('CellCircuit');
end
end

function tf = isfield_path(S, parts)
tf = true;
for k=1:numel(parts)
    if ~isstruct(S) || ~isfield(S, parts{k})
        tf = false;
        return
    end
    S = S.(parts{k});
end
end
