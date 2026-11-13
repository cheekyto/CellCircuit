function passed = test_recordings_table(code_dir)

% Test of the reading of the recordings table by "run_pipeline" (November 2026), without data. A .csv table whose file names are full of
% "_" and whose optional columns are empty was read with "_" as the delimiter (no column found: "The recordings table must have the
% columns folder, file_name, group"). It checks that such a table is read column by column, as .csv and as .xlsx:
%   1. with every row excluded (include = 0), the run stops because there is no recording to run (the table itself was read)
%   2. a wrong value in a later row (crop_ellipse with 3 numbers) is reported for that row (the values are in their own columns)
%   3. (November 2026) a table saved by the window "CellCircuit" for a few recordings of a long list: the unticked rows without a group
%      or a file name are skipped (they stopped the run: "Row 1 of the recordings table: the columns file_name and group cannot be
%      empty"), and the ticked row is reached (its folder is missing: the run reports it); a ticked row without a group still stops the run
%   4. (November 2026) photo of the chip: "recording_picture" finds the photo in the folder of a recording (not the figures of the
%      pipeline; none chosen when there are several), "make_recordings_table" fills the column picture and leaves out the datasets
%      written by "SelectElectrodes" (not a recording whose name starts with the name of another one); a recording with a picture but
%      no picture_organoid stops the run when "SelectElectrodes" is run (not otherwise); a wrong picture_3D_side, or a missing picture
%      with picture_organoid, stops the run with its row
%
% Example (from a terminal):
%   matlab -batch "addpath('<repo>/tests'); test_recordings_table('<repo>')"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dir:             folder with the code to be tested

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% passed:               =true if all the checks passed (each check is printed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(char(code_dir));
passed = true;
folder = tempname;      mkdir(folder);      cleanup = onCleanup(@() rmdir(folder, 's'));
header = "folder,file_name,group,to_discard,organoid_label,crop_ellipse,crop_excel,file_name_3Brain,source,include";
names = ["101_lineA_100div_wt_drugs", "102_lineA_100div_wt_drugs", "103_lineA_100div_het_drugs"];
rows = strings(numel(names), 1);
for r=1:numel(names)
    rows(r) = sprintf("%s,%s,wt,,,,,,Kilosort,0", fullfile(folder, extractBefore(names(r), 4) + "_rec_100Hz_filter"), names(r));
end

csv1 = fullfile(folder, 'excluded.csv');
writelines([header; rows], csv1);
passed = report("csv, names with '_', empty optional columns, every row excluded: no recording to run", ...
    stops_with(csv1, "No recording to run")) && passed;

rows2 = rows;   rows2(2) = replace(rows2(2), ",,,,,,Kilosort,0", ",,org1,1 2 3,,,Kilosort,1");
csv2 = fullfile(folder, 'wrong_row.csv');
writelines([header; rows2], csv2);
passed = report("csv: a wrong crop_ellipse in row 2 is reported for that row", stops_with(csv2, "Row 2 of the recordings table: crop_ellipse")) ...
    && passed;

xlsx1 = fullfile(folder, 'excluded.xlsx');
writetable(readtable(csv1, 'Delimiter', ',', 'TextType', 'string'), xlsx1);
passed = report("xlsx, the same table: no recording to run", stops_with(xlsx1, "No recording to run")) && passed;

% 3. the table of the window: every row, the unticked ones without a group (and one without a file name), one ticked row
W = table([0;0;1], [string(fullfile(folder,'recA')); string(fullfile(folder,'recB')); string(fullfile(folder,'recC_missing'))], ...
    ["recA_curated"; ""; "recC_curated"], [""; ""; "groupA"], ["";"";""], ["Kilosort";"Kilosort";"Kilosort"], ["";"";""], ["";"";""], ...
    ["";"";""], ["";"";""], 'VariableNames', {'include', 'folder', 'file_name', 'group', 'to_discard', 'source', 'organoid_label', ...
    'crop_ellipse', 'crop_excel', 'file_name_3Brain'});
for ext = [".xlsx" ".csv"]
    file = fullfile(folder, "window" + ext);
    writetable(W, file);
    passed = report(sprintf("%s of the window: unticked rows without group or file name skipped, the ticked row reached", ext), ...
        reaches(file, "recC_missing does not exist")) && passed;
end
W.include(1) = 1;   W.group(1) = "";                                        % a ticked row without a group
file = fullfile(folder, "window_group.xlsx");
writetable(W, file);
passed = report("a ticked row without a group still stops the run, with its row", ...
    stops_with(file, "Row 1 of the recordings table: the columns file_name and group cannot be empty")) && passed;

% 4. photo of the chip (column picture)
recP = fullfile(folder, 'data', 'recP');    recQ = fullfile(folder, 'data', 'recQ');
mkdir(recP);    mkdir(recQ);
imwrite(uint8(zeros(8, 8, 3)), fullfile(recP, 'chip photo.jpg'));
for f = ["recP_org1_crop_map.png", "recP_org1_crop_photo.png", "recP_plot.png", "._chip photo.jpg"]
    writelines("x", fullfile(recP, f));                                     % figures of the pipeline, hidden file of macOS
end
writelines("x", fullfile(recQ, 'a.jpg'));   writelines("x", fullfile(recQ, 'b.HEIC'));
[pic, cand] = recording_picture(recP, "recP");
passed = report("recording_picture: the photo, not the figures of the pipeline", pic == "chip photo.jpg" && isscalar(cand)) && passed;
[pic, cand] = recording_picture(recQ, "recQ");
passed = report("recording_picture: several photos, none chosen", pic == "" && numel(cand) == 2) && passed;
save(fullfile(recP, 'recP_QCmanual.mat'), 'pic');     save(fullfile(recQ, 'recQ_QCmanual.mat'), 'pic');
spikeTimestamp_BeforeCropping = 1;                                          % datasets written by "SelectElectrodes" (and later stages)
save(fullfile(recP, 'recP_orgA_QCmanual.mat'), 'pic', 'spikeTimestamp_BeforeCropping');
save(fullfile(recP, 'recP_orgA_curated_QCmanual.mat'), 'pic');
save(fullfile(recQ, 'recQ_day2_QCmanual.mat'), 'pic');                      % a recording whose name starts with another one (not cropped)
evalc('T = make_recordings_table(fullfile(folder, ''data''), fullfile(folder, ''built.xlsx''), ''group'', "g1")');
passed = report("make_recordings_table: column picture filled (one photo), picture_organoid empty", ...
    isequal(T.picture, ["chip photo.jpg"; ""; ""]) && all(T.picture_organoid == "")) && passed;
passed = report("make_recordings_table: the datasets written by SelectElectrodes left out, not a recording with a longer name", ...
    isequal(T.file_name, ["recP"; "recQ"; "recQ_day2"])) && passed;
P = table(string(recP), "recP", "g1", "orgA", "chip photo.jpg", "", "", 'VariableNames', ...
    {'folder', 'file_name', 'group', 'organoid_label', 'picture', 'picture_organoid', 'picture_3D_side'});
file = fullfile(folder, "picture.xlsx");    writetable(P, file);
passed = report("a picture without picture_organoid stops the run of SelectElectrodes (before it starts)", ...
    stops_with(file, "have a picture but no picture_organoid", "SelectElectrodes")) && passed;
passed = report("... but not the other stages", reaches(file, "recP")) && passed;
P.picture_organoid = "1";   P.picture_3D_side = "up";     writetable(P, file);
passed = report("a wrong picture_3D_side stops the run with its row", stops_with(file, "Row 1 of the recordings table: picture_3D_side")) && passed;
P.picture_3D_side = "top";  P.picture = "missing.jpg";    writetable(P, file);
passed = report("a missing picture (with picture_organoid) stops the run with its row", ...
    stops_with(file, "Row 1 of the recordings table: the picture")) && passed;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = reaches(table_file, expected)

% =true if "run_pipeline" reads the table and reaches the recordings: its messages contain the expected text, whether it ends with an
% error (the recordings that failed) or not, and the error is not about the table

ok = false;
try
    out = evalc('run_pipeline(table_file, ''stages'', "IntegrateDataSets")');
    ok = contains(out, expected);
catch ME
    ok = ~contains(ME.message, "recordings table") && ~contains(ME.message, "No recording to run");
    if ~ok
        fprintf("        (error: %s)\n", ME.message);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = stops_with(table_file, expected, stage)

% =true if "run_pipeline" stops with an error that contains the expected text (the stages are never reached)

ok = false;
try
    if nargin < 3, stage = "IntegrateDataSets"; end
    evalc('run_pipeline(table_file, ''stages'', stage)');
catch ME
    ok = contains(ME.message, expected);
    if ~ok
        fprintf("        (error: %s)\n", ME.message);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = report(name, ok)

if ok
    fprintf("  PASS  %s\n", name);
else
    fprintf("  FAIL  %s\n", name);
end

end
