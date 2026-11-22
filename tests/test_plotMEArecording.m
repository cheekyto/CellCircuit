function passed = test_plotMEArecording(code_dir)

% Function-level test of "plotMEArecording" (November 2026) on a synthetic dataset (6 units, a baseline and a TTX condition), whose
% spike trains give known answers. It checks that
%   1. the dataset is read: units, electrodes (positionXY + 1), conditions named by condition and dose as the compiled tables
%      ("label_condition"), or by their file names without the parts they all share (the baseline found by its name, also when it is not
%      listed first; two conditions with the same name numbered), window of the recording; the structure of "data" = the file
%   2. raster: every spike of the units and window chosen is drawn, in the rows of the units (4-5 of them labelled with their numbers;
%      axis "Neuron"), with the start of each condition and its name above the plot (on a second row when it would overlap), no title;
%      population firing rate = spikes of all the units / number of units / bin (Hz), as blue bars (one per bin), with the bin chosen; the time axis in s, min or ms; the window
%      of a condition (by name or number) and a window given; the preview in a panel (raster plot as an image, nothing saved)
%   3. bursts: the planted bursts are found and marked (first spike), none in a window without bursts; at most MAX_BURST_UNITS units
%   4. map: rate of each electrode (sum of its units for "mean", largest for "median"), x to the right and y upwards, NaN without
%      units, automatic or given color limit, BORDER
%   5. outline of the organoid: from the cropping (outline from the photo -also for the analyzed dataset made from the cropped one-, or
%      ellipse), none without a photo, a photo that cannot be read
%      (warning, no outline), no outline when not asked
%   6. waveforms: the Kilosort templates of the units chosen (in their order), and a clear error without them
%   7. video: one frame per step, rate of each electrode in each step = direct count, MP4 file with that number of frames and the frame
%      rate VIDEO_SPEED / VIDEO_STEP, with or without the raster plot and the population rate above the map, progress reported
%   8. the plots are saved next to the dataset ("auto"), and wrong choices stop with a clear message
% No reference recording is needed
%
% Example (from a terminal):
%   matlab -batch "addpath('<repo>/tests'); test_plotMEArecording('<repo>')"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dir:             folder with the code to be tested

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% passed:               =true if all the checks passed (each check is printed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(char(code_dir));
passed = true;
root = tempname;    mkdir(root);
cleanup = onCleanup(@() rmdir(root, 's'));
folder = fullfile(root, 'rec_D3');      mkdir(folder);
hide = struct('visible', false);

% synthetic dataset: spike times in sec (saved in msec), electrodes 0-based (positionXY)
burst_starts = [11 23 35];
u3 = sort([0.5:3:89.9, reshape(burst_starts + (0:5)'*0.010, 1, [])]);       % isolated spikes every 3 s + 3 bursts of 6 spikes 10 ms apart
trains = {0.1:0.2:89.95, 0.25:0.5:89.95, u3, 0.7:1:59.95, 40.05:0.1:49.96, []};
pos = {[10 20], [10 20], [0 63], [63 0], [30 30], [5 5]};                    % units 1 and 2 on the same electrode; 3 top-left; 4 bottom-right
clusters = [3 5 8 11 13 17];
st = struct('timeStamp', {}, 'nUnits', {}, 'positionXY', {}, 'cluster', {}, 't_step', {});
for k=1:6
    st(k) = struct('timeStamp', trains{k}(:)*1000, 'nUnits', numel(trains{k}), 'positionXY', pos{k}, 'cluster', clusters(k), 't_step', 0.05);
end
spikeTimestamp = st;
rec_conditions = {'rec_D3_ttx_00.brw', [60 90]; 'rec_D3.brw', [0 60]};     % (TTX listed first: the baseline is found by its file name)
rec_conditions = rec_conditions([2 1], :);
file = fullfile(folder, 'rec_D3_curated_QCmanual.mat');
save(file, 'spikeTimestamp', 'rec_conditions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Dataset

[f, D] = plotMEArecording(file, "data");
ok = isempty(f) && D.n_units == 6 && D.name == "rec_D3_curated" && D.recording == "rec_D3" && isequal([D.conditions.name], ["Baseline", "+TTX"]) ...
    && isequal([D.conditions.file_label], ["Baseline", "ttx_00"]) && isequal([D.conditions.start; D.conditions.stop], [0 60; 60 90]) ...
    && D.base_cond == 1 && isequal(D.t_range, [0 90]) && isequal(D.position(3,:), [1 64]) && isequal(D.position(4,:), [64 1]);
passed = report("data: units, electrodes (positionXY + 1), conditions named by condition and dose (as compiled) or by file name, window", ok) && passed;
rec_conditions = {'rec_D3_ttx_00.brw', [60 90]; 'rec_D3.brw', [0 60]; 'rec_D3_ttx_00.brw', [90 120]};
save(fullfile(root, 'order_QCmanual.mat'), 'spikeTimestamp', 'rec_conditions');
[~, D2] = plotMEArecording(fullfile(root, 'order_QCmanual.mat'), "data");
[~, D2f] = plotMEArecording(fullfile(root, 'order_QCmanual.mat'), "data", struct('condition_names', "file"));
ok = isequal([D2.conditions.name], ["+TTX (1)", "Baseline", "+TTX (3)"]) && isequal([D2f.conditions.name], ["ttx_00 (1)", "Baseline", "ttx_00 (3)"]) ...
    && D2.base_cond == 2 && isequal(D2.t_range, [0 120]);
passed = report("data: the baseline found by its file name when it is not first; two conditions with the same name numbered", ok) && passed;
rec_conditions = {'r.brw', [0 10]; 'r_100div_gabazine_20uM_00.brw', [10 20]; 'r_100div_D-AP5_100uM_DNQX_5uM_00.brw', [20 30]; 'r_100div_TTX_2uM_00.brw', [30 40]};
save(fullfile(root, 'drugs_QCmanual.mat'), 'spikeTimestamp', 'rec_conditions');
[~, D3] = plotMEArecording(fullfile(root, 'drugs_QCmanual.mat'), "data");
ok = isequal([D3.conditions.treatment_label], ["Baseline", "+Gabazine 20uM", "+DNQX+D-AP5 100uM5uM", "+TTX 2uM"]) && ...
    isequal([D3.conditions.file_label], ["Baseline", "gabazine_20uM", "D-AP5_100uM_DNQX_5uM", "TTX_2uM"]);
for c=2:4                                                                   % (the same labels as "Analyze_spikeTrain", by "label_condition")
    [nm, ds] = label_condition(rec_conditions, c, 1, {'CGP', 'Gabazine', 'DNQX', 'D-AP5', 'TTX'}, 0.6);
    ok = ok && D3.conditions(c).treatment_label == string(nm) + " " + string(ds);
end
passed = report("data: condition and dose = label_condition with TREATMENT of pipeline_parameters.m; file names without their shared parts", ok) && passed;
rec_conditions = {'rec_D3.brw', [0 60]; 'rec_D3_ttx_00.brw', [60 90]};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Raster plot

[f, R] = plotMEArecording(file, "raster", hide);
n_all = sum(cellfun(@(t) sum(t >= 0 & t < 90), trains));
axs = findobj(f, 'Type', 'axes');
lines_c = findobj(f, 'Tag', 'condition');
names_c = findobj(f, 'Tag', 'condition_name');
rate_bars = findobj(f, 'Type', 'histogram');                                % (revised November 2026: bars, one per bin)
edges = 0:0.2:90;
direct = histcounts([trains{:}], edges)/6/0.2;                              % (spikes of the 6 units / number of units / bin, user's choice)
ax_r = axs(arrayfun(@(a) contains(string(a.YLabel.String), "Neuron"), axs));
ax_p = axs(arrayfun(@(a) contains(string(a.YLabel.String), "Population"), axs));
ok = R.n_spikes == n_all && numel(axs) == 2 && max(abs(R.rate - direct)) < 1e-9 && numel(R.rate_time) == 450 && ...
    numel(lines_c) == 2 && all([lines_c.Value] == 60) && numel(names_c) == 2 && isequal(sort(string({names_c.String})), ["+TTX", "Baseline"]) && ...
    all(arrayfun(@(h) h.Position(2) == ax_r.YLim(2) && strcmp(h.VerticalAlignment, 'bottom'), names_c)) && ...
    isscalar(ax_r) && isscalar(ax_p) && string(ax_p.YLabel.String) == "Population rate (Hz)" && isscalar(rate_bars) && isequal(rate_bars.FaceColor, [0 0.35 0.80]) && ...
    max(abs(rate_bars.BinEdges - edges)) < 1e-9 && max(abs(rate_bars.BinCounts - direct)) < 1e-9 && ...
    all(arrayfun(@(a) isempty(a.Title.String), axs)) && numel(ax_r.YTick) <= 5;
passed = report("raster, whole recording: every spike, population rate = direct count (Hz, blue bars, one per bin), names of the conditions above the plot, no title", ok) && passed;
close(f);
[f, R] = plotMEArecording(D, "raster", struct('visible', false, 'units', [3 1], 'condition', "+TTX", 'population_rate', false));
ax = findobj(f, 'Type', 'axes');
line_r = findobj(ax, 'Type', 'line');
n_exp = sum(trains{3} >= 60 & trains{3} < 90) + sum(trains{1} >= 60 & trains{1} < 90);
rows3 = line_r.YData(line_r.XData >= 60);
ok = numel(ax) == 1 && isequal(R.window, [60 90]) && R.n_spikes == n_exp && isequal(R.units, [3 1]) && ...
    isequal(string(ax.YTickLabel(:))', ["3", "1"]) && sum(rows3 == 1) == sum(trains{3} >= 60 & trains{3} < 90);
passed = report("raster: units chosen in the order given (rows labelled 3, 1), window of a condition by name, no population rate", ok) && passed;
close(f);
[f, R2] = plotMEArecording(D, "raster", struct('visible', false, 'units', [3 1], 'condition', 2, 'population_rate', false));
close(f);
[f, R2f] = plotMEArecording(D, "raster", struct('visible', false, 'units', [3 1], 'condition', "ttx_00", 'condition_names', "file", 'population_rate', false));
names_f = string(get(findobj(f, 'Tag', 'condition_name'), 'String'));
close(f);
[f, R3] = plotMEArecording(file, "raster", struct('visible', false, 'window', [10 20]));
close(f);
ok = isequal(R2.window, R.window) && R2.n_spikes == R.n_spikes && isequal(R2f.window, R.window) && isequal(names_f, "ttx_00") && isequal(R3.window, [10 20]) && R3.n_spikes == sum(cellfun(@(t) sum(t >= 10 & t < 20), trains));
passed = report("raster: condition by number = by name (also by file name); a window given; the structure of ""data"" = the file", ok) && passed;
% bin of the population rate, unit of the time axis
[f, R4] = plotMEArecording(D, "raster", struct('visible', false, 'rate_bin', 1, 'time_unit', "min"));
axs = findobj(f, 'Type', 'axes');
line_r = findobj(f, 'Type', 'line');
lines_c = findobj(f, 'Tag', 'condition');
ok = numel(R4.rate) == 90 && max(abs(R4.rate - histcounts([trains{:}], 0:1:90)/6)) < 1e-9 && max(abs(R4.rate_time - (0.5:1:89.5))) < 1e-9 && ...
    all(arrayfun(@(a) max(abs(a.XLim - [0 1.5])) < 1e-12, axs)) && abs(max(line_r.XData) - max([trains{:}])/60) < 1e-12 && ...
    all(abs([lines_c.Value] - 1) < 1e-12) && any(arrayfun(@(a) string(a.XLabel.String) == "Time (min)", axs));
close(f);
[f, ~] = plotMEArecording(D, "raster", struct('visible', false, 'window', [10 11], 'time_unit', "ms"));
ax = findobj(f, 'Type', 'axes');
ok = ok && all(arrayfun(@(a) isequal(a.XLim, [10000 11000]), ax)) && any(arrayfun(@(a) string(a.XLabel.String) == "Time (ms)", ax));
close(f);
passed = report("raster: bin of the population rate chosen (1 s), time axis in min or ms (window still in sec)", ok) && passed;
passed = report("raster: a wrong time unit or bin stops with a message", ...
    throws(@() plotMEArecording(D, "raster", struct('visible', false, 'time_unit', "h")), "time_unit") && ...
    throws(@() plotMEArecording(D, "raster", struct('visible', false, 'rate_bin', 0)), "rate_bin")) && passed;
% 4-5 neuron numbers on the rows; names of conditions that would overlap on two rows
st34 = repmat(st(1), 1, 34);
spikeTimestamp = st34;
rec_conditions = {'r.brw', [0 60]; 'r_ttx_00.brw', [60 61]; 'r_cgp_00.brw', [61 90]};
save(fullfile(root, 'many_QCmanual.mat'), 'spikeTimestamp', 'rec_conditions');
spikeTimestamp = st;    rec_conditions = {'rec_D3.brw', [0 60]; 'rec_D3_ttx_00.brw', [60 90]};
[f, ~] = plotMEArecording(fullfile(root, 'many_QCmanual.mat'), "raster", hide);
ax_r = findobj(f, 'Type', 'axes', '-function', @(a) contains(string(a.YLabel.String), "Neuron"));
names_c = findobj(f, 'Tag', 'condition_name');
two_rows = arrayfun(@(h) iscell(h.String) && numel(h.String) == 2, names_c);
ok = isequal(string(ax_r.YTickLabel(:))', ["1", "10", "20", "30"]) && string(ax_r.YLabel.String) == "Neuron" && numel(names_c) == 3 && sum(two_rows) == 1;
close(f);
[f, ~] = plotMEArecording(fullfile(root, 'many_QCmanual.mat'), "raster", struct('visible', false, 'units', [2 9 4 30 31 7 12 3]));
ax_r = findobj(f, 'Type', 'axes', '-function', @(a) contains(string(a.YLabel.String), "Neuron"));
ok = ok && numel(ax_r.YTick) == 5 && ax_r.YTick(1) == 1 && ax_r.YTick(end) == 8;
close(f);
passed = report("raster: 4-5 neuron numbers (1 10 20 30 for 34 neurons); a name that would overlap the one before goes on a second row", ok) && passed;
% preview (window CellCircuit): drawn in a panel, the raster plot as an image
u = uifigure('Visible', 'off');     p = uipanel(u);
[fp, P] = plotMEArecording(D, "raster", struct('parent', p, 'preview', true, 'units', [3 1], 'window', [10 40], 'output', "auto"));
img = findobj(p, 'Type', 'image');
ok = fp == p && isscalar(img) && isempty(findobj(p, 'Type', 'line')) && P.n_spikes == sum(trains{3} >= 10 & trains{3} < 40) + ...
    sum(trains{1} >= 10 & trains{1} < 40) && strlength(P.file) == 0 && ~isfile(fullfile(folder, 'rec_D3_curated_raster_10-40s.png')) && ...
    nnz(any(img.CData(1,:,:) == 0, 3)) == numel(unique(floor((trains{3}(trains{3} >= 10 & trains{3} < 40) - 10)/30*900)));
[~, V0] = plotMEArecording(D, "video", struct('parent', p, 'preview', true, 'window', [40 44]));
ok = ok && strlength(V0.file) == 0 && V0.n_frames == 10 && ~isfile(fullfile(folder, 'rec_D3_curated_video_40-44s.mp4')) && ...
    numel(findobj(p, 'Type', 'image')) == 2;                                % (first frame of the map, and the raster plot as an image)
delete(u);
passed = report("preview: drawn in the panel given, the raster plot as an image (a dark column per bin with spikes), nothing saved, video not written", ok) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Bursts

[f, B] = plotMEArecording(D, "bursts", struct('visible', false, 'units', [1 3], 'condition', "Baseline"));
first = findobj(f, 'Type', 'line', 'DisplayName', 'first spike');
ok = B.n_bursts(2) == 3 && max(abs(B.bursts{2}(:,1)' - burst_starts)) < 1e-9 && max(abs(B.bursts{2}(:,2)' - (burst_starts + 0.05))) < 1e-9 ...
    && B.n_bursts(1) == 0 && numel(first.XData) == 3 && all(first.YData == 2.3);
passed = report("bursts: the 3 planted bursts found and marked in the row of their unit, none in a regular train", ok) && passed;
close(f);
[f, B] = plotMEArecording(D, "bursts", struct('visible', false, 'units', 3, 'window', [60 90]));
close(f);
passed = report("bursts: none in a window without bursts", B.n_bursts == 0) && passed;
passed = report("bursts: more units than MAX_BURST_UNITS stop with a message", ...
    throws(@() plotMEArecording(D, "bursts", struct('visible', false, 'MAX_BURST_UNITS', 2)), "at most 2 units")) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Electrode map

[f, M] = plotMEArecording(D, "map", struct('visible', false, 'condition', "Baseline"));
rate = @(k) sum(trains{k} >= 0 & trains{k} < 60)/60;
im = findobj(f, 'Type', 'image');
ax = findobj(f, 'Type', 'axes');
ok = abs(M.map(21, 11) - (rate(1) + rate(2))) < 1e-9 && abs(M.map(64, 1) - rate(3)) < 1e-9 && abs(M.map(1, 64) - rate(4)) < 1e-9 && ...
    abs(M.map(31, 31) - rate(5)) < 1e-9 && M.map(6, 6) == 0 && nnz(~isnan(M.map)) == 5 && strcmp(ax.YDir, 'normal') && isempty(ax.Title.String) && ...
    isequaln(im.CData, M.map) && isequal(im.AlphaData, ~isnan(M.map)) && ...
    abs(M.color_limit - prctile(M.map(~isnan(M.map)), 99)) < 1e-9 && isequal(ax.CLim, [0 M.color_limit]);
passed = report("map: rate of each electrode (sum of its units), y upwards, NaN without units, automatic color limit, no title", ok) && passed;
close(f);
[f, M2] = plotMEArecording(D, "map", struct('visible', false, 'condition', "Baseline", 'rate_measure', "median", 'color_limit', 3, 'BORDER', 1));
ax = findobj(f, 'Type', 'axes');
ok = abs(M2.map(21, 11) - 1/0.2) < 1e-6 && isnan(M2.map(64, 1)) && isnan(M2.map(1, 64)) && isequal(ax.CLim, [0 3]) && M2.map(6, 6) == 0;
passed = report("map: ""median"" = largest 1/median interval of the units of the electrode, color limit given, BORDER", ok) && passed;
close(f);
[f, ~] = plotMEArecording(D, "map", struct('visible', false, 'MAP_COLORMAP', "hot"));
ax = findobj(f, 'Type', 'axes');
ok = isequal(ax.Colormap, hot(256)) && throws(@() plotMEArecording(D, "map", struct('visible', false, 'MAP_COLORMAP', "rainbowz")), "name of a colormap");
close(f);
passed = report("map: colors chosen (hot), a wrong name stops with a message", ok) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Outline of the organoid

[f, ~] = plotMEArecording(D, "map", hide);
[~, D3] = plotMEArecording(D, "data", struct('outline', true));
ok = isempty(findobj(f, 'Tag', 'outline')) && D3.outline_source == "none (no photo of the chip in the folder)" && D3.photo_searched;
passed = report("outline: none without a cropping or a photo of the chip", ok) && passed;
close(f);
crop_outline = [20 20; 40 22; 38 44; 22 40];
save(fullfile(folder, 'rec_D3_org1_QCmanual.mat'), 'spikeTimestamp', 'rec_conditions', 'crop_outline');
[f, ~] = plotMEArecording(fullfile(folder, 'rec_D3_org1_QCmanual.mat'), "map", hide);
h = findobj(f, 'Tag', 'outline');
ok = isscalar(h) && isequal(h.XData(:), crop_outline([1:end 1], 1)) && isequal(h.YData(:), crop_outline([1:end 1], 2));
close(f);
save(fullfile(folder, 'rec_D3_org1_curated_QCmanual.mat'), 'spikeTimestamp', 'rec_conditions');      % (analyzed: without the outline)
[~, D4] = plotMEArecording(fullfile(folder, 'rec_D3_org1_curated_QCmanual.mat'), "data");
ok = ok && isequal(D4.outlines, {crop_outline}) && D4.outline_source == "cropping (photo)";
passed = report("outline: the one saved by the cropping from the photo, also for the analyzed dataset made from it", ok) && passed;
[f, ~] = plotMEArecording(fullfile(folder, 'rec_D3_org1_QCmanual.mat'), "map", struct('visible', false, 'empty_color', "white"));
ax = findobj(f, 'Type', 'axes');      h = findobj(f, 'Tag', 'outline');
ok = isequal(ax.Color, [1 1 1]) && isequal(h.Color, [0 0 0]);
close(f);
u = uifigure('Visible', 'off');     p = uipanel(u);
plotMEArecording(fullfile(folder, 'rec_D3_org1_QCmanual.mat'), "video", struct('parent', p, 'preview', true, 'window', [40 44], 'empty_color', "white"));
img = findobj(p, 'Type', 'image');
ax_map = img(arrayfun(@(i) size(i.CData, 1) == 64 && size(i.CData, 2) == 64 && ismatrix(i.CData), img)).Parent;
ok = ok && isequal(ax_map.Color, [1 1 1]) && isequal(findobj(p, 'Tag', 'outline').Color, [0 0 0]) && ...
    throws(@() plotMEArecording(D, "map", struct('visible', false, 'empty_color', "gray")), "black"" or ""white");
delete(u);
passed = report("map and video: electrodes without units in white (outline in black), a wrong color stops with a message", ok) && passed;
ellipse_cropped = [32 30 10 5];
save(fullfile(folder, 'rec_D3_org2_QCmanual.mat'), 'spikeTimestamp', 'rec_conditions', 'ellipse_cropped');
[f, ~] = plotMEArecording(fullfile(folder, 'rec_D3_org2_QCmanual.mat'), "map", hide);
h = findobj(f, 'Tag', 'outline');
ok = isscalar(h) && abs(min(h.XData) - 22) < 1e-9 && abs(max(h.XData) - 42) < 1e-9 && abs(max(h.YData) - 35) < 1e-9;
close(f);
[f, ~] = plotMEArecording(fullfile(folder, 'rec_D3_org2_QCmanual.mat'), "map", struct('visible', false, 'outline', false));
ok = ok && isempty(findobj(f, 'Tag', 'outline'));
passed = report("outline: the ellipse of the cropping; none when not asked", ok) && passed;
close(f);
fid = fopen(fullfile(root, 'not_a_photo.jpg'), 'w');    fprintf(fid, 'text');    fclose(fid);
warning('off', 'plotMEArecording:photo');
[f, M3] = plotMEArecording(file, "map", struct('visible', false, 'photo', fullfile(root, 'not_a_photo.jpg')));
warning('on', 'plotMEArecording:photo');
ok = isempty(findobj(f, 'Tag', 'outline')) && M3.outline_source == "none (photo not read)";
passed = report("outline: a photo that cannot be read gives a warning and no outline", ok) && passed;
close(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Waveforms

passed = report("waveforms: a clear error without the Kilosort templates", ...
    throws(@() plotMEArecording(D, "waveforms", struct('visible', false, 'units', [2 4])), "Kilosort templates")) && passed;
rows = cell(1, 6);
fid = fopen(fullfile(folder, 'metadata.csv'), 'w');
fprintf(fid, ',unitType,cluster_id,new_cluster_id,waveform_template,peak_channel\n');
for k=1:6
    rows{k} = single([zeros(1, 21) -exp(-((1:61) - 21 - k).^2/20)*k]);
    fprintf(fid, '%d,%d,%d,%d,"[%s]",%d\n', k-1, 1 + (k == 4), clusters(k), clusters(k), sprintf(' %.8g', rows{k}), k);
end
fclose(fid);
waveforms = rows;
save(fullfile(folder, 'waveforms.mat'), 'waveforms');
[f, W] = plotMEArecording(D, "waveforms", struct('visible', false, 'units', [4 2]));
ax = findobj(f, 'Type', 'axes');
titles = arrayfun(@(a) string(a.Title.String), ax);
ok = isequal(W.templates, double(vertcat(rows{[4 2]}))) && isequal(W.unit_type(:)', [2 1]) && abs(W.sampling_rate - 20000) < 1e-6 && ...
    numel(ax) == 2 && any(startsWith(titles, "Unit 4 (x 64, y 1, MUA)")) && any(startsWith(titles, "Unit 2 (x 11, y 21, good)"));
passed = report("waveforms: the Kilosort template of each unit chosen, in its order, with its electrode and type", ok) && passed;
close(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. Video

done = [];
[f, V] = plotMEArecording(D, "video", struct('visible', false, 'window', [40 44], 'progress', @(x) record_progress(x)));
edges = 40:0.4:44;
direct = nan(64, 64, 10);
for k=1:6
    p = pos{k} + 1;
    c = histcounts(trains{k}, edges)/0.4;
    old = squeeze(direct(p(2), p(1), :))';     old(isnan(old)) = 0;
    direct(p(2), p(1), :) = reshape(old + c, 1, 1, []);
end
reader = VideoReader(char(V.file));
ok = isempty(f) && V.n_frames == 10 && V.frame_rate == 10 && isfile(V.file) && V.file == string(fullfile(folder, 'rec_D3_curated_video_40-44s.mp4')) && ...
    max(abs(double(V.frames(~isnan(direct))) - direct(~isnan(direct)))) < 1e-4 && isequal(isnan(V.frames), isnan(direct)) && ...
    all(abs(V.frames(31, 31, :) - 10) < 1e-4) && max(abs(V.frame_time - (40.2:0.4:43.8))) < 1e-9 && reader.NumFrames == 10 && ...
    abs(reader.FrameRate - 10) < 1e-6 && mod(reader.Width, 2) == 0 && mod(reader.Height, 2) == 0 && done(end) == 1;
passed = report("video: one frame per step, rates = direct count, MP4 with 10 frames at 10 frames/s next to the dataset, progress", ok) && passed;
first_frame = read(reader, 1);
ok = mean(first_frame, 'all') < 128;                                        % (black background)
height_with_plots = reader.Height;
clear reader
out = fullfile(root, 'only_map.mp4');
[~, V2] = plotMEArecording(D, "video", struct('visible', false, 'window', [40 41.2], 'raster', false, 'population_rate', false, ...
    'output', out, 'VIDEO_STEP', 0.2, 'VIDEO_SPEED', 2));
reader = VideoReader(out);
ok = ok && V2.n_frames == 6 && V2.frame_rate == 10 && reader.NumFrames == 6 && reader.Height < height_with_plots;
clear reader
passed = report("video: the map alone, file and step given, frame rate = speed / step", ok) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8. Saving, and wrong choices

[f, R] = plotMEArecording(D, "raster", struct('visible', false, 'output', "auto"));
close(f);
ok = R.file == string(fullfile(folder, 'rec_D3_curated_raster_0-90s.png')) && isfile(R.file);
[f, R] = plotMEArecording(D, "map", struct('visible', false, 'output', fullfile(root, 'map.fig')));
close(f);
g = openfig(char(R.file), 'invisible');
ok = ok && isfile(R.file) && ~isempty(findobj(g, 'Type', 'image'));
close(g);
passed = report("the plots are saved next to the dataset (auto, PNG) or in the file given (.fig)", ok) && passed;
ok = throws(@() plotMEArecording(D, "raster", struct('visible', false, 'units', [1 7])), "units 1 to 6") && ...
    throws(@() plotMEArecording(D, "raster", struct('visible', false, 'window', [100 120])), "outside the recording") && ...
    throws(@() plotMEArecording(D, "raster", struct('visible', false, 'condition', "XYZ")), "no condition ""XYZ""") && ...
    throws(@() plotMEArecording(D, "histogram"), "the plot must be") && ...
    throws(@() plotMEArecording(D, "map", struct('visible', false, 'rate_measure', "max")), "rate_measure") && ...
    throws(@() plotMEArecording(D, "raster", struct('typo', 1)), "not defined");
passed = report("wrong choices stop with a clear message (unit, window, condition, plot, rate measure, input name)", ok) && passed;

close all
if passed
    fprintf("\nALL CHECKS PASSED\n");
else
    fprintf("\nSOME CHECKS FAILED\n");
end

    function record_progress(x)
        done(end+1) = x;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = throws(fn, text)
% true when fn stops with an error whose message contains text
try
    fn();
    ok = false;
catch err
    ok = contains(err.message, text);
    if ~ok
        fprintf('        (message: %s)\n', err.message);
    end
end
end

function ok = report(name, ok)
if ok
    fprintf("  PASS  %s\n", name);
else
    fprintf("  FAIL  %s\n", name);
end
end
