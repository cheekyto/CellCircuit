function [fig, result] = plotMEArecording(data, plot_type, inputs)

% VISUALIZATION OF A RECORDING (revised November 2026: active function, moved out of the folder "legacy" and driven by the tab
% "4. Visualization" of the window "CellCircuit"; "makeVideo_MEA" merged as the plot "video"; user's choice: the old options 3 -shuffling-,
% 4 -hypersynchronized units-, 7 -autocorrelation- and 8 -raw raster plot- removed, the old options 0 and 1 merged)
% Plots of the units of a dataset of the pipeline (file "<name>_QCmanual.mat": the units after the quality control; also the curated,
% integrated and cropped datasets), for the units and the time window chosen:
%   "raster"     raster plot (one row per unit, numbered as in the dataset), with the start and the name of each recording condition;
%                optional: the population firing rate below it (old options 0 and 1)
%   "bursts"     raster plot in which the bursts of each unit are colored, with marks at the first spike, last spike, peak and middle of
%                each burst ("BurstDetection"), to check the burst detection by eye (old option 2)
%   "map"        electrode map: the 64 x 64 electrodes, each colored by the firing rate of its units in the window, with the outline of
%                the organoid (old option 5, repaired: units after the quality control, automatic color scale, outline of the cropping)
%   "waveforms"  spike waveform of each unit: its Kilosort template on its peak electrode ("loadKilosortTemplates"; old option 6, repaired)
%   "video"      MP4 video of the electrode map, step by step (firing rate of each electrode in each step; "makeVideo_MEA"), with the
%                raster plot and/or the population firing rate of the window above it and a cursor that moves with the video
%   "data"       no plot: reads the dataset (units, conditions, outline of the organoid) once, to draw several plots without reading it again
% Electrode coordinates: x = positionXY(1)+1 (to the right) and y = positionXY(2)+1 (upwards), from 1 to 64, as in "SelectElectrodes" and
% "OrganoidsFromPhoto" (photo turned so that the "3D" text of the chip is on the right). Outline of the organoid (map and video): the one
% saved by the cropping ("SelectElectrodes": outline from the photo, or ellipse), otherwise the organoids found in the photo of the chip in
% the folder of the recording ("OrganoidsFromPhoto"), otherwise none
% (user's choices, 14 November 2026) no titles; 4-5 unit numbers on the raster plot; the names of the conditions above the plot; the bin of
% the population firing rate (200 ms by default) and the unit of the time axis (s by default) can be chosen; the population rate in blue;
% a quick preview in a panel (parent, preview) for the tab of the window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   data:           dataset file ("<name>_QCmanual.mat", with spikeTimestamp and rec_conditions), or the structure returned by "data"
%   plot_type:      "raster" (default), "bursts", "map", "waveforms", "video" or "data"
%   inputs:         (optional) structure that overwrites any of the choices below and of the PARAMETERS, e.g., struct('units', 1:20,
%                   'condition', "Baseline", 'population_rate', false)

units = [];                     % units to plot: their numbers in spikeTimestamp (as in the raster plot); [] (or -1) = all
window = [];                    % [start end] time window, in sec; [] = the window of "condition", or the whole recording
condition = 0;                  % recording condition whose window is plotted (when window = []): its number, or its name (e.g., "Baseline");
                                % =0, the whole recording
population_rate = true;         % (raster, video) population firing rate: below the raster plot, or above the electrode map of the video
raster = true;                  % (video) raster plot of the window above the electrode map
condition_lines = true;         % (raster, bursts, video) start and name of each recording condition
rate_measure = "mean";          % (map) firing rate of each unit in the window: "mean" = number of spikes / length of the window; "median" =
                                % 1 / median of the inter-spike intervals; "fast" = 1 / 10th percentile of the intervals (firing in bursts)
color_limit = [];               % (map, video) rate at the top of the color scale, in Hz; [] = automatic (AUTO_PERCENTILE)
empty_color = "black";          % (map, video) color of the electrodes without units: "black" (default) or "white" (user's request,
                                % November 2026); the outline of the organoid is drawn in the other one
outline = true;                 % (map, video) outline of the organoid(s)
photo = "";                     % (outline) photo of the chip (file in the folder of the recording, or full path); "" = the one found in the
                                % folder ("recording_picture"); only used when the dataset has no outline of its cropping
photo_text_side = "";           % (outline) side of the photo where the "3D" text is (as the column picture_3D_side); "" = found
photo_n_organoids = 0;          % (outline) number of organoids in the photo (as the column picture_n_organoids); 0 = found
output = "";                    % file to save: .png, .pdf or .fig for the plots, .mp4 for the video; "auto" = next to the dataset, with the
                                % name of the plot and of the window; "" = the plot is not saved (the video is always saved: "" = "auto")
visible = true;                 % =false, the figure is created without being shown (e.g., to save it)
progress = [];                  % (video) function called with the share of the frames done (0 to 1), e.g., to show the progress
condition_names = "treatment";  % (raster, bursts, video) names of the recording conditions: "treatment" = condition and dose, as the compiled
                                % tables group them (e.g., "+Gabazine 20uM"; "label_condition"); "file" = from the file names (e.g., "gabazine_20uM")
rate_bin = 0.2;                 % (raster, video) bin of the population firing rate, in sec (user's choice: 200 ms by default, as the old option 1)
time_unit = "s";                % (raster, bursts, video) unit of the time axis: "ms", "s" (default) or "min" (window is always in sec)
parent = [];                    % container to draw in (e.g., a panel of the window "CellCircuit"), instead of a new figure
preview = false;                % =true, quick preview (window "CellCircuit"): the raster plot is drawn as an image (one column of pixels per
                                % time bin, fast for any number of spikes), the plot is not saved, and the video is not written (its first
                                % frame is drawn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   fig:            the figure, or "parent" when it is given ([] for "data", and for "video": its figure is closed once the video is written)
%   result:         structure with what was drawn:
%       all plots       .units (units drawn), .window ([start end], in sec), .file (file saved, "" if none)
%       raster          .n_spikes (spikes drawn), .rate_time and .rate (population firing rate: center of each bin, in sec, and spikes of
%                       the units drawn in the bin / number of units / bin, in Hz)
%       bursts          .n_bursts (bursts of each unit), .bursts{k} ([first spike, last spike, peak, middle] of each burst of the unit k, in sec)
%       map             .map (MATRIX_SIZE x MATRIX_SIZE, map(y, x) = rate of the electrode (x, y), in Hz; NaN = no unit drawn there),
%                       .unit_rate (rate of each unit drawn), .color_limit, .outline_source (where the outline comes from)
%       waveforms       .templates (one row per unit), .sampling_rate (Hz), .unit_type (bombcell: 0 noise, 1 good, 2 MUA, 3 non-somatic)
%       video           .file ("" for a preview), .n_frames, .frame_rate (frames per sec), .frames (MATRIX_SIZE x MATRIX_SIZE x n_frames, as .map), .frame_time
%                       (center of each step, in sec), .color_limit, .outline_source
%       data            the dataset: .file, .folder, .name, .recording, .spikeTimestamp, .n_units, .position ([x y] electrode of each
%                       unit), .conditions (.name -as condition_names-, .file, .start, .stop -in sec-, .file_label and .treatment_label of
%                       each condition), .base_cond, .condition_names, .t_range, .outlines (cell
%                       array of [x y] polygons, in electrodes), .outline_names, .outline_source (where the outline comes from), .photo_searched
% The electrode of several units: the sum of their rates ("mean", and the video: firing of the electrode), or the largest ("median", "fast")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;            % factor to rescale from msec to sec
RATE_COLOR = [0 0.35 0.80];     % color of the population firing rate (user's choice: blue)
CONDITION_COLOR = [0.85 0 0];   % color of the start and the name of each recording condition
VIDEO_STEP = 0.4;               % (video) length of each step, in sec: each frame shows the firing rate of the electrodes in one step (as makeVideo_MEA)
VIDEO_SPEED = 4;                % (video) speed, times the real time (0.4-s steps at 10 frames per sec, as makeVideo_MEA)
FAST_QUANTILE = 0.1;            % (rate_measure "fast") quantile of the inter-spike intervals
AUTO_PERCENTILE = 99;           % automatic color limit: this percentile of the rates of the electrodes with units (video: of the rates above 0)
BORDER = 0;                     % (map, video) electrodes at each edge of the array that are not shown (makeVideo_MEA hid 3-4, ~4% of 64)
BURST_METHOD = 1;               % (bursts) burst detection method of "BurstDetection" (=1, Maximal Interval)
MAX_BURST_UNITS = 100;          % (bursts) most units in one plot
MAX_WAVEFORM_UNITS = 36;        % (waveforms) most units in one plot (one panel each)
MAP_COLORMAP = "jet";           % (map, video) colors of the rates (as makeVideo_MEA); (user's request, November 2026) any colormap of
                                % MATLAB can be chosen, e.g., "turbo", "parula", "hot", "gray" (the tab of the window offers a few)
MATRIX_SIZE = 64;               % dimension of the electrode array (64*64=4096)
VIDEO_WIDTH = 760;              % (video) width of the frames, in points (the map is about this size; each plot above it adds VIDEO_PANEL)
VIDEO_PANEL = 150;              % (video) height of the raster plot and of the population firing rate, in points

if nargin < 2 || strlength(string(plot_type)) == 0, plot_type = "raster"; end
if nargin < 3, inputs = struct(); end
plot_type = lower(string(plot_type));
if ~ismember(plot_type, ["raster", "bursts", "map", "waveforms", "video", "data"])
    error('plotMEArecording: the plot must be "raster", "bursts", "map", "waveforms", "video" or "data" (not "%s")', plot_type);
end
outline_asked = isfield(inputs, 'outline') && inputs.outline;               % ("data": the photo is searched only when asked)
override_inputs(inputs);                    % overwrite the choices / PARAMETERS above with the ones provided as input (if any)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dataset

if isstruct(data)
    D = data;
else
    D = read_dataset(string(data), K_MSEC_SEC);
end
D = name_conditions(D, condition_names);
fig = [];
if plot_type == "data"
    if outline_asked
        D = photo_outline(D, photo, photo_text_side, photo_n_organoids);
    end
    result = D;
    return
end
if outline && ismember(plot_type, ["map", "video"])
    D = photo_outline(D, photo, photo_text_side, photo_n_organoids);
end

units = unit_list(D, units);
[t0, t1] = time_window(D, window, condition);
T = spike_times(D, units, t0, t1, K_MSEC_SEC);                             % spike times of each unit in the window, in sec
nu = numel(units);
[tf, tu] = time_factor(time_unit);                                         % (the axes show the time in "time_unit": t * tf)
if ~(isnumeric(rate_bin) && isscalar(rate_bin) && rate_bin > 0)
    error('plotMEArecording: rate_bin (bin of the population firing rate) must be a time above 0, in sec (e.g., 0.2)');
end
result = struct('units', units, 'window', [t0 t1], 'file', "");
default_name = sprintf('%s_%s_%s-%ss', D.name, plot_type, num_text(t0), num_text(t1));

switch plot_type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case "raster"           % raster plot (+ population firing rate) (old options 0 and 1)

        fig = new_figure(visible, parent, "Raster plot - " + D.name, [80 80 1100 650]);
        if population_rate
            tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
            ax = nexttile(tl, [2 1]);
        else
            tl = tiledlayout(fig, 1, 1, 'Padding', 'compact');
            ax = nexttile(tl);
        end
        result.n_spikes = draw_raster(ax, T, units, [0 0 0], [1 1 1], t0, t1, tf, preview);
        xlim(ax, [t0 t1]*tf);
        if condition_lines
            draw_conditions(ax, D, t0, t1, tf, CONDITION_COLOR, tl);
        end
        [result.rate_time, result.rate, rate_edges] = population(T, t0, t1, rate_bin);
        if population_rate
            ax.XTickLabel = {};
            ax2 = nexttile(tl);
            stairs(ax2, rate_edges*tf, [result.rate, result.rate(end)], 'Color', RATE_COLOR, 'LineWidth', 0.9);   % (a line: still visible
                                                                            %  with thousands of bins, unlike bars)
            xlim(ax2, [t0 t1]*tf);
            if condition_lines
                draw_conditions(ax2, D, t0, t1, tf, CONDITION_COLOR);
            end
            ylabel(ax2, 'Population rate (Hz)');      xlabel(ax2, "Time (" + tu + ")");
            linkaxes([ax ax2], 'x');
        else
            xlabel(ax, "Time (" + tu + ")");
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case "bursts"           % raster plot + bursts of each unit (old option 2)

        if nu > MAX_BURST_UNITS
            error('plotMEArecording: the burst plot shows at most %d units (%d chosen): choose fewer units (e.g., 1:50)', MAX_BURST_UNITS, nu);
        end
        COLORS = [0 0.60 0.90; 0.90 0.50 0; 0 0.65 0.30; 0.60 0.30 0.80; 0.85 0.65 0];   % colors of the bursts (one after the other)
        fig = new_figure(visible, parent, "Bursts - " + D.name, [80 80 1100 650]);
        tl = tiledlayout(fig, 1, 1, 'Padding', 'compact');
        ax = nexttile(tl);
        hold(ax, 'on');     box(ax, 'on');
        result.n_bursts = zeros(nu, 1);     result.bursts = cell(nu, 1);
        t_spk = cell(nu, 1);    r_spk = cell(nu, 1);    c_spk = cell(nu, 1);    marks = cell(nu, 1);
        for k=1:nu
            s = T{k};
            lab = zeros(numel(s), 1);
            b = zeros(0, 4);
            if numel(s) >= 3
                B = BurstDetection(s(:)/K_MSEC_SEC, BURST_METHOD);          % (in msec, as in the pipeline)
                lab = B.spikes(:,3);
                if B.nBurst > 0
                    b = [[B.IndivBurst.start]' [B.IndivBurst.last]' [B.IndivBurst.peak]' [B.IndivBurst.median]'] * K_MSEC_SEC;
                end
            end
            result.n_bursts(k) = size(b, 1);     result.bursts{k} = b;
            t_spk{k} = s(:);     r_spk{k} = k*ones(numel(s), 1);     c_spk{k} = lab(:);
            marks{k} = [b, k*ones(size(b,1), 1)];
        end
        t_spk = vertcat(t_spk{:})*tf;      r_spk = vertcat(r_spk{:});      c_spk = vertcat(c_spk{:});      marks = vertcat(marks{:});
        ms = marker_size(nu);
        plot(ax, t_spk(c_spk == 0), r_spk(c_spk == 0), '|', 'Color', [0 0 0], 'MarkerSize', ms);     % spikes outside the bursts (black)
        for ci=1:size(COLORS, 1)
            in = c_spk > 0 & mod(c_spk - 1, size(COLORS, 1)) + 1 == ci;
            plot(ax, t_spk(in), r_spk(in), '|', 'Color', COLORS(ci,:), 'MarkerSize', ms, 'LineWidth', 1.5);
        end
        if isempty(marks)
            marks = nan(1, 5);                                              % (no burst: the marks are still listed in the legend)
        end
        marks(:, 1:4) = marks(:, 1:4)*tf;
        h = gobjects(1, 4);
        h(1) = plot(ax, marks(:,1), marks(:,5) + 0.3, '*', 'Color', [0.85 0 0], 'DisplayName', 'first spike');
        h(2) = plot(ax, marks(:,2), marks(:,5) + 0.3, '+', 'Color', [0 0 0.9], 'DisplayName', 'last spike');
        h(3) = plot(ax, marks(:,3), marks(:,5) + 0.3, 'o', 'Color', [0 0 0], 'DisplayName', 'peak (shortest intervals)');
        h(4) = plot(ax, marks(:,4), marks(:,5) + 0.3, 'o', 'Color', [0.85 0 0.85], 'DisplayName', 'middle spike');
        legend(ax, h, 'Location', 'southoutside', 'Orientation', 'horizontal', 'AutoUpdate', 'off', 'Box', 'off');   % (below: the names
        unit_ticks(ax, units);                                              %  of the conditions are above the plot)
        xlim(ax, [t0 t1]*tf);      xlabel(ax, "Time (" + tu + ")");
        if condition_lines
            draw_conditions(ax, D, t0, t1, tf, CONDITION_COLOR, tl);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case "map"              % electrode map, color-coded firing rate (old option 5)

        [M, result.unit_rate] = electrode_map(D, units, T, t0, t1, rate_measure, FAST_QUANTILE, MATRIX_SIZE, BORDER);
        lim = color_limit;
        if isempty(lim)
            lim = auto_limit(M(~isnan(M)), AUTO_PERCENTILE);
        end
        result.map = M;     result.color_limit = lim;     result.outline_source = D.outline_source;
        fig = new_figure(visible, parent, "Electrode map - " + D.name, [100 100 700 620]);
        tl = tiledlayout(fig, 1, 1, 'Padding', 'compact');
        ax = nexttile(tl);
        [~, ink] = draw_map(ax, M, lim, MAP_COLORMAP, MATRIX_SIZE, empty_color);
        cb = colorbar(ax);
        cb.Label.String = rate_label(rate_measure);
        if outline
            draw_outlines(ax, D, ink);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case "waveforms"        % spike waveforms: Kilosort templates (old option 6)

        if nu > MAX_WAVEFORM_UNITS
            error('plotMEArecording: the waveform plot shows at most %d units (%d chosen): choose fewer units', MAX_WAVEFORM_UNITS, nu);
        end
        try
            [W, unit_type, fs] = loadKilosortTemplates(D.folder, D.spikeTimestamp(units));
        catch ME
            error(['plotMEArecording: the spike waveforms are the Kilosort templates written next to the spike times (metadata.csv and ' ...
                'waveforms.mat in the folder of the recording; Kilosort units only). They could not be read: %s'], ME.message);
        end
        result.templates = W;   result.sampling_rate = fs;      result.unit_type = unit_type;
        TYPE_TEXT = ["noise", "good", "MUA", "non-somatic"];
        fig = new_figure(visible, parent, "Spike waveforms - " + D.name, [80 80 1000 700]);
        tl = tiledlayout(fig, 'flow', 'TileSpacing', 'compact', 'Padding', 'compact');
        for k=1:nu
            ax = nexttile(tl);
            [~, trough] = min(W(k,:));
            plot(ax, ((1:size(W,2)) - trough)/fs/K_MSEC_SEC, W(k,:), 'k', 'LineWidth', 1.2);
            type_k = "";
            if unit_type(k) >= 0 && unit_type(k) <= 3
                type_k = ", " + TYPE_TEXT(unit_type(k) + 1);
            end
            title(ax, sprintf('Unit %d (x %s, y %s%s)', units(k), num_text(D.position(units(k),1)), num_text(D.position(units(k),2)), type_k));
            box(ax, 'off');
        end
        xlabel(tl, 'Time from the trough (ms)');     ylabel(tl, 'Kilosort template (arbitrary units)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case "video"            % video of the electrode map (makeVideo_MEA), with the raster plot and/or the population rate

        n_frames = floor((t1 - t0)/VIDEO_STEP + 1e-9);
        if n_frames < 1
            error('plotMEArecording: the window (%g s) is shorter than one step of the video (%g s)', t1 - t0, VIDEO_STEP);
        end
        edges = t0 + (0:n_frames)*VIDEO_STEP;
        frame_time = edges(1:end-1) + VIDEO_STEP/2;
        F = nan(MATRIX_SIZE, MATRIX_SIZE, n_frames, 'single');              % F(y, x, k): rate of the electrode (x, y) in the step k, in Hz
        for k=1:nu
            p = D.position(units(k), :);
            if any(isnan(p)) || any(p < 1) || any(p > MATRIX_SIZE)
                continue
            end
            c = single(histcounts(T{k}, edges))/VIDEO_STEP;
            old = squeeze(F(p(2), p(1), :))';
            old(isnan(old)) = 0;
            F(p(2), p(1), :) = reshape(old + c, 1, 1, []);
        end
        F = hide_border(F, BORDER);
        lim = color_limit;
        if isempty(lim)
            lim = auto_limit(F(F > 0), AUTO_PERCENTILE);
        end
        file = string(output);
        if strlength(file) == 0 || file == "auto"
            file = string(fullfile(D.folder, default_name + ".mp4"));
        end
        profile = 'MPEG-4';
        if ~ispc && ~ismac                                                  % (MATLAB writes MPEG-4 only on Windows and macOS)
            profile = 'Motion JPEG AVI';
            file = regexprep(file, '\.mp4$', '.avi');
        end
        n_panels = double(raster) + double(population_rate);
        dark = [0.92 0.92 0.92];
        fig = new_figure(visible, parent, "Video - " + D.name, [60 60 VIDEO_WIDTH VIDEO_WIDTH*0.78 + n_panels*VIDEO_PANEL]);
        set_background(fig, [0 0 0]);
        tl = tiledlayout(fig, n_panels + 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        cursors = gobjects(0);
        ax_time = [];                                                       % (the lowest plot above the map carries the time axis)
        if raster
            axr = nexttile(tl);     ax_time = axr;
            draw_raster(axr, T, units, dark, [0 0 0], t0, t1, tf, preview);
            xlim(axr, [t0 t1]*tf);     dark_axes(axr, dark);
            if condition_lines
                draw_conditions(axr, D, t0, t1, tf, [1 0.35 0.35]);         % (dashed lines only: no names in the video)
            end
            cursors(end+1) = xline(axr, frame_time(1)*tf, 'r', 'LineWidth', 1.5);
            if population_rate, axr.XTickLabel = {}; end
        end
        if population_rate
            axp = nexttile(tl);     ax_time = axp;
            [rt, rr] = population(T, t0, t1, rate_bin);
            plot(axp, rt*tf, rr, 'Color', [0.45 0.75 1], 'LineWidth', 1);
            xlim(axp, [t0 t1]*tf);     dark_axes(axp, dark);     ylabel(axp, 'Rate (Hz)');
            if condition_lines && ~raster
                draw_conditions(axp, D, t0, t1, tf, [1 0.35 0.35]);
            end
            cursors(end+1) = xline(axp, frame_time(1)*tf, 'r', 'LineWidth', 1.5);
        end
        if ~isempty(ax_time)
            xlabel(ax_time, "Time (" + tu + ")", 'Color', dark);
        end
        axm = nexttile(tl, [3 1]);
        [im, ink] = draw_map(axm, F(:,:,1), lim, MAP_COLORMAP, MATRIX_SIZE, empty_color);
        dark_axes(axm, dark);
        axm.Color = 1 - ink;                                                % (the electrodes without units: empty_color)
        cb = colorbar(axm, 'Color', dark);
        cb.Label.String = 'Firing rate (Hz)';
        if outline
            draw_outlines(axm, D, ink);
        end
        frame_title = @(k) title(tl, sprintf('Time: %s %s', time_text(frame_time(k)*tf, tu), tu), ...   % (user's choice: the time only)
            'Color', [1 0.35 0.35], 'FontWeight', 'bold');
        result.n_frames = n_frames;     result.frame_rate = VIDEO_SPEED/VIDEO_STEP;     result.frames = F;
        result.frame_time = frame_time;     result.color_limit = lim;     result.outline_source = D.outline_source;
        if preview                                                          % (preview: the first frame only, no file)
            frame_title(1);
            return
        end
        v = VideoWriter(char(file), profile);
        v.FrameRate = VIDEO_SPEED/VIDEO_STEP;
        if strcmp(profile, 'MPEG-4'), v.Quality = 95; end
        open(v);
        closer = onCleanup(@() close(v));
        frame_size = [];
        for k=1:n_frames
            im.CData = F(:,:,k);
            im.AlphaData = ~isnan(F(:,:,k));
            for j=1:numel(cursors)
                cursors(j).Value = frame_time(k)*tf;
            end
            frame_title(k);
            drawnow;
            img = getframe(fig).cdata;
            if isempty(frame_size)
                frame_size = 2*floor([size(img,1) size(img,2)]/2);          % (even sizes, for the MPEG-4 encoder)
            end
            img = fit_frame(img, frame_size);
            writeVideo(v, img);
            if ~isempty(progress) && (k == n_frames || mod(k, max(1, round(n_frames/20))) == 0)
                progress(k/n_frames);
            end
        end
        clear closer
        if isempty(parent)
            close(fig);
            fig = [];
        end
        result.file = file;
        return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the plot

file = string(output);
if strlength(file) > 0 && ~preview
    if file == "auto"
        file = string(fullfile(D.folder, default_name + ".png"));
    end
    if endsWith(lower(file), ".fig")
        fig.CreateFcn = 'set(gcbo, ''Visible'', ''on'')';                   % (a hidden figure is saved so that its file opens visible)
        savefig(fig, char(file));
    else
        exportgraphics(fig, char(file), 'Resolution', 200);
    end
    result.file = file;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   LOCAL FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = read_dataset(file, K_MSEC_SEC)

% reads the units, the recording conditions and the outline of the cropping (if any) of a dataset file

if ~isfile(file)
    error('plotMEArecording: the file %s does not exist', file);
end
vars = string({whos('-file', file).name});
if ~ismember("spikeTimestamp", vars)
    error('plotMEArecording: %s has no variable spikeTimestamp (choose a "_QCmanual.mat" file of the pipeline)', file);
end
want = intersect(["spikeTimestamp", "rec_conditions", "crop_outline", "ellipse_cropped"], vars, 'stable');
S = load(file, want{:});
[folder, base] = fileparts(file);
name = regexprep(string(base), '_QCmanual$', '');
st = S.spikeTimestamp;
n = numel(st);

position = nan(n, 2);                                                       % electrode of each unit: x = positionXY(1)+1, y = positionXY(2)+1
t_last = 0;
for i=1:n
    if isfield(st, 'positionXY') && numel(st(i).positionXY) >= 2
        p = double(st(i).positionXY(:))';
        position(i,:) = p(1:2) + 1;
    end
    if ~isempty(st(i).timeStamp)
        t_last = max(t_last, double(max(st(i).timeStamp))*K_MSEC_SEC);
    end
end

conditions = struct('name', {}, 'file', {}, 'start', {}, 'stop', {}, 'file_label', {}, 'treatment_label', {});
base_cond = 1;
if isfield(S, 'rec_conditions') && ~isempty(S.rec_conditions)
    rc = S.rec_conditions;
    nc = size(rc, 1);
    stems = strings(nc, 1);
    for c=1:nc
        aux = split(string(rc{c,1}), ".b");                                 % (file name without the extension .brw / .bxr)
        stems(c) = aux(1);
    end
    base_cond = baselineCondition(rc(:,1));
    for c=1:nc
        w = double(rc{c,2});
        if c == base_cond
            label = "Baseline";
        elseif startsWith(stems(c), stems(base_cond))
            label = regexprep(extractAfter(stems(c), strlength(stems(base_cond))), '^[_\-\s]+', '');
        else
            label = stems(c);
        end
        if strlength(label) == 0
            label = "condition " + c;
        end
        conditions(c) = struct('name', label, 'file', string(rc{c,1}), 'start', w(1), 'stop', w(2), 'file_label', label, ...
            'treatment_label', label); %#ok<AGROW>
    end
    others = setdiff(1:nc, base_cond);                                      % (user's request) the part of the file name that every other
    if numel(others) >= 2                                                   % condition shares (e.g., "100div_" and "_00") is left out
        labels = shorter_names([conditions(others).file_label]);
        for c=1:numel(others)
            conditions(others(c)).file_label = labels(c);
        end
    end
    [treatment, ratio] = treatment_settings(string(folder), regexprep(name, '_(curated|integrated)$', ''));
    for c=1:nc                                                              % (user's request) the condition and dose by which
        try                                                                 % "Compile_MEAproject" groups the recordings ("label_condition",
            [label, dose] = label_condition(rc, c, base_cond, treatment, ratio);   % as "Analyze_spikeTrain")
            label = string(label);
            if strlength(string(dose)) > 0
                label = label + " " + string(dose);
            end
            conditions(c).treatment_label = label;
        catch                                                               % (file names that do not follow the rule: the file label)
            conditions(c).treatment_label = conditions(c).file_label;
        end
    end
else
    conditions(1) = struct('name', "Recording", 'file', name, 'start', 0, 'stop', t_last, 'file_label', "Recording", 'treatment_label', "Recording");
end

D = struct('file', file, 'folder', string(folder), 'name', name, 'recording', regexprep(name, '_(curated|integrated)$', ''), ...
    'spikeTimestamp', st, 'n_units', n, 'position', position, 'conditions', conditions, 'base_cond', base_cond, ...
    't_range', [min([conditions.start]) max([conditions.stop])], 'outlines', {{}}, 'outline_names', strings(0,1), ...
    'outline_source', "none", 'photo_searched', false, 'condition_names', "");
D = name_conditions(D, "treatment");

D = cropping_outline(D, S);
if isempty(D.outlines) && D.recording ~= name                               % analyzed dataset ("_curated", "_integrated"): the outline of the
    source = fullfile(folder, D.recording + "_QCmanual.mat");               % dataset it was made from (a cropped one), not copied by
    if isfile(source)                                                       % "IntegrateDataSets"
        crop_vars = intersect(["crop_outline", "ellipse_cropped"], string({whos('-file', source).name}), 'stable');
        if ~isempty(crop_vars)
            D = cropping_outline(D, load(source, crop_vars{:}));
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = name_conditions(D, mode)

% names of the conditions shown: "treatment" = condition and dose, as grouped by "Compile_MEAproject" (e.g., "+Gabazine 20uM"); "file" =
% the part of the file name after the name of the baseline, without what all the conditions share (e.g., "gabazine_20uM"). Two conditions
% with the same name get their numbers (e.g., "+TTX (2)" and "+TTX (3)")

mode = lower(string(mode));
switch mode
    case "treatment"
        names = [D.conditions.treatment_label];
    case "file"
        names = [D.conditions.file_label];
    otherwise
        error('plotMEArecording: condition_names must be "treatment" (condition and dose) or "file" (file name), not "%s"', mode);
end
for c=1:numel(names)
    if sum(names == names(c)) > 1
        D.conditions(c).name = names(c) + " (" + c + ")";
    else
        D.conditions(c).name = names(c);
    end
end
D.condition_names = mode;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [treatment, ratio] = treatment_settings(folder, recording)

% Analyze_spikeTrain.TREATMENT and RATIO_MATCH of the parameter set used by the last analysis of the recording ("<file_name>_parameters.mat"
% written by "run_pipeline", in its folder; the name of a cropped or analyzed dataset is shortened until that file is found), otherwise
% of pipeline_parameters.m

P = [];
parts = split(recording, "_");
for k=numel(parts):-1:1
    file = fullfile(folder, strjoin(parts(1:k), "_") + "_parameters.mat");
    if isfile(file)
        S = load(file);
        if isfield(S, 'parameters_stage') && isfield(S.parameters_stage, 'Analyze_spikeTrain')
            P = S.parameters_stage.Analyze_spikeTrain;
        elseif isfield(S, 'parameters')
            P = S.parameters;
        end
        break
    end
end
if ~(isstruct(P) && isfield(P, 'Analyze_spikeTrain') && all(isfield(P.Analyze_spikeTrain, {'TREATMENT', 'RATIO_MATCH'})))
    P = load_parameters(fullfile(fileparts(mfilename('fullpath')), 'pipeline_parameters.m'));
end
treatment = P.Analyze_spikeTrain.TREATMENT;
ratio = P.Analyze_spikeTrain.RATIO_MATCH;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function names = shorter_names(names)

% names of the conditions without the start and the end that all of them share, cut at a "_" (e.g., "100div_TTX_2uM_00" and
% "100div_CGP_20uM_00" -> "TTX_2uM" and "CGP_20uM"); unchanged when they are all the same or a name would be empty

names = string(names);
if numel(unique(names)) < 2
    return
end
c = char(names(1));
for k=2:numel(names)                                                        % start shared by all, up to its last "_"
    n = char(names(k));     m = min(numel(c), numel(n));
    c = c(1:find([c(1:m) ~= n(1:m), true], 1) - 1);
end
cut = find(c == '_', 1, 'last');
if ~isempty(cut)
    trimmed = extractAfter(names, cut);
    if all(strlength(trimmed) > 0), names = trimmed; end
end
r = fliplr(char(names(1)));
for k=2:numel(names)                                                        % end shared by all, from its first "_"
    n = fliplr(char(names(k)));     m = min(numel(r), numel(n));
    r = r(1:find([r(1:m) ~= n(1:m), true], 1) - 1);
end
cut = find(r == '_', 1, 'last');                                            % (r is reversed: its last "_" is the first one of the end)
if ~isempty(cut)
    trimmed = extractBefore(names, strlength(names) - cut + 1);
    if all(strlength(trimmed) > 0), names = trimmed; end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = cropping_outline(D, S)

% outline saved by the cropping ("SelectElectrodes"), if any

if isfield(S, 'crop_outline') && ~isempty(S.crop_outline)                   % option 3: outline of the organoid in the photo
    D.outlines = {double(S.crop_outline)};
    D.outline_names = D.recording;
    D.outline_source = "cropping (photo)";
elseif isfield(S, 'ellipse_cropped') && isnumeric(S.ellipse_cropped) && numel(S.ellipse_cropped) == 4 && ~any(isnan(S.ellipse_cropped))
    e = double(S.ellipse_cropped);                                          % option 2: ellipse [x0 y0 radiusx radiusy]
    a = linspace(0, 2*pi, 121)';
    D.outlines = {[e(1) + e(3)*cos(a), e(2) + e(4)*sin(a)]};
    D.outline_names = D.recording;
    D.outline_source = "cropping (ellipse)";
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = photo_outline(D, photo, text_side, n_organoids)

% outline of the organoids found in the photo of the chip ("OrganoidsFromPhoto"), when the dataset has no outline of its cropping (once)

if ~isempty(D.outlines) || D.photo_searched
    return
end
D.photo_searched = true;
photo = string(photo);
if strlength(photo) == 0
    photo = recording_picture(D.folder, D.recording);
    if strlength(photo) == 0
        D.outline_source = "none (no photo of the chip in the folder)";
        return
    end
end
if isempty(regexp(char(photo), '^(/|~|[A-Za-z]:[\\/]|\\\\)', 'once'))
    photo = string(fullfile(D.folder, photo));                              % (a file in the folder of the recording)
end
[~, pname, pext] = fileparts(photo);
try
    res = OrganoidsFromPhoto(photo, struct('FIGURE', false, 'TEXT_SIDE', lower(string(text_side)), 'N_ORGANOIDS', double(n_organoids)));
catch ME
    warning('plotMEArecording:photo', 'No outline: the organoids could not be found in %s (%s)', photo, ME.message);
    D.outline_source = "none (photo not read)";
    return
end
if res.n_organoids == 0
    D.outline_source = sprintf("none (no organoid found in %s%s)", pname, pext);
    return
end
D.outlines = {res.organoids.outline};
D.outline_names = "organoid " + string(1:res.n_organoids)';
D.outline_source = sprintf("photo %s%s, %d organoid(s)", pname, pext, res.n_organoids);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function units = unit_list(D, units)

% units chosen (numbers in spikeTimestamp), checked; [] or -1 = all

if isempty(units) || isequal(units, -1)
    units = 1:D.n_units;
    return
end
units = double(units(:))';
if any(units ~= round(units)) || any(units < 1) || any(units > D.n_units)
    bad = units(units ~= round(units) | units < 1 | units > D.n_units);
    error('plotMEArecording: the dataset has units 1 to %d (unit %s does not exist)', D.n_units, num_text(bad(1)));
end
[~, first] = unique(units, 'stable');
units = units(first);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t0, t1, label] = time_window(D, window, condition)

% time window: given, the window of a condition, or the whole recording

if ~isempty(window)
    if ~isnumeric(window) || numel(window) ~= 2 || window(2) <= window(1)
        error('plotMEArecording: the window must be [start end], in sec, with end > start');
    end
    t0 = double(window(1));     t1 = double(window(2));
    label = sprintf('%s-%s s', num_text(t0), num_text(t1));
elseif ~(isnumeric(condition) && isscalar(condition) && condition == 0) && ~(isstring(condition) && strlength(condition) == 0)
    names = [D.conditions.name];
    if isnumeric(condition)
        c = condition;
        if c < 1 || c > numel(D.conditions) || c ~= round(c)
            error('plotMEArecording: the recording has conditions 1 to %d', numel(D.conditions));
        end
    else
        c = find(strcmpi(names, string(condition)), 1);
        if isempty(c)
            error('plotMEArecording: no condition "%s" (conditions: %s)', string(condition), strjoin(names, ', '));
        end
    end
    t0 = D.conditions(c).start;     t1 = D.conditions(c).stop;
    label = sprintf('%s (%s-%s s)', names(c), num_text(t0), num_text(t1));
else
    t0 = D.t_range(1);      t1 = D.t_range(2);
    label = sprintf('whole recording (%s-%s s)', num_text(t0), num_text(t1));
end
if t1 <= D.t_range(1) || t0 >= D.t_range(2)
    error('plotMEArecording: the window %s-%s s is outside the recording (%s-%s s)', num_text(t0), num_text(t1), ...
        num_text(D.t_range(1)), num_text(D.t_range(2)));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = spike_times(D, units, t0, t1, K_MSEC_SEC)

% spike times of each unit in [t0, t1), in sec

T = cell(numel(units), 1);
for k=1:numel(units)
    s = double(D.spikeTimestamp(units(k)).timeStamp(:))*K_MSEC_SEC;
    T{k} = s(s >= t0 & s < t1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n_spikes = draw_raster(ax, T, units, color, background, t0, t1, tf, as_image)

% raster plot, one row per unit (row k = units(k)); the time in the unit of the axis (t * tf). As ticks in a single plot call (as
% "Analyze_spikeTrain"), or, for a preview, as an image (one column of pixels per time bin: fast for any number of spikes)

nu = numel(T);
counts = cellfun(@numel, T);
n_spikes = sum(counts);
if as_image
    n_px = 900;                                                             % (columns of the image)
    edges = linspace(t0, t1, n_px + 1);
    M = false(max(nu, 1), n_px);
    for k=1:nu
        M(k, :) = histcounts(T{k}, edges) > 0;
    end
    img = repmat(reshape(background, 1, 1, 3), size(M, 1), n_px);
    for ch=1:3
        layer = img(:, :, ch);      layer(M) = color(ch);      img(:, :, ch) = layer;
    end
    centers = (edges(1:end-1) + edges(2:end))/2;
    image(ax, centers([1 end])*tf, [1 max(nu, 1)], img);
    ax.YDir = 'normal';
else
    rows = repelem((1:nu)', counts(:));
    t = vertcat(T{:})*tf;
    if isempty(t)
        t = zeros(0, 1);    rows = zeros(0, 1);
    end
    plot(ax, t, rows, '|', 'Color', color, 'MarkerSize', marker_size(nu));
end
unit_ticks(ax, units);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function unit_ticks(ax, units)

% (user's choice) 4-5 unit numbers on the rows: round numbers (1, 10, 20, 30...) for consecutive units, otherwise evenly spaced rows

nu = numel(units);
ylim(ax, [0 nu + 1]);
if nu <= 5
    rows = 1:nu;
elseif all(diff(units) == 1)
    nice = [1 2 5]' * 10.^(0:6);   nice = nice(:)';
    [~, i] = min(abs(log(nice) - log(nu/4)));
    step = nice(i);
    tick_units = ceil(units(1)/step)*step : step : units(end);
    if numel(tick_units) < 5 && (isempty(tick_units) || tick_units(1) - units(1) >= step/2)
        tick_units = [units(1), tick_units];                                % (the first unit, when the first round number is far from it)
    end
    rows = tick_units - units(1) + 1;
else
    rows = unique(round(linspace(1, nu, 5)));
end
ax.YTick = rows;
ax.YTickLabel = string(units(rows));
ylabel(ax, 'Neuron');                                                       % (user's choice: "Neuron")

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ms = marker_size(nu)
ms = max(2, min(10, 400/max(nu, 1)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_conditions(ax, D, t0, t1, tf, color, tl)

% dashed line at the start of each recording condition in the window; with "tl" (the layout of the plot), the name of each condition above
% the plot (user's choice: clearly visible, out of the spikes), at its start, on a second row when it would overlap the name before it

names_above = nargin >= 7;
x = [];     names = strings(0);
for c=1:numel(D.conditions)
    s = D.conditions(c).start;      e = D.conditions(c).stop;
    if e <= t0 || s >= t1
        continue
    end
    if s > t0                                                               % (no line on the left edge of the plot)
        xline(ax, s*tf, '--', 'Color', color, 'LineWidth', 1.2, 'Tag', 'condition');
    end
    x(end+1) = max(s, t0)*tf; %#ok<AGROW>
    names(end+1) = D.conditions(c).name; %#ok<AGROW>
end
if ~names_above || isempty(names)
    return
end
y = ax.YLim(2);
right = -inf(1, 3);                                                         % right end of the last name of each row
n_rows = 1;
for c=1:numel(names)
    h = text(ax, x(c), y, names(c), 'Color', color, 'FontWeight', 'bold', 'FontSize', 11, 'Interpreter', 'none', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Clipping', 'off', 'Tag', 'condition_name');
    row = find(x(c) > right, 1);                                            % (the first row where it does not overlap)
    if isempty(row), row = 3; end
    if row > 1
        h.String = [{char(names(c))}; repmat({' '}, row - 1, 1)];          % (blank lines below the name lift it by one row each)
    end
    ext = h.Extent;
    right(row) = max(right(row), ext(1) + ext(3));
    n_rows = max(n_rows, row);
end
title(tl, repmat({' '}, n_rows, 1), 'FontSize', 11);                        % (room for the names above the plot)

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [centers, rate, edges] = population(T, t0, t1, bin)

% population firing rate (user's choice, November 2026): the non-normalized population signal of the network-burst detection, i.e., the
% spikes of all the units drawn in each bin / number of units / length of the bin (Hz per unit, labelled "Hz"; as the old option 1). Not
% the normalized one (proportion of the spikes of each unit in each bin, averaged across units, without unit)

n_bins = floor((t1 - t0)/bin + 1e-9);
if n_bins < 1
    centers = (t0 + t1)/2;      bin = t1 - t0;      edges = [t0 t1];
else
    edges = t0 + (0:n_bins)*bin;
    centers = edges(1:end-1) + bin/2;
end
counts = histcounts(vertcat(T{:}), edges);
rate = counts / max(numel(T), 1) / bin;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M, unit_rate] = electrode_map(D, units, T, t0, t1, measure, quant, MATRIX_SIZE, border)

% rate of each electrode: M(y, x), NaN where no unit is drawn

measure = lower(string(measure));
if ~ismember(measure, ["mean", "median", "fast"])
    error('plotMEArecording: rate_measure must be "mean", "median" or "fast" (not "%s")', measure);
end
M = nan(MATRIX_SIZE, MATRIX_SIZE);
unit_rate = zeros(numel(units), 1);
for k=1:numel(units)
    s = T{k};
    switch measure
        case "mean"
            r = numel(s)/(t1 - t0);
        case "median"
            r = 0;
            if numel(s) >= 3, r = 1/median(diff(s)); end
        case "fast"
            r = 0;
            if numel(s) >= 3, r = 1/quantile(diff(s), quant); end
    end
    unit_rate(k) = r;
    p = D.position(units(k), :);
    if any(isnan(p)) || any(p < 1) || any(p > MATRIX_SIZE)
        continue
    end
    old = M(p(2), p(1));
    if isnan(old)
        M(p(2), p(1)) = r;
    elseif measure == "mean"
        M(p(2), p(1)) = old + r;                                            % (firing of the electrode: sum of its units)
    else
        M(p(2), p(1)) = max(old, r);
    end
end
M = hide_border(M, border);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = hide_border(M, border)
if border > 0
    n = size(M, 1);
    edge = [1:border, n-border+1:n];
    M(edge, :, :) = NaN;
    M(:, edge, :) = NaN;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lim = auto_limit(values, pct)
lim = 1;
if ~isempty(values)
    lim = double(prctile(double(values(:)), pct));
end
if ~(lim > 0)
    lim = 1;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [im, ink] = draw_map(ax, M, lim, cmap, MATRIX_SIZE, empty_color)

% electrode map: x to the right, y upwards; the electrodes without units in empty_color ("black" or "white"); ink = the other color
% (for the outline of the organoid)

switch lower(string(empty_color))
    case "black"
        background = [0 0 0];
    case "white"
        background = [1 1 1];
    otherwise
        error('plotMEArecording: empty_color (electrodes without units) must be "black" or "white", not "%s"', empty_color);
end
ink = 1 - background;

im = imagesc(ax, [1 MATRIX_SIZE], [1 MATRIX_SIZE], M);
im.AlphaData = ~isnan(M);
hold(ax, 'on');
ax.Color = background;
ax.YDir = 'normal';
axis(ax, 'image');
xlim(ax, [0.5 MATRIX_SIZE + 0.5]);      ylim(ax, [0.5 MATRIX_SIZE + 0.5]);
clim(ax, [0 lim]);
try
    colors = feval(char(cmap), 256);
catch
    colors = [];
end
if ~(isnumeric(colors) && size(colors, 2) == 3)
    error('plotMEArecording: MAP_COLORMAP must be the name of a colormap of MATLAB (e.g., "jet", "turbo", "parula", "hot"), not "%s"', cmap);
end
colormap(ax, colors);
ax.XTick = [1 16 32 48 64];     ax.YTick = [1 16 32 48 64];
xlabel(ax, 'Electrode x');      ylabel(ax, 'Electrode y');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_outlines(ax, D, color)

% outline of each organoid (dashed), with its name when there are several

for k=1:numel(D.outlines)
    xy = D.outlines{k};
    plot(ax, xy([1:end 1], 1), xy([1:end 1], 2), '--', 'Color', color, 'LineWidth', 1.5, 'Tag', 'outline');
    if numel(D.outlines) > 1
        text(ax, mean(xy(:,1)), max(xy(:,2)) + 1.5, D.outline_names(k), 'Color', color, 'HorizontalAlignment', 'center', ...
            'Clipping', 'on');
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function txt = rate_label(measure)
switch lower(string(measure))
    case "median"
        txt = '1 / median interval between spikes (Hz)';
    case "fast"
        txt = '1 / 10th percentile of the intervals (Hz)';
    otherwise
        txt = 'Firing rate (Hz)';
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dark_axes(ax, color)
ax.Color = [0 0 0];     ax.XColor = color;      ax.YColor = color;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function img = fit_frame(img, sz)
% every frame of the video with the size of the first one (cropped, or padded with black)
out = zeros(sz(1), sz(2), 3, 'like', img);
h = min(sz(1), size(img,1));    w = min(sz(2), size(img,2));
out(1:h, 1:w, :) = img(1:h, 1:w, :);
img = out;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fig = new_figure(visible, parent, name, position)
% a new figure, or the container given (emptied first)
if ~isempty(parent)
    delete(parent.Children);
    fig = parent;
    set_background(fig, [1 1 1]);
    return
end
vis = 'on';
if ~visible, vis = 'off'; end
fig = figure('Visible', vis, 'Name', char(name), 'NumberTitle', 'off', 'Color', [1 1 1], 'Position', position);
end

function set_background(container, color)
% background of a figure (Color) or of a panel (BackgroundColor)
if isprop(container, 'BackgroundColor')
    container.BackgroundColor = color;
else
    container.Color = color;
end
end

function [tf, tu] = time_factor(time_unit)
% factor from sec to the unit of the time axis, and its name
tu = lower(string(time_unit));
switch tu
    case "ms"
        tf = 1000;
    case "s"
        tf = 1;
    case "min"
        tf = 1/60;
    otherwise
        error('plotMEArecording: time_unit must be "ms", "s" or "min" (not "%s")', time_unit);
end
end

function txt = time_text(t, tu)
% time of a frame of the video, with the digits that its unit needs
switch tu
    case "ms"
        txt = sprintf('%.0f', t);
    case "min"
        txt = sprintf('%.2f', t);
    otherwise
        txt = sprintf('%.1f', t);
end
end

function txt = num_text(x)
txt = sprintf('%g', x);
end
