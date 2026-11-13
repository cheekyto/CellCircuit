function passed = test_OrganoidsFromPhoto(code_dir)

% Function-level test of "OrganoidsFromPhoto" and of the cropping from the photo in "SelectElectrodes" (option 3) on synthetic photos of
% a chip (November 2026): a dark board, the big circle of the well, the package around the chip (dark square with its rows of leads and
% its rounded opening, as in the photos), the light chip surface with the strip on the left of the electrode array and the bands above and below it, and pale round organoids at known electrodes. It checks that
%   1. the array is found (corners within 1.5 electrodes), and each organoid (number, center within 1.5 electrodes, radius within 10%)
%   2. the electrodes kept: all the electrodes under the organoid, the ones within MARGIN x radius of its outline, none farther; two
%      touching organoids are split, no electrode within SHARED_GAP of the other organoid is kept, and no electrode is kept twice
%   3. the same photo turned (TEXT_SIDE "top", "bottom", "left") gives the same organoids; the corners given by hand (ARRAY_CORNERS) too
%   4. an elongated organoid is not split; N_ORGANOIDS forces the number; a small pale spot is not an organoid; a photo without organoid
%      stops with an error; three touching organoids are found with N_ORGANOIDS = 3; a thin pale strip along the edge of the array,
%      touching the organoid, is left out (MIN_WIDTH); the chip is found when the largest pale spot is not on it, and in a distant
%      photo (array of 50 pixels)
%   5. side of the "3D" text: found from the pattern of the chip, with the landmarks of the board (connector, triangle, text) in the 4
%      turns of a photo; a wrong side given is reported (warning) and the run goes on with it; without landmarks, the side is found
%      from the pattern of the chip (no warning), in the 4 turns
%   6. "SelectElectrodes" (option 3, its default) keeps the clusters on the electrodes kept and saves the 2 figures (electrode map, copy
%      of the photo); the recordings table built from the folder lists the recording, not its cropped dataset, and the cropping run
%      again from its row by "run_pipeline" with another CROPPING_SAFETY writes the cropped dataset again
%
% Example (from a terminal):
%   matlab -batch "addpath('<repo>/tests'); test_OrganoidsFromPhoto('<repo>')"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dir:             folder with the code to be tested

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% passed:               =true if all the checks passed (each check is printed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(char(code_dir));
passed = true;
N = 64;                                 % electrodes per side
quiet = struct('FIGURE', false, 'TEXT_SIDE', "right");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1-2. One organoid, and two touching organoids

[I1, chip] = synthetic_chip([33 30 22]);                                    % one organoid: center (33, 30), radius 22 electrodes
r = OrganoidsFromPhoto(I1, quiet);
passed = report("array found (corners within 1.5 electrodes)", corner_error(r, chip) < 1.5) && passed;
passed = report("one organoid: center within 1.5 electrodes, radius within 10%", r.n_organoids == 1 ...
    && norm(r.organoids(1).centroid - [33 30]) < 1.5 && abs(r.organoids(1).radius/22 - 1) < 0.10) && passed;
passed = report("equivalent ellipse close to the organoid", all(abs(r.organoids(1).ellipse - [33 30 22 22]) < [1.5 1.5 2.2 2.2])) && passed;
[EX, EY] = ndgrid(1:N, 1:N);
d = hypot(EX - 33, EY - 30);
under = d <= 22 - 1.5;                                                      % surely under the organoid
far = d > 22*1.10 + 1.5;                                                    % surely beyond the margin (MARGIN = 0.10)
crop = r.organoids(1).crop;
passed = report("electrodes kept: all those under the organoid, none beyond the margin", all(crop(under)) && ~any(crop(far))) && passed;
passed = report("electrodes kept: the list matches the map", isequal(sortrows(r.organoids(1).electrodes), sortrows([EX(crop) EY(crop)]))) && passed;

[I2, chip2] = synthetic_chip([24 22 15; 42 44 15]);                         % two touching organoids (distance 28.4, radii 15)
r2 = OrganoidsFromPhoto(I2, quiet);
ok = r2.n_organoids == 2 && corner_error(r2, chip2) < 1.5;
if ok                                                                       % numbered from the top: (42, 44) first
    ok = norm(r2.organoids(1).centroid - [42 44]) < 2 && norm(r2.organoids(2).centroid - [24 22]) < 2;
end
passed = report("two touching organoids are split (numbered from the top)", ok) && passed;
if r2.n_organoids == 2
    both = r2.organoids(1).crop & r2.organoids(2).crop;
    d1 = hypot(EX - 42, EY - 44) - 15;  d2 = hypot(EX - 24, EY - 22) - 15;  % distance to the edge of each organoid
    gap_ok = ~any(r2.organoids(1).crop(d2 <= 2 - 1)) && ~any(r2.organoids(2).crop(d1 <= 2 - 1));    % SHARED_GAP = 2 (1 electrode of tolerance)
    [a1, b1] = find(r2.organoids(1).crop);  [a2, b2] = find(r2.organoids(2).crop);
    D = hypot(a1 - a2', b1 - b2');                                          % distances between the electrodes kept for the two organoids
    contact = [24 22] + (15 + (hypot(18, 22) - 30)/2)*[18 22]/hypot(18, 22);  % point where the two organoids touch
    near1 = hypot(a1 - contact(1), b1 - contact(2)) <= 6;   near2 = hypot(a2 - contact(1), b2 - contact(2)) <= 6;
    passed = report("shared border: no electrode within SHARED_GAP of the other organoid, none kept twice", gap_ok && ~any(both(:))) && passed;
    passed = report(sprintf("shared border: the electrodes kept are %.1f electrodes apart across the border (>= 2 x SHARED_GAP), %.1f elsewhere (>= SHARED_GAP)", ...
        min(D(near1, near2), [], 'all'), min(D(:))), min(D(near1, near2), [], 'all') >= 2*2 - 0.5 && min(D(:)) >= 2) && passed;
    passed = report("shared border: the electrodes far from it are kept", all(r2.organoids(1).crop(hypot(EX - 42, EY - 44) <= 10)) ...
        && all(r2.organoids(2).crop(hypot(EX - 24, EY - 22) <= 10))) && passed;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Photo turned; corners given by hand

turned = {"top", 1; "bottom", -1; "left", 2};                               % rot90 of the reference photo -> side of the "3D" text
for t = 1:size(turned, 1)
    rt = OrganoidsFromPhoto(rot90(I2, turned{t,2}), struct('FIGURE', false, 'TEXT_SIDE', turned{t,1}));
    ok = rt.n_organoids == 2 && corner_error(rt, chip2) < 1.5;
    if ok
        ok = norm(rt.organoids(1).centroid - r2.organoids(1).centroid) < 1 && norm(rt.organoids(2).centroid - r2.organoids(2).centroid) < 1 ...
            && nnz(xor(rt.organoids(1).crop, r2.organoids(1).crop)) <= 0.03*nnz(r2.organoids(1).crop);
    end
    passed = report(sprintf('photo turned ("3D" text on the %s): same array and organoids', turned{t,1}), ok) && passed;
end
rm = OrganoidsFromPhoto(I1, struct('FIGURE', false, 'TEXT_SIDE', "right", 'ARRAY_CORNERS', chip.corners));
passed = report("corners of the array given by hand (ARRAY_CORNERS): same organoid", rm.array.source == "manual" ...
    && norm(rm.organoids(1).centroid - [33 30]) < 1.5) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Elongated organoid, forced number, small spot, no organoid

I3 = synthetic_chip([32 32 24], 'elongation', 1.6);                         % one elongated organoid (1.6 x 1)
r3 = OrganoidsFromPhoto(I3, quiet);
passed = report("an elongated organoid is not split", r3.n_organoids == 1) && passed;
r3b = OrganoidsFromPhoto(I3, struct('FIGURE', false, 'TEXT_SIDE', "right", 'N_ORGANOIDS', 2));
passed = report("N_ORGANOIDS = 2 splits it in two", r3b.n_organoids == 2) && passed;
I4 = synthetic_chip([30 34 20; 56 8 2.5]);                                  % an organoid and a small pale spot (radius 2.5 electrodes)
r4 = OrganoidsFromPhoto(I4, quiet);
passed = report("a small pale spot is not an organoid (MIN_ORGANOID)", r4.n_organoids == 1 && norm(r4.organoids(1).centroid - [30 34]) < 1.5) && passed;
try
    OrganoidsFromPhoto(synthetic_chip(zeros(0, 3)), quiet);
    stopped = false;
catch ME
    stopped = contains(ME.message, 'no organoid');
end
passed = report("a photo without organoid stops with an error", stopped) && passed;
[r5, figs] = OrganoidsFromPhoto(I2, struct('TEXT_SIDE', "right", 'SHOW_ORGANOID', 2, 'UNITS', [24 22; 42 44; 5 60], 'LABEL', "orgB"));
passed = report("SHOW_ORGANOID: 2 figures (copy of the photo, electrode map)", numel(figs) == 2 && all(isgraphics(figs)) && r5.n_organoids == 2) && passed;
close(figs);

I6 = synthetic_chip([24 40 13; 46 42 13; 35 21 13]);                       % three touching organoids
r6 = OrganoidsFromPhoto(I6, struct('FIGURE', false, 'TEXT_SIDE', "right", 'N_ORGANOIDS', 3));
ok = r6.n_organoids == 3;
if ok
    c6 = reshape([r6.organoids.centroid], 2, [])';
    truth = [24 40; 46 42; 35 21];
    ok = all(arrayfun(@(t) min(hypot(c6(:,1) - truth(t,1), c6(:,2) - truth(t,2))), 1:3) < 2);
end
passed = report("three touching organoids found with N_ORGANOIDS = 3 (centers within 2 electrodes)", ok) && passed;
[I9, chip9] = synthetic_chip([33 30 22], 'decoy', true);                   % the first guess of the place of the chip is wrong
r9 = OrganoidsFromPhoto(I9, quiet);
passed = report("the chip is found when the first guess (largest pale spot) is not on it (centers of the circles of the well)", ...
    hypot(r9.coarse.organoid(1) - chip9.decoy(1), r9.coarse.organoid(2) - chip9.decoy(2)) < 20 && corner_error(r9, chip9) < 1.5 ...
    && r9.n_organoids == 1 && norm(r9.organoids(1).centroid - [33 30]) < 1.5) && passed;
[I10, chip10] = synthetic_chip([33 30 22], 'side', 50);                     % a distant photo: array of 50 pixels
r10 = OrganoidsFromPhoto(I10, quiet);
passed = report(sprintf("a distant photo (array of 50 pixels): array found (corners within %.1f electrodes)", corner_error(r10, chip10)), ...
    corner_error(r10, chip10) < 2.5 && r10.n_organoids == 1) && passed;
I8 = synthetic_chip([48 40 18], 'strip', [65 67.5 2 50]);                   % a thin pale strip (2.5 electrodes) along the right edge of
one = struct('FIGURE', false, 'TEXT_SIDE', "right", 'N_ORGANOIDS', 1);    % the array, outside it, touching the organoid
r8 = OrganoidsFromPhoto(I8, one);
one.MIN_WIDTH = 0;      r8b = OrganoidsFromPhoto(I8, one);
along = EX >= 62 & EY <= 14;                                                % along the strip, beyond the margin of the organoid
passed = report("a thin pale strip along the edge of the array is left out (MIN_WIDTH): no electrode kept along it, all those under the organoid kept", ...
    ~any(r8.organoids(1).crop(along)) && all(r8.organoids(1).crop(hypot(EX - 48, EY - 40) <= 18 - 1.5)) && any(r8b.organoids(1).crop(along))) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Side of the "3D" text: landmarks of the board

[I7, chip7] = synthetic_chip([24 22 15; 42 44 15], 'landmarks', true);     % with the connector, the triangle and the text
turns = {"right", 0; "top", 1; "bottom", -1; "left", 2};                     % rot90 of the reference photo -> side of the "3D" text
for t = 1:size(turns, 1)
    rt = OrganoidsFromPhoto(rot90(I7, turns{t,2}), struct('FIGURE', false));
    ok = rt.text_side == turns{t,1} && rt.orientation.landmarks == turns{t,1} && rt.n_organoids == 2 && corner_error(rt, chip7) < 1.5;
    passed = report(sprintf('side of the "3D" text found from the landmarks (%s; %s)', turns{t,1}, rt.orientation.how), ok) && passed;
end
lastwarn('');
rw = OrganoidsFromPhoto(rot90(I7, 1), struct('FIGURE', false, 'TEXT_SIDE', "bottom"));     % (the text is on the top)
[msg, id] = lastwarn();
passed = report("a wrong side given: mismatch reported (warning), the given side used, the rest run", id == "OrganoidsFromPhoto:orientation" ...
    && contains(msg, "MISMATCH") && rw.text_side == "bottom" && ~isempty(rw.organoids)) && passed;
lastwarn('');
rr = OrganoidsFromPhoto(rot90(I7, 1), struct('FIGURE', false, 'TEXT_SIDE', "top"));
passed = report("the right side given: agrees with the landmarks, no warning", isempty(lastwarn()) && contains(rr.orientation.message, "agree")) && passed;
for t = 1:size(turns, 1)                                                   % no landmarks: side from the pattern of the chip
    lastwarn('');
    rn = OrganoidsFromPhoto(rot90(I2, turns{t,2}), struct('FIGURE', false));
    ok = rn.orientation.landmarks == "" && rn.text_side == turns{t,1} && rn.orientation.pattern == turns{t,1} && isempty(lastwarn()) ...
        && rn.n_organoids == 2 && corner_error(rn, chip2) < 1.5;
    passed = report(sprintf('no landmark: side of the "3D" text found from the pattern of the chip (%s)', turns{t,1}), ok) && passed;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. "SelectElectrodes", option 3 (cropping from the photo; its default)

folder = tempname;      mkdir(folder);
cleanup = onCleanup(@() rmdir(folder, 's'));
imwrite(rot90(I7, -1), fullfile(folder, 'chip.png'));                       % the "3D" text at the bottom, found from the landmarks
r2 = OrganoidsFromPhoto(I7, struct('FIGURE', false, 'TEXT_SIDE', "right"));
xy = [42 44; 40 40; 30 50; 24 22; 20 18; 33 33; 5 60; 60 5; 55 50];      % electrodes of the clusters (1-based)
spikeTimestamp = struct('timeStamp', num2cell(1:size(xy, 1)), 'cluster', num2cell(100 + (1:size(xy, 1))), ...
    'positionXY', num2cell((xy - 1)', 1));                                  % positionXY is 0-based
rec_conditions = [];    minimum_activity = [];  idx_abnormalActivity = [];  good_cluster = [];  n_discarded_HS = 0;
n_discarded_manual = 0;     n_clusterOriginal = size(xy, 1);    n_discarded_QCgeneral = 0;
save(fullfile(folder, 'syn_QCmanual.mat'), 'spikeTimestamp', 'rec_conditions', 'minimum_activity', 'idx_abnormalActivity', 'good_cluster', ...
    'n_discarded_HS', 'n_discarded_manual', 'n_clusterOriginal', 'n_discarded_QCgeneral');
evalc(['SelectElectrodes(struct(''system'', 0, ''folder'', folder, ''file_data'', ''syn'', ''organoid_label'', ''orgA'', ' ...
    '''file_photo'', ''chip.png'', ''photo_organoid'', 1))']);
out = load(fullfile(folder, 'syn_orgA_QCmanual.mat'));
kept = cell2mat(arrayfun(@(t) t.positionXY(:)' + 1, out.spikeTimestamp(:), 'UniformOutput', false));
expected = xy(r2.organoids(1).crop(sub2ind([N N], xy(:,1), xy(:,2))), :);
passed = report("SelectElectrodes (option 3, the default): the clusters on the electrodes kept", isequal(sortrows(kept), sortrows(expected)) ...
    && numel(out.spikeTimestamp_BeforeCropping) == size(xy, 1) && out.option == 3 && out.photo_organoid == 1) && passed;
passed = report("SelectElectrodes (option 3): outline, electrodes kept and the 2 figures saved", ~isempty(out.crop_outline) ...
    && ~isempty(out.crop_electrodes) && isfile(fullfile(folder, 'syn_orgA_crop_map.png')) && isfile(fullfile(folder, 'syn_orgA_crop_photo.png'))) && passed;
close all
% (November 2026) the table built from the folder lists the recording, not its cropped dataset; the cropping is run again from that row
% (run_pipeline, another CROPPING_SAFETY): the cropped dataset is written again
evalc('B = make_recordings_table(folder, fullfile(folder, ''built.xlsx''))');
passed = report("make_recordings_table: the recording, not its cropped dataset", isequal(B.file_name, "syn") && B.picture == "chip.png") && passed;
B.group = "g1";     B.organoid_label = "orgA";      B.picture_organoid = "1";
writetable(B, fullfile(folder, 'recordings.xlsx'));
evalc(['run_pipeline(fullfile(folder, ''recordings.xlsx''), ''stages'', "SelectElectrodes", ''stop_on_error'', 1, ' ...
    '''changes'', {''SelectElectrodes.CROPPING_SAFETY'', 1.4})']);
out2 = load(fullfile(folder, 'syn_orgA_QCmanual.mat'));
passed = report("cropping run again (run_pipeline, CROPPING_SAFETY 1.4): the cropped dataset written again, with more electrodes", ...
    size(out2.crop_electrodes, 1) > size(out.crop_electrodes, 1) && numel(out2.spikeTimestamp_BeforeCropping) == size(xy, 1)) && passed;
close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I, chip] = synthetic_chip(organoids, varargin)

% synthetic photo of a chip in the reference orientation ("3D" text on the right); organoids: one row per organoid, [x y radius] in
% electrodes; 'elongation': ratio of the axes of the organoids (along x); 'strip': a thin pale strip, [x1 x2 y1 y2] in electrodes
% (e.g., along an edge of the array, outside it, as the package can look in a photo); 'decoy': a pale spot on the well, outside the
% chip, larger than the organoids (chip.decoy: its center); 'side': side of the array (pixels; default 150). chip.corners: corners of the array (bottom-left, bottom-right,
% top-right, top-left), pixels

opts = struct('elongation', 1, 'landmarks', false, 'strip', [], 'decoy', false, 'side', 150);
for k = 1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
n = 64;     s = opts.side;                                                 % side of the array, pixels
H = 1100;   W = 1150;
R = 3.25*s;     C = [W/2 + 10, H/2 - 5];                                   % big circle of the well
if opts.landmarks                                                           % a smaller chip on a larger photo, with the board around
    s = 120;    R = 3.25*s;     H = 1600;   W = 1800;   C = [1100 800];
end
A = C + [-0.11 0.05]*s;                                                     % center of the array
b = 0.175*s;
[X, Y] = meshgrid(1:W, 1:H);
I = repmat(reshape([0.10 0.10 0.11], 1, 1, 3), H, W);                       % dark board
if opts.landmarks                                                           % gold landmarks (reference orientation, in radii of the circle)
    u = (X - C(1))/R;   v = (Y - C(2))/R;     gold = [0.86 0.70 0.30];
    I = paint(I, abs(v) <= 1.5 & ((u >= -2.60 & u <= -2.42) | (u >= -2.25 & u <= -2.07)), gold);    % connector: 2 columns, left
    I = paint(I, u >= 1.0 & v <= -1.3 & v >= -1.9 & (u - 1.0) <= (-1.3 - v), gold);                 % triangle: top-right corner
    for m = 0:4                                                             % "3D" text: marks along the right edge, below the center
        I = paint(I, u >= 1.40 & u <= 1.55 & v >= 0.40 + 0.2*m & v <= 0.52 + 0.2*m, gold);
    end
end
dc = hypot(X - C(1), Y - C(2));
I = paint(I, dc < R, [0.55 0.38 0.22]);                                     % floor of the well (gold traces)
I = paint(I, abs(dc - R) < 4, [0.12 0.10 0.10]);                            % big circle (dark line)
% package around the array (November 2026: as in the photos; electrodes from the center of the array, e = pixels per electrode)
e = s/n;    xe = (X - A(1))/e;  ye = (Y - A(2))/e;                          % (ye downwards)
lines_x = mod(xe, 2.2) < 1.1;   lines_y = mod(ye, 2.2) < 1.1;               % leads: gold lines every 2.2 electrodes
lead = [0.80 0.66 0.40];    plastic = [0.14 0.12 0.11];
side_leads = abs(xe) >= 67 & abs(xe) <= 105 & abs(ye) <= 54;               % left and right: horizontal lines
ends_leads = (ye <= -67 & ye >= -105 & abs(xe) <= 60) | (ye >= 67 & ye <= 105 & abs(xe) <= 40);   % above and below: vertical lines
in_plastic = abs(xe) <= 87 & abs(ye) <= 85;                                 % on the plastic of the package up to 87 electrodes, then on the board
I = paint(I, (side_leads | ends_leads) & in_plastic, plastic);
I = paint(I, (side_leads & lines_y) | (ends_leads & lines_x), lead);
I = paint(I, hypot(xe + 96, ye + 92) <= 9 | hypot(xe - 94, ye - 91) <= 9, [0.08 0.08 0.09]);   % 2 holes of the board (top-left, bottom-right)
rounded = @(half, rc) max(abs(xe) - (half - rc), 0).^2 + max(abs(ye) - (half - rc), 0).^2 <= rc^2;
I = paint(I, abs(xe) <= 67 & abs(ye) <= 67, [0.12 0.11 0.11]);              % package
I = paint(I, rounded(57, 10), [0.55 0.52 0.50]);                            % rim of its opening
I = paint(I, rounded(54, 8), [0.20 0.16 0.13]);                             % dark frame
I = paint(I, abs(xe) <= 32 & abs(ye) <= 32 + 64*0.175, [0.62 0.44 0.26]);   % chip: bands above and below the array
I = paint(I, xe >= -32 - 64*0.175 & xe <= -32 & abs(ye) <= 32, [0.48 0.33 0.18]);     % strip on the left of the array
I = paint(I, xe >= 43 & xe <= 47 & abs(ye) <= 31, [0.05 0.05 0.05]);       % slot on the right
I = paint(I, abs(xe) <= 32 & abs(ye) <= 32, [0.72 0.52 0.31]);              % electrode array
chip.decoy = [];
if opts.decoy                                                               % a pale spot on the well, larger than the organoids
    chip.decoy = C + [0, 2.3*s];
    I = paint(I, hypot(X - chip.decoy(1), Y - chip.decoy(2)) <= 0.5*s, [0.93 0.88 0.86]);
end
px = s/n;
for k = 1:size(organoids, 1)
    cx = A(1) - s/2 + (organoids(k,1) - 0.5)*px;    cy = A(2) + s/2 - (organoids(k,2) - 0.5)*px;
    rx = organoids(k,3)*px*sqrt(opts.elongation);   ry = organoids(k,3)*px/sqrt(opts.elongation);
    rho = hypot((X - cx)/rx, (Y - cy)/ry);
    I = paint(I, rho <= 1, reshape([0.93 0.86 0.83], 1, 1, 3) - 0.06*cos(6*atan2(Y - cy, X - cx)).*(rho > 0.5));    % pale, with lobes
end
if ~isempty(opts.strip)
    sx = A(1) - s/2 + (opts.strip(1:2) - 0.5)*px;    sy = A(2) + s/2 - (opts.strip(3:4) - 0.5)*px;
    I = paint(I, X >= sx(1) & X <= sx(2) & Y <= sy(1) & Y >= sy(2), [0.92 0.85 0.82]);
end
rng(3);
I = min(max(I + 0.015*randn(size(I)), 0), 1);                               % noise
chip.corners = A + s/2*[-1 1; 1 1; 1 -1; -1 -1];
chip.side = s;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = paint(I, mask, rgb)

% color rgb (3 values, or an image) on the pixels of mask
if numel(rgb) == 3, rgb = reshape(rgb, 1, 1, 3); end
I = I.*~mask + rgb.*mask;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = corner_error(r, chip)

% largest distance between the corners found and the true ones, in electrodes
e = max(hypot(r.array.corners(:,1) - chip.corners(:,1), r.array.corners(:,2) - chip.corners(:,2)))/(chip.side/64);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = report(name, ok)

if ok
    fprintf("  PASS  %s\n", name);
else
    fprintf("  FAIL  %s\n", name);
end

end
