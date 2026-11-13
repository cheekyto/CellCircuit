function [result, fig] = OrganoidsFromPhoto(photo, inputs)

% function to find, in a photo of a 3Brain chip (64 x 64 = 4096 electrodes), the electrode array and the organoids on it, and to give,
% for each organoid, the electrodes under it: its outline on the electrode map, the electrodes to keep for "SelectElectrodes", and the
% equivalent ellipse ("x0 y0 radiusx radiusy", as "crop_ellipse"), without setting the ellipse by hand (November 2026)
% Coordinates of the electrodes: as in "SelectElectrodes", x = positionXY(1)+1 and y = positionXY(2)+1, from 1 to 64. In the photo turned
% so that the "3D" text of the chip (outside the big circle) is on the right, read vertically (reference orientation), electrode (1,1) is
% the bottom-left corner of the array, x grows to the right and y upwards
% Only the photo is used (no activity of the recording). Steps:
%   1. organoids, coarse: pale (bright, unsaturated) regions surrounded by the tan surface of the chip (a first guess of where the chip is)
%   2. chip (November 2026, revised): the pattern of the edges around the electrode array, the same on every chip (the dark square of the
%      package and its 4 rows of leads, its rounded opening, the chip with the strip and the bands around the array, all centered on the
%      array), is searched at all sizes (array from 2% to 32% of the photo) around the organoid and around the centers of the concentric
%      circles of the well; it gives the place, size and tilt of the array. The pattern ("OrganoidsFromPhoto_chip.mat", next to this
%      function) is the mean of the edges of 49 photos of chips, straightened on the electrode coordinates
%   3. side of the "3D" text: the turn of the photo where the pattern of the chip (not the same in the 4 turns: strip on the left of the
%      array, slot on its right) fits best; it is checked against the landmarks of the board around the well (gold, outside the big
%      circle): the 2 columns of contacts of the connector (on the left of the well when the text is on the right), the triangle
%      (top-right corner) and the "3D" text itself (along the right edge, below the well), by priority connector > triangle > text.
%      When TEXT_SIDE is given, it is used, and a mismatch with the pattern or the landmarks is reported (warning)
%   4. electrode array: the pattern of the chip fitted in the photo turned to the reference orientation (size within 1%, tilt within
%      0.5 degree); in the reference orientation, the array (A) is the square of the chip surface with the strip on its left and the bands
%      above and below it. A poor fit (below MIN_FIT) is reported (warning): then give the corners of the array (ARRAY_CORNERS)
%   5. organoids, fine: pale (unsaturated) regions on the array, without their parts narrower than MIN_WIDTH (strips along the edges of
%      the array), at least MIN_ORGANOID of its area; touching organoids are split (Gaussian mixture of their pixels), and a part is split
%      again while the border between its two halves is a narrow neck (or until there are N_ORGANOIDS parts)
%   6. electrodes of each organoid: under it or within MARGIN x its radius from its outline (as CROPPING_SAFETY for an ellipse), but not
%      within SHARED_GAP electrodes of another organoid, nor of the border that two touching organoids share (conservative on that border)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   photo:          file of the photo (.jpg, .png, .tif...; .heic on a Mac), or the image itself (RGB array)
%   inputs:         (optional) structure that overwrites any of the PARAMETERS below, e.g., struct('TEXT_SIDE', "top")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   result:         structure with
%       .text_side          side of the photo where the "3D" text is (given, or found)
%       .orientation        how the side was set: .pattern (side where the pattern of the chip fits best), .landmarks (side given by
%                           the landmarks, "" if none visible), .how (which landmarks), .message (agreement or mismatch), .scores (fit
%                           of the pattern of the chip in the 4 turns)
%       .circle             [x y radius] of the big circle of the well, in pixels of the photo in the reference orientation ([] if the
%                           chip was not searched)
%       .array              .corners (4 x 2: bottom-left, bottom-right, top-right, top-left corners of the array, in pixels of the photo in
%                           the reference orientation), .side (pixels), .angle (degrees), .source ("pattern" or "manual"), .fit (fit
%                           of the pattern of the chip, normalized correlation of the edges: usually 0.5 to 0.75)
%       .n_organoids        number of organoids found
%       .organoids          structure array, one per organoid, numbered from the top of the array (reference orientation) down:
%           .outline            outline of the organoid, [x y] in electrodes (closed polygon; it can go beyond the array)
%           .centroid           [x y] center of the organoid, in electrodes
%           .area               area of the organoid, in electrodes (squares of 1 x 1)
%           .radius             equivalent radius, sqrt(area/pi), in electrodes
%           .ellipse            equivalent ellipse [x0 y0 radiusx radiusy] (centroid, 2 x standard deviation along x and y), as "crop_ellipse"
%           .crop               MATRIX_SIZE x MATRIX_SIZE logical, crop(x, y) = true for the electrodes to keep (x, y from 1 to MATRIX_SIZE)
%           .electrodes         [x y] of the electrodes to keep
%       .view               the photo around the array, turned and straightened so that its axes are the electrode coordinates:
%                           .image (uint8 RGB), .x and .y (electrode coordinates of its columns and rows; .y grows upwards)
%       .labels             organoid of each pixel of .view.image (0 = none)
%   fig:            figure: the photo (whole, then around the array, with the grid of electrodes and the organoids) and the electrodes
%                   kept for each organoid ([] when FIGURE = false). With SHOW_ORGANOID = k: 2 figures, [photo, map]: (1) a copy of the photo
%                   with the outline of organoid k and the electrodes kept (dashed), (2) the electrode map: electrodes with active units
%                   (UNITS) kept or not kept, electrodes without a unit, with the same outlines

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TEXT_SIDE = "";                 % side of the photo where the "3D" text of the chip is ("right", "top", "left" or "bottom"); the photo is
                                % turned so that the text is on the right (reference orientation). ="" (or "auto"): found from the photo
N_ORGANOIDS = 0;                % number of organoids on the array (=0, found from the photo: touching organoids are split when the border
                                % between them is a narrow neck; when they overlap a lot in the photo, give their number)
MARGIN = 0.10;                  % electrodes kept around an organoid, as a share of its radius (0.10 = the default CROPPING_SAFETY of 1.10)
SHARED_GAP = 2;                 % electrodes within this distance (in electrodes) of another organoid are not kept (border shared by two organoids)
MAX_ANGLE = 10;                 % largest tilt of the chip in the photo, in degrees
MIN_FIT = 0.30;                 % (November 2026) a fit of the pattern of the chip below this is reported (array maybe misplaced; in the
                                % photos, 0.5 to 0.75 when the array was right, below 0.25 when it was not)
ARRAY_CORNERS = [];             % (optional) 4 x 2 corners of the array set by hand (bottom-left, bottom-right, top-right, top-left; pixels of
                                % the photo in the reference orientation): the array is not searched
MATRIX_SIZE = 64;               % dimension of the electrode array (64*64=4096)
VIEW_BORDER = 0.20;             % photo shown around the array, as a share of its side
MIN_ORGANOID = 0.04;            % smallest organoid, as a share of the area of the array (smaller pale spots are reflections or dust)
MIN_ORGANOID_SHARE = 0.15;      % smallest organoid, as a share of the largest one
ORGANOID_BORDER = 6;            % organoids are outlined up to this distance (electrodes) outside the array (beyond: the rim of the package)
MIN_WIDTH = 8;                  % (November 2026) the parts of the pale regions narrower than this (electrodes) are not organoid: strips
                                % along the edges of the array, where the package or the bands look pale (=0: kept)
FIGURE = true;                  % =true, draw the figure
SHOW_ORGANOID = 0;              % >0: instead, the 2 figures of the cropping of this organoid (number), as saved by "SelectElectrodes"
UNITS = zeros(0, 2);            % (with SHOW_ORGANOID) [x y] electrodes of the active units (positionXY + 1), shown on the electrode map
LABEL = "";                     % (with SHOW_ORGANOID) name of the organoid, for the titles

if nargin < 2, inputs = struct(); end
override_inputs(inputs);
TEXT_SIDE = lower(string(TEXT_SIDE));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Photo, and the turns that can put it in the reference orientation

I = read_photo(photo);                                                      % as the photo is shown
sides_all = ["right", "top", "bottom", "left"];     turns_all = [0, -1, 1, 2];  % turn (rot90) that puts the "3D" text on the right
given = strlength(TEXT_SIDE) > 0 && TEXT_SIDE ~= "auto";
if given && ~ismember(TEXT_SIDE, sides_all)
    error('TEXT_SIDE must be "right", "top", "left", "bottom" or "" (side of the photo where the "3D" text of the chip is; "": found)');
end
[H, W, ~] = size(I);

result = struct('text_side', TEXT_SIDE, 'orientation', [], 'circle', [], 'array', [], 'n_organoids', 0, 'organoids', [], 'view', [], ...
    'labels', [], 'coarse', []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Organoids, coarse (photo reduced to about 1000 pixels): a first guess of the place of the chip

f = max(1, round(max(H, W)/1000));                                          % reduction factor
Ic = block_mean(I, f);
[pale_c, tan_c] = pale_and_tan(Ic);
[L, n] = label_mask(pale_c);
area = accumarray(L(L > 0), 1, [n 1]);
best = 0;   best_area = 0;  best_score = 0;
min_area = max(30, 2e-4*numel(L));      max_area = 0.10*numel(L);
for c = find(area >= min_area & area <= max_area)'
    [r, q] = find(L == c);
    r_eq = sqrt(area(c)/pi);
    ring = ring_mask(size(L), r, q, max(2, round(0.4*r_eq)));
    tan_share = mean(tan_c(ring));                                          % share of tan surface around the region
    if tan_share >= 0.75 && area(c) > best_area                             % the largest pale region surrounded by the tan surface
        best_area = area(c);    best_score = tan_share;     best = c;
    end
end
c0 = [];                                                                    % center of the organoid (coarse pixels; none found: [])
if best > 0
    [r, q] = find(L == best);
    c0 = [mean(q) mean(r)];
    result.coarse = struct('organoid', [(c0 - 0.5)*f + 0.5, sqrt(area(best)/pi)*f], 'tan_share', best_score);    % (pixels of the photo as shown)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2-3. Chip: the pattern around the electrode array (photo as shown); side of the "3D" text (pattern, landmarks of the board)

Y = mean(I, 3);
tm = chip_template();
chip = struct('center', [], 'side', NaN, 'angle', 0, 'fit', -inf);
if isempty(ARRAY_CORNERS) || ~given
    seeds = ring_centers(I);                                                % centers of the concentric circles of the well
    if ~isempty(c0), seeds = [(c0 - 0.5)*f + 0.5; seeds]; end              % and the organoid
    pyr = pyramid(Y);
    for j = 1:size(seeds, 1)
        c = find_chip(pyr, seeds(j,:), tm.sym, MAX_ANGLE);
        if c.fit > chip.fit, chip = c; end
    end
end
circle = [];    landmarks = "";     how = "";   pattern = "";
fits = nan(1, 4);   fitted = cell(1, 4);
if isfinite(chip.fit)
    cc = (chip.center - 0.5)/f + 0.5;                                        % center of the chip (reduced photo)
    circle = find_circle(Ic, cc);                                           % big circle of the well, around the chip
    if isempty(circle) || hypot(circle(1) - cc(1), circle(2) - cc(2)) > chip.side/f
        circle = [cc, 3.25*chip.side/f];                                    % (not found: its usual size)
    end
    [landmarks, how] = find_landmarks(Ic, circle);                          % side of the "3D" text given by the landmarks of the board
    circle(1:2) = (circle(1:2) - 0.5)*f + 0.5;   circle(3) = circle(3)*f;    % pixels of the photo (as shown)
    for k = 1:4                                                             % fit of the pattern in each turn of the photo
        fitted{k} = fit_pattern(rot90(Y, turns_all(k)), turn_point(chip.center, turns_all(k), H, W), chip.side, chip.angle, tm.fine, true);
        fits(k) = fitted{k}.fit;
    end
    [~, k] = max(fits);     pattern = sides_all(k);
end
orientation = struct('pattern', pattern, 'landmarks', landmarks, 'how', how, 'message', "", 'scores', []);
sure = false;                                                               % the pattern fits clearly best in one turn
if any(isfinite(fits))
    orientation.scores = array2table(fits, 'VariableNames', cellstr(sides_all));
    sf = sort(fits, 'descend');
    sure = sf(1) >= MIN_FIT && sf(1) - sf(2) >= 0.05*sf(1);                % (in the photos, the best fit was 21% to 42% above the next)
end
lm_text = "";
if strlength(landmarks) > 0, lm_text = sprintf('the landmarks of the board (%s) put it on the %s', how, landmarks); end
if given
    side_used = TEXT_SIDE;
    bad = strings(0, 1);    agree = strings(0, 1);
    if strlength(pattern) > 0
        if pattern ~= TEXT_SIDE && sure, bad(end+1) = sprintf('the pattern of the chip around the array puts it on the %s', pattern); %#ok<AGROW>
        elseif pattern == TEXT_SIDE, agree(end+1) = "the pattern of the chip"; end %#ok<AGROW>
    end
    if strlength(landmarks) > 0                                             % (the landmarks count when the pattern does not settle it)
        if landmarks ~= TEXT_SIDE && ~(pattern == TEXT_SIDE && sure), bad(end+1) = lm_text; %#ok<AGROW>
        elseif landmarks == TEXT_SIDE, agree(end+1) = sprintf("the landmarks of the board (%s)", how); end %#ok<AGROW>
    end
    if ~isempty(bad)
        orientation.message = sprintf(['ORIENTATION MISMATCH: the "3D" text was given on the %s of the photo, but %s. The given side is ' ...
            'used: check TEXT_SIDE and the figure'], TEXT_SIDE, strjoin(bad, ', and '));
        warning('OrganoidsFromPhoto:orientation', '%s', orientation.message);
    elseif isempty(agree)
        orientation.message = sprintf('"3D" text on the %s (given; not checked)', TEXT_SIDE);
    else
        orientation.message = sprintf('"3D" text on the %s (given; %s agree)', TEXT_SIDE, strjoin(agree, " and "));
        if strlength(landmarks) > 0 && landmarks ~= TEXT_SIDE
            orientation.message = [orientation.message, sprintf(' (%s, likely misread)', lm_text)];
        end
    end
elseif strlength(pattern) > 0
    side_used = pattern;
    if sure
        if strlength(landmarks) == 0
            orientation.message = sprintf('"3D" text on the %s (found from the pattern of the chip around the array; no landmark of the board visible)', pattern);
        elseif landmarks == pattern
            orientation.message = sprintf('"3D" text on the %s (found from the pattern of the chip around the array; the landmarks agree: %s)', pattern, how);
        else
            orientation.message = sprintf('"3D" text on the %s (found from the pattern of the chip around the array; %s, likely misread)', pattern, lm_text);
        end
    elseif strlength(landmarks) > 0                                         % the pattern does not settle it: the landmarks
        side_used = landmarks;
        orientation.message = sprintf(['"3D" text on the %s, found from the landmarks of the board (%s): the pattern of the chip does not ' ...
            'settle it (best on the %s). Check the figure, or give the side of the "3D" text'], landmarks, how, pattern);
        warning('OrganoidsFromPhoto:orientation', '%s', orientation.message);
    else
        orientation.message = sprintf(['"3D" text on the %s, found from the pattern of the chip, which does not settle it (no landmark of ' ...
            'the board visible). Check the figure, or give the side of the "3D" text'], pattern);
        warning('OrganoidsFromPhoto:orientation', '%s', orientation.message);
    end
else
    error(['OrganoidsFromPhoto: the chip was not found in the photo; give the side of the "3D" text (TEXT_SIDE) and the corners of the ' ...
        'electrode array (ARRAY_CORNERS)']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Electrode array: the pattern of the chip fitted in the photo turned to the reference orientation

k = find(sides_all == side_used);   k_turn = turns_all(k);
I = rot90(I, k_turn);
if ~isempty(circle), circle(1:2) = turn_point(circle(1:2), k_turn, H, W); end
result.text_side = side_used;
if ~isempty(ARRAY_CORNERS)
    corners = ARRAY_CORNERS;    source = "manual";      match = NaN;
    side = mean(sqrt(sum(diff([corners; corners(1,:)]).^2, 2)));
    angle = atan2d(corners(2,2) - corners(1,2), corners(2,1) - corners(1,1));
else
    a = fit_pattern(rot90(Y, k_turn), fitted{k}.center, fitted{k}.side, fitted{k}.angle, tm.fine, false);
    corners = a.corners;    side = a.side;  angle = a.angle;    match = a.fit;  source = "pattern";
    if match < MIN_FIT
        warning('OrganoidsFromPhoto:array', ['The electrode array may be misplaced (fit of the pattern of the chip %.2f, usually 0.5 to ' ...
            '0.75): check the figure, or give the corners of the array (ARRAY_CORNERS)'], match);
    end
end
clear Y
result.circle = circle;
result.orientation = orientation;
[H, W, ~] = size(I);
result.array = struct('corners', corners, 'side', side, 'angle', angle, 'source', source, 'fit', match);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Organoids, fine: photo straightened on the electrode coordinates

px = max(2, min(5, round(side/MATRIX_SIZE)));                               % pixels per electrode in the view (as in the photo, from 2 to 5)
step = 1/px;
lo = 0.5 - VIEW_BORDER*MATRIX_SIZE;     hi = MATRIX_SIZE + 0.5 + VIEW_BORDER*MATRIX_SIZE;
xv = lo + step/2 : step : hi;           yv = xv;                            % electrode coordinates of the view
[XV, YV] = meshgrid(xv, yv);
[UV, VV] = electrode_to_photo(XV, YV, corners, MATRIX_SIZE);
view = zeros([size(XV) 3]);
smooth = max(0, side/MATRIX_SIZE/px/2);                                     % smoothing before sampling (when the photo has more pixels)
for ch = 1:3
    plane = I(:,:,ch);
    if smooth > 0.6, plane = gauss_smooth(plane, smooth); end
    view(:,:,ch) = interp2(plane, UV, VV, 'linear', 0);
end
result.view = struct('image', uint8(255*min(max(view, 0), 1)), 'x', xv, 'y', yv);

[V, mn, sat] = color_features(view);
on_array = XV >= 0.5 & XV <= MATRIX_SIZE + 0.5 & YV >= 0.5 & YV <= MATRIX_SIZE + 0.5;
lit = V > 0.2;                                                              % (the saturation of dark pixels is not reliable)
near_array = XV >= 0.5 - ORGANOID_BORDER & XV <= MATRIX_SIZE + 0.5 + ORGANOID_BORDER & YV >= 0.5 - ORGANOID_BORDER & YV <= MATRIX_SIZE + 0.5 + ORGANOID_BORDER;
th = otsu(sat(on_array & lit));
if ~(mean(sat(on_array & lit & sat < th)) < 0.85*mean(sat(on_array & lit & sat >= th)))     % (no organoid: no part of the array
    % clearly less saturated; in the photos, the organoids are 0.4 to 0.75 times as saturated as the chip)
    error('OrganoidsFromPhoto: no organoid found on the electrode array');
end
pale = lit & near_array & sat < th;                                         % organoids: much less saturated than the tan surface of the chip
r_close = max(1, round(1.5*px));                                            % seams between the lobes of an organoid (1.5 electrodes)
pale = erode_mask(dilate_mask(pale, r_close), r_close);
pale = majority(pale, max(1, round(px/2)));
if MIN_WIDTH > 0                                                            % (November 2026) thin pale strips along the edges of the
    body = open_mask(pale, max(1, round(MIN_WIDTH/2*px)));                  % array (package, bands) removed, before the holes are filled
    pale = pale & dilate_mask(body, px);                                    % (a strip and the organoid enclose part of the array);
end                                                                         % elsewhere, the outline is kept (within 1 electrode)
pale = fill_holes(pale);
[Lp, n] = label_mask(pale);
sizes = accumarray(Lp(Lp > 0), 1, [n 1]);
on_share = accumarray(Lp(Lp > 0), double(on_array(Lp > 0)), [n 1])./max(sizes, 1);
% organoids: at least MIN_ORGANOID of the array (and MIN_ORGANOID_SHARE of the largest one), mostly on the array (smaller pale spots:
% reflections, dust)
keep = sizes >= max(MIN_ORGANOID*(MATRIX_SIZE*px)^2, MIN_ORGANOID_SHARE*max(sizes)) & on_share >= 0.5;
labels = zeros(size(Lp));
parts = {};
for c = find(keep)'
    parts = [parts, {find(Lp == c)}]; %#ok<AGROW>
end
if isempty(parts)
    error('OrganoidsFromPhoto: no organoid found on the electrode array');
end
% number of organoids given: the N_ORGANOIDS largest regions; a single region holding several organoids is split directly into them
if N_ORGANOIDS > 0
    if numel(parts) > N_ORGANOIDS
        [~, o] = sort(cellfun(@numel, parts), 'descend');
        parts = parts(o(1:N_ORGANOIDS));
    elseif numel(parts) == 1 && N_ORGANOIDS >= 3
        parts = split_k(parts{1}, N_ORGANOIDS, size(Lp));
    end
end
% split touching organoids (number given: the part with the narrowest neck, until there are N_ORGANOIDS)
neck = nan(1, numel(parts));    halves = cell(1, numel(parts));
for p = 1:numel(parts), [neck(p), halves{p}] = try_split(parts{p}, size(Lp)); end
while true
    if N_ORGANOIDS > 0
        if numel(parts) >= N_ORGANOIDS, break, end
        [~, p] = min(neck);                                                 % the part with the narrowest neck
        if isnan(neck(p)), break, end
    else
        p = find(neck < 0.90, 1);                                           % a narrow neck: two organoids
        if isempty(p) || numel(parts) >= 8, break, end
    end
    new = halves{p};
    parts(p) = [];      neck(p) = [];       halves(p) = [];
    for h = 1:2
        parts{end+1} = new{h}; %#ok<AGROW>
        [neck(end+1), halves{end+1}] = try_split(new{h}, size(Lp)); %#ok<AGROW>
    end
end

if N_ORGANOIDS == 0                                                         % (after splitting; when the number is given, all are kept)
    min_part = MIN_ORGANOID*(MATRIX_SIZE*px)^2;
    parts = parts(cellfun(@(idx) numel(idx) >= min_part && mean(on_array(idx)) >= 0.5, parts));
end
if isempty(parts)
    error('OrganoidsFromPhoto: no organoid found on the electrode array');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Electrodes of each organoid

n_org = numel(parts);
cy = cellfun(@(idx) mean(YV(idx)), parts);      cx = cellfun(@(idx) mean(XV(idx)), parts);
[~, order] = sortrows([-round(cy(:)/4), cx(:)]);                                  % numbered from the top down (then from the left)
parts = parts(order);
for k = 1:n_org, labels(parts{k}) = k; end
[EX, EY] = ndgrid(1:MATRIX_SIZE, 1:MATRIX_SIZE);                            % electrodes, (x, y)
dist = inf(MATRIX_SIZE, MATRIX_SIZE, n_org);                                % distance from each electrode to each organoid, in electrodes
organoids = struct('outline', cell(1, n_org), 'centroid', [], 'area', [], 'radius', [], 'ellipse', [], 'crop', [], 'electrodes', []);
for k = 1:n_org
    mask = labels == k;
    xk = XV(mask);  yk = YV(mask);
    organoids(k).centroid = [mean(xk) mean(yk)];
    organoids(k).area = numel(xk)*step^2;
    organoids(k).radius = sqrt(organoids(k).area/pi);
    organoids(k).ellipse = [mean(xk) mean(yk) 2*std(xk) 2*std(yk)];
    organoids(k).outline = mask_outline(mask, xv, yv);
    edge_pix = mask & ~erode_mask(mask);
    dist(:,:,k) = min_distance([EX(:) EY(:)], [XV(edge_pix) YV(edge_pix)], MATRIX_SIZE);
    under = interp2(xv, yv, double(mask), EX, EY, 'nearest', 0) > 0;
    d = dist(:,:,k);    d(under) = 0;      dist(:,:,k) = d;
end
for k = 1:n_org
    others = dist(:,:,setdiff(1:n_org, k));
    crop = dist(:,:,k) <= MARGIN*organoids(k).radius;
    if ~isempty(others)
        crop = crop & min(others, [], 3) > SHARED_GAP;                      % not within SHARED_GAP of another organoid
        mask = labels == k;
        border = mask & (circshift(labels, 1, 1) > 0 | circshift(labels, -1, 1) > 0 | circshift(labels, 1, 2) > 0 | circshift(labels, -1, 2) > 0) ...
            & ~(circshift(labels, 1, 1) == k & circshift(labels, -1, 1) == k & circshift(labels, 1, 2) == k & circshift(labels, -1, 2) == k);
        border = border & (ismember(circshift(labels, 1, 1), setdiff(1:n_org, k)) | ismember(circshift(labels, -1, 1), setdiff(1:n_org, k)) ...
            | ismember(circshift(labels, 1, 2), setdiff(1:n_org, k)) | ismember(circshift(labels, -1, 2), setdiff(1:n_org, k)));
        if any(border(:))                                                   % border shared with a touching organoid: not within SHARED_GAP
            crop = crop & min_distance([EX(:) EY(:)], [XV(border) YV(border)], MATRIX_SIZE) > SHARED_GAP;   % of it (also around its ends)
        end
    end
    organoids(k).crop = crop;
    organoids(k).electrodes = [EX(crop) EY(crop)];
end
result.n_organoids = n_org;
result.organoids = organoids;
result.labels = labels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure

fig = [];
if FIGURE && SHOW_ORGANOID > 0                                              % cropping of one organoid (as saved by "SelectElectrodes")
    k = SHOW_ORGANOID;
    if k > n_org
        error('OrganoidsFromPhoto: organoid %d asked, but %d organoid(s) were found in the photo', k, n_org);
    end
    colors = lines(max(n_org, 1));      o = organoids(k);
    name = string(LABEL);   if strlength(name) == 0, name = "organoid " + k; end
    % (1) copy of the photo with the cropping
    fig = figure('Name', 'Cropping on the photo', 'Color', 'w', 'Position', [50 50 1150 540]);
    subplot(1, 2, 1)
    small = block_mean_uint8(I, f);
    image(small, 'XData', ((1:size(small,2)) - 0.5)*f + 0.5, 'YData', ((1:size(small,1)) - 0.5)*f + 0.5);
    axis image off;     hold on
    plot(corners([1:4 1],1), corners([1:4 1],2), 'y-', 'LineWidth', 1.2);
    [u, v] = electrode_to_photo(o.outline(:,1), o.outline(:,2), corners, MATRIX_SIZE);
    plot(u, v, '-', 'Color', colors(k,:), 'LineWidth', 1.2);
    title(sprintf('Photo ("3D" text on the right)'), 'FontSize', 9);
    subplot(1, 2, 2)
    draw_view(result, colors, true, k);
    title(sprintf('%s: organoid %d of %d (dashed: electrodes kept)', name, k, n_org), 'FontSize', 9, 'Interpreter', 'none');
    % (2) electrode map with the active units
    fig(2) = figure('Name', 'Cropping on the electrode map', 'Color', 'w', 'Position', [80 80 640 560]);
    units = round(UNITS);
    units = units(all(units >= 1 & units <= MATRIX_SIZE, 2), :);
    category = zeros(MATRIX_SIZE);                                          % rows = y: 0 = no active unit, 1 = not kept, 2 = kept
    for i = 1:size(units, 1)
        category(units(i,2), units(i,1)) = 1 + o.crop(units(i,1), units(i,2));
    end
    image(1:MATRIX_SIZE, 1:MATRIX_SIZE, category + 1, 'CDataMapping', 'scaled');
    colormap(gca, [0.92 0.92 0.92; 0.25 0.45 0.85; 0.95 0.45 0.10]);
    set(gca, 'YDir', 'normal');     axis image;     hold on
    for j = setdiff(1:n_org, k)
        plot(organoids(j).outline(:,1), organoids(j).outline(:,2), '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8);
    end
    plot(o.outline(:,1), o.outline(:,2), 'k-', 'LineWidth', 1.5);
    draw_crop_outline(o.crop, [0 0 0]);
    xlim([0.5 MATRIX_SIZE + 0.5]);  ylim([0.5 MATRIX_SIZE + 0.5]);
    set(gca, 'XTick', [1 8:8:MATRIX_SIZE], 'YTick', [1 8:8:MATRIX_SIZE], 'FontSize', 8);
    xlabel('x (electrodes)');   ylabel('y (electrodes)');
    clim([0.5 3.5]);    cb = colorbar;  cb.Ticks = [1 2 3];
    cb.TickLabels = {'no active unit', 'active, not kept', 'active, kept'};
    n_kept = nnz(category == 2);    n_units = nnz(category > 0);
    title({sprintf('%s: electrodes with active units kept (%d of %d)', name, n_kept, n_units), ...
        'solid: outline of the organoid; dashed: electrodes kept'}, 'FontSize', 9, 'Interpreter', 'none');
elseif FIGURE
    fig = figure('Name', 'Organoids on the electrode array', 'Color', 'w', 'Position', [50 50 1500 520]);
    colors = lines(max(n_org, 1));
    % whole photo (reference orientation) with the circle and the array
    subplot(1, 3, 1)
    small = block_mean_uint8(I, f);
    image(small, 'XData', ((1:size(small,2)) - 0.5)*f + 0.5, 'YData', ((1:size(small,1)) - 0.5)*f + 0.5);
    axis image off;     hold on
    if ~isempty(circle)
        t = linspace(0, 2*pi, 200);
        plot(circle(1) + circle(3)*cos(t), circle(2) + circle(3)*sin(t), 'c-', 'LineWidth', 1);
    end
    plot(corners([1:4 1],1), corners([1:4 1],2), 'y-', 'LineWidth', 1.5);
    title(sprintf('Photo turned ("3D" text was on the %s)\narray: %.0f px, %.1f deg (%s, fit %.2f)', result.text_side, side, angle, source, ...
        result.array.fit), 'FontSize', 9);
    % photo straightened on the electrode coordinates
    subplot(1, 3, 2)
    draw_view(result, colors, true, 0);
    title(sprintf('Organoids on the array (%d)', n_org), 'FontSize', 9);
    % electrodes kept for each organoid
    subplot(1, 3, 3)
    draw_crops(result, colors);
    txt = arrayfun(@(k) sprintf('%d: ellipse [%.0f %.0f %.0f %.0f], %d electrodes', k, organoids(k).ellipse, nnz(organoids(k).crop)), 1:n_org, 'UniformOutput', false);
    title([{'Electrodes kept'}; txt(:)], 'FontSize', 8);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   LOCAL FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = read_photo(photo)
% RGB image in [0, 1], turned as the photo is shown (orientation of the camera, EXIF)
if isnumeric(photo) || islogical(photo)
    I = to_double(photo);   return
end
photo = char(photo);
if ~isfile(photo), error('OrganoidsFromPhoto: the photo "%s" does not exist', photo); end
[~, ~, ext] = fileparts(photo);
if any(strcmpi(ext, {'.heic', '.heif'}))                                    % MATLAB reads HEIC only with a support package: converted by macOS
    if ~ismac, error('OrganoidsFromPhoto: save the photo "%s" as .jpg (HEIC photos are read only on a Mac)', photo); end
    tmp = [tempname '.jpg'];
    [status, msg] = system(sprintf('sips -s format jpeg "%s" --out "%s"', photo, tmp));
    if status ~= 0, error('OrganoidsFromPhoto: the photo "%s" could not be converted to .jpg (%s)', photo, msg); end
    I = to_double(imread(tmp));     delete(tmp);
    I = rot90(I, heif_rotation(photo));                                     % (the conversion keeps the pixels as stored)
    return
end
I = to_double(imread(photo));
try
    info = imfinfo(photo);
    if isfield(info, 'Orientation') && ~isempty(info(1).Orientation)
        switch info(1).Orientation                                          % EXIF orientation: how to show the stored pixels
            case 2, I = fliplr(I);
            case 3, I = rot90(I, 2);
            case 4, I = flipud(I);
            case 5, I = permute(I, [2 1 3]);
            case 6, I = rot90(I, -1);
            case 7, I = rot90(permute(I, [2 1 3]), 2);
            case 8, I = rot90(I, 1);
        end
    end
catch
end
end

function k = heif_rotation(photo)
% number of 90-degree anticlockwise turns that show a HEIC/HEIF photo as it was taken: property "irot" of the image, in the box "meta" at
% the start of the file (a photo taken with the phone held sideways is stored turned, e.g., irot = 3: 90 degrees clockwise to show it)
k = 0;
fid = fopen(photo, 'r', 'ieee-be');
if fid < 0, return, end
closer = onCleanup(@() fclose(fid));
pos = 0;
for n = 1:20                                                                % top-level boxes: "ftyp", "meta", "mdat"...
    if fseek(fid, pos, 'bof') ~= 0, return, end
    sz = fread(fid, 1, 'uint32');   type = fread(fid, [1 4], 'uint8=>char');
    if isempty(sz) || numel(type) < 4, return, end
    head = 8;
    if sz == 1, sz = fread(fid, 1, 'uint64');   head = 16; end
    if strcmp(type, 'meta')
        b = fread(fid, min(sz - head, 2^20), 'uint8=>uint8')';
        i = strfind(char(b), 'irot');
        if ~isempty(i) && i(1) + 4 <= numel(b)
            k = mod(double(b(i(1) + 4)), 4);
        end
        return
    end
    if sz < 8, return, end
    pos = pos + sz;
end
end

function I = to_double(I)
if isa(I, 'uint8'), I = double(I)/255;
elseif isa(I, 'uint16'), I = double(I)/65535;
else, I = double(I);
end
if size(I, 3) == 1, I = repmat(I, 1, 1, 3); end
I = I(:,:,1:3);
end

function B = block_mean(I, f)
% mean of blocks of f x f pixels
if f == 1, B = I; return, end
[H, W, C] = size(I);
h = floor(H/f);     w = floor(W/f);
B = zeros(h, w, C);
for ch = 1:C
    X = I(1:h*f, 1:w*f, ch);
    X = reshape(sum(reshape(X, f, []), 1), h, w*f);                         % sum over the rows of each block
    X = reshape(sum(reshape(X', f, []), 1), w, h)';                         % then over its columns
    B(:,:,ch) = X/f^2;
end
end

function B = block_mean_uint8(I, f)
B = uint8(255*min(max(block_mean(I, f), 0), 1));
end

function [V, mn, sat, hue] = color_features(I)
R = I(:,:,1);   G = I(:,:,2);   B = I(:,:,3);
V = max(I, [], 3);      mn = min(I, [], 3);
sat = (V - mn)./max(V, 1e-6);
if nargout > 3
    d = max(V - mn, 1e-6);
    hue = zeros(size(V));
    is_r = V == R;      is_g = ~is_r & V == G;      is_b = ~is_r & ~is_g;
    hue(is_r) = 60*mod((G(is_r) - B(is_r))./d(is_r), 6);
    hue(is_g) = 60*((B(is_g) - R(is_g))./d(is_g) + 2);
    hue(is_b) = 60*((R(is_b) - G(is_b))./d(is_b) + 4);
    hue(hue > 330) = hue(hue > 330) - 360;                                  % reds close to 0
end
end

function [pale, tan] = pale_and_tan(I)
% pale: bright and unsaturated (organoids); tan: orange surface of the chip
[V, ~, sat, hue] = color_features(I);
pale = sat < 0.45 & V > 0.45;
tan = sat >= 0.35 & hue >= -10 & hue <= 55 & V > 0.15;
end

function t = otsu(v)
% threshold of Otsu (two classes)
v = v(:);   v = v(isfinite(v));
if isempty(v) || min(v) == max(v), t = mean(v); return, end
edges = linspace(min(v), max(v), 257);
h = histcounts(v, edges)';      p = h/sum(h);       c = (edges(1:end-1) + edges(2:end))'/2;
w = cumsum(p);      m = cumsum(p.*c);       mt = m(end);
sb = (mt*w - m).^2./max(w.*(1 - w), 1e-12);
[~, i] = max(sb);   t = c(i);
end

function [L, n] = label_mask(mask)
% connected regions (4 neighbors)
[H, W] = size(mask);
L = zeros(H, W);    n = 0;
idx = find(mask);
queue = zeros(numel(idx), 1);
for s = idx'
    if L(s) > 0, continue, end
    n = n + 1;      L(s) = n;
    queue(1) = s;   head = 1;   tail = 1;
    while head <= tail
        p = queue(head);    head = head + 1;
        r = mod(p - 1, H) + 1;
        nb = [p - 1*(r > 1), p + 1*(r < H), p - H*(p > H), p + H*(p <= H*(W - 1))];
        for q = nb
            if q ~= p && mask(q) && L(q) == 0
                L(q) = n;   tail = tail + 1;    queue(tail) = q;
            end
        end
    end
end
end

function ring = ring_mask(sz, r, q, d)
% pixels within d of the region (rows r, columns q), outside it
region = false(sz);     region(sub2ind(sz, r, q)) = true;
rows = max(1, min(r) - d):min(sz(1), max(r) + d);
cols = max(1, min(q) - d):min(sz(2), max(q) + d);
grown = false(sz);
grown(rows, cols) = dilate_mask(region(rows, cols), d);
ring = grown & ~region;
end

function M = dilate_mask(M, d)
[a, b] = meshgrid(-d:d);
M = conv2(double(M), double(a.^2 + b.^2 <= d^2), 'same') > 0.5;
end

function M = erode_mask(M, d)
if nargin < 2, d = 1; end
[a, b] = meshgrid(-d:d);
K = double(a.^2 + b.^2 <= d^2);
M = conv2(double(M), K, 'same') > sum(K(:)) - 0.5;
end

function M = open_mask(M, d)
M = dilate_mask(erode_mask(M, d), d);
end

function M = majority(M, d)
[a, b] = meshgrid(-d:d);
K = double(a.^2 + b.^2 <= d^2);
M = conv2(double(M), K/sum(K(:)), 'same') > 0.5;
end

function M = fill_holes(M)
% background regions that do not touch the border become part of the mask
[Lb, ~] = label_mask(~M);
border = unique([Lb(1,:), Lb(end,:), Lb(:,1)', Lb(:,end)']);
M = M | (Lb > 0 & ~ismember(Lb, border));
end

function X = gauss_smooth(X, sigma)
k = ceil(3*sigma);  g = exp(-(-k:k).^2/(2*sigma^2));    g = g/sum(g);
X = conv2(g, g, X, 'same')./conv2(g, g, ones(size(X)), 'same');
end

function [gx, gy] = gradients(Y)
s = [1 2 1]'/4;     dk = [1 0 -1]/2;
gx = conv2(s, dk, Y, 'same');       gy = conv2(dk', s', Y, 'same');         % brighter to the right / downwards = positive
gx(:, [1 end]) = 0;     gy([1 end], :) = 0;
end

function circle = find_circle(I, c)
% big circle of the well, [x y radius] in pixels of I ([] if not found). Hough transform of the edges, centers near the organoid c;
% the well shows two concentric circles, the inner one (the big circle) and an outer ring about 1.4 times larger: when the strongest
% circle has a strong one about 1.4 times smaller, that one is the big circle
circle = [];
Y = gauss_smooth(mean(I, 3), 1.5);
[gx, gy] = gradients(Y);
g = hypot(gx, gy);
[H, W] = size(Y);
m = min(H, W);
r_min = 0.08*m;     r_max = 0.6*m;      c_max = 0.12*m;                     % radius of the circle, distance of its center from the organoid
[q, r] = meshgrid(1:W, 1:H);
dc = hypot(q - c(1), r - c(2));
zone = dc >= r_min - c_max & dc <= r_max + c_max;
gz = sort(g(zone));
strong = zone & g > gz(ceil(0.9*numel(gz)));                                % the 10% strongest edges
x = q(strong);  y = r(strong);  w = g(strong);
ux = gx(strong)./w;     uy = gy(strong)./w;
n_acc = 2*ceil(c_max) + 1;      acc = zeros(n_acc);
for rad = r_min:r_max
    for sgn = [-1 1]
        ax = round(x + sgn*rad*ux - c(1)) + ceil(c_max) + 1;
        ay = round(y + sgn*rad*uy - c(2)) + ceil(c_max) + 1;
        ok = ax >= 1 & ax <= n_acc & ay >= 1 & ay <= n_acc;
        acc = acc + accumarray([ay(ok) ax(ok)], 1, [n_acc n_acc]);
    end
end
acc = gauss_smooth(acc, 1.5);
[~, i] = max(acc(:));
[ay, ax] = ind2sub(size(acc), i);
center = [ax ay] - ceil(c_max) - 1 + c;
% radius: radial edges at the same distance from the center (per unit of circumference)
d = hypot(x - center(1), y - center(2));
radial = abs((x - center(1)).*ux + (y - center(2)).*uy)./max(d, 1e-6) > 0.95;
sel = radial & d >= r_min & d <= r_max;
if nnz(sel) < 20, return, end
bins = (floor(r_min):ceil(r_max) + 1)';
h = accumarray(floor(d(sel)) - bins(1) + 1, w(sel), [numel(bins) 1])./bins;
h = conv(h, [1 2 3 2 1]'/9, 'same');
peaks = find(h(2:end-1) >= h(1:end-2) & h(2:end-1) >= h(3:end)) + 1;
if isempty(peaks), return, end
[~, i] = max(h(peaks));     p1 = peaks(i);
inner = peaks(bins(peaks) >= bins(p1)/1.5 & bins(peaks) <= bins(p1)/1.3 & h(peaks) >= 0.4*h(p1));
if ~isempty(inner)
    [~, i] = max(h(inner));     p1 = inner(i);
end
rad = bins(p1);
% least-squares circle through the edges at that radius
sel = radial & abs(d - rad) <= max(2, 0.03*rad);
A = [x(sel) y(sel) ones(nnz(sel), 1)];
b = x(sel).^2 + y(sel).^2;
p = A\b;
xc = p(1)/2;    yc = p(2)/2;    R = sqrt(p(3) + xc^2 + yc^2);
if hypot(xc - center(1), yc - center(2)) > 0.05*rad || abs(R - rad) > 0.05*rad
    xc = center(1);     yc = center(2);     R = rad;
end
circle = [xc yc R];
end

function tm = chip_template()
% pattern of the edges around the electrode array, from "OrganoidsFromPhoto_chip.mat" (next to this function): mean absolute gradients
% (along x and along y) of 49 photos of chips straightened on the electrode coordinates (1.2 sides around the array, centered on it), at
% 1 pixel per electrode (fine: reference orientation) and 0.5 (coarse: mean of the 4 turns, the same in any orientation)
persistent cache
if isempty(cache)
    file = fullfile(fileparts(mfilename('fullpath')), 'OrganoidsFromPhoto_chip.mat');
    if ~isfile(file)
        error('OrganoidsFromPhoto: the pattern of the chip is missing (%s; it comes with the code)', file);
    end
    P = load(file);
    cache.fine = make_template(double(P.fine_gx), double(P.fine_gy), double(P.fine_px), "fine");
    cx = double(P.coarse_gx);   cy = double(P.coarse_gy);                   % |gx| and |gy| swap at each turn of 90 degrees
    cache.sym = make_template((cx + rot90(cy, 1) + rot90(cx, 2) + rot90(cy, 3))/4, (cy + rot90(cx, 1) + rot90(cy, 2) + rot90(cx, 3))/4, ...
        double(P.coarse_px), "sym");
end
tm = cache;
end

function t = make_template(TX, TY, px, id)
% template for the normalized correlation (2 channels, |gx| and |gy|): zero mean, weighted; the inside of the array (organoids) is left out
[h, w] = size(TX);
c = [(w + 1)/2, (h + 1)/2];                                                 % center of the array (pixels of the template)
[X, Yy] = meshgrid(((1:w) - c(1))/px, ((1:h) - c(2))/px);                  % electrodes from the center of the array
Wm = double(~(abs(X) <= 28 & abs(Yy) <= 28));
Wm([1 end], :) = 0;     Wm(:, [1 end]) = 0;
Nw = sum(Wm(:));    Tz = cell(1, 2);    tn = 0;     Tc = {TX, TY};
for k = 1:2
    mu = sum(Wm(:).*Tc{k}(:))/Nw;
    Tz{k} = Wm.*(Tc{k} - mu);   tn = tn + sum(Wm(:).*(Tc{k}(:) - mu).^2);
end
t = struct('id', id, 'size', [h w], 'W', Wm, 'Tz', {Tz}, 'tn', tn, 'Nw', Nw, 'px', px, 'center', c);
end

function C = ring_centers(I)
% centers of the concentric circles of the well, [x y] pixels of I: Hough votes of the strong edges for all the radii from 0.08 to 0.55 x
% the photo (photo reduced to about 500 pixels); the 3 best peaks, at least 0.1 x the photo apart
[H, W, ~] = size(I);
f = max(1, round(max(H, W)/500));
Y = gauss_smooth(block_mean(mean(I, 3), f), 1.5);
[gx, gy] = gradients(Y);    g = hypot(gx, gy);
[h, w] = size(Y);   m = min(h, w);
gs = sort(g(:));    strong = g > gs(ceil(0.9*numel(gs)));
[r, q] = find(strong);  wg = g(strong);     ux = gx(strong)./wg;    uy = gy(strong)./wg;
acc = zeros(h, w);
for rad = round(0.08*m):round(0.55*m)
    for sgn = [-1 1]
        ax = round(q + sgn*rad*ux);     ay = round(r + sgn*rad*uy);
        ok = ax >= 1 & ax <= w & ay >= 1 & ay <= h;
        acc = acc + accumarray([ay(ok) ax(ok)], 1, [h w]);
    end
end
acc = gauss_smooth(acc, 2);
[qq, rr] = meshgrid(1:w, 1:h);
C = zeros(0, 2);
for k = 1:3
    [~, i] = max(acc(:));   [ay, ax] = ind2sub(size(acc), i);
    C(end+1, :) = ([ax ay] - 0.5)*f + 0.5; %#ok<AGROW>
    acc(hypot(qq - ax, rr - ay) < 0.1*m) = 0;
end
end

function pyr = pyramid(Y)
% Y reduced by 1, 2, 4, 8... (mean of blocks), down to about 128 pixels
pyr = struct('Y', {Y}, 'f', 1);
while min(size(pyr(end).Y)) >= 256
    pyr(end+1) = struct('Y', {block_mean(pyr(end).Y, 2)}, 'f', 2*pyr(end).f); %#ok<AGROW>
end
end

function chip = find_chip(pyr, seed, t, max_angle)
% coarse search of the pattern of the chip (orientation-free template t) at all the sizes of the array from 2% to 32% of the photo (steps
% of 4%), its center within 0.8 x its side of seed; the tilt is the dominant orientation of the edges around seed (modulo 90 degrees)
[H, W] = size(pyr(1).Y);    m = min(H, W);
l = find([pyr.f] <= max(1, 0.4*m/400), 1, 'last');
B = pyr(l).Y;   fl = pyr(l).f;  sl = (seed - 0.5)/fl + 0.5;     rad = 0.2*m/fl;
x0 = max(1, round(sl(1) - rad));    x1 = min(size(B, 2), round(sl(1) + rad));
y0 = max(1, round(sl(2) - rad));    y1 = min(size(B, 1), round(sl(2) + rad));
[gx, gy] = gradients(B(y0:y1, x0:x1));
ang = angle_of(sum((gx(:).^2 + gy(:).^2).*exp(4i*atan2(gy(:), gx(:)))))/4*180/pi;
if abs(ang) > max_angle, ang = 0; end
chip = struct('center', seed, 'side', NaN, 'angle', ang, 'fit', -inf);
for s = exp(log(0.02*m):log(1.04):log(0.32*m))
    [v, c] = ncc_at(pyr, seed, s, ang, 0.8, t);
    if v > chip.fit, chip.fit = v;  chip.center = c;    chip.side = s; end
end
end

function a = fit_pattern(Y, center, side, ang0, t, quick)
% the pattern of the chip (template t, reference orientation) fitted around center (pixels of Y): sizes within 8% of side (steps of 1%),
% tilt within 1.5 degrees of ang0 (steps of 0.5), center within 0.12 x side; quick: sizes within 4% (3) and the tilt ang0, center within
% 0.15 x side (fit of a turn of the photo). a: .fit, .center, .side, .angle, .corners (bottom-left, bottom-right, top-right, top-left)
if quick
    angles = ang0;  sizes = side*[0.96 1 1.04];     search = 0.15;
else
    angles = ang0 + (-1.5:0.5:1.5);     sizes = side*exp(log(0.92):0.01:log(1.08));     search = 0.12;
end
r = ceil(2.6*max(sizes));
x0 = max(1, round(center(1)) - r);  x1 = min(size(Y, 2), round(center(1)) + r);
y0 = max(1, round(center(2)) - r);  y1 = min(size(Y, 1), round(center(2)) + r);
pyr = pyramid(Y(y0:y1, x0:x1));
c_loc = center - [x0 y0] + 1;
a = struct('fit', -inf, 'center', center, 'side', side, 'angle', ang0, 'corners', []);
for ang = angles
    for s = sizes
        [v, c] = ncc_at(pyr, c_loc, s, ang, search, t);
        if v > a.fit, a.fit = v;    a.center = c + [x0 y0] - 1;     a.side = s;     a.angle = ang; end
    end
end
ca = cosd(a.angle);     sa = sind(a.angle);
cr = a.side/2*[-1 1; 1 1; 1 -1; -1 -1];                                     % (y downwards)
a.corners = a.center + [ca*cr(:,1) - sa*cr(:,2), sa*cr(:,1) + ca*cr(:,2)];
end

function [best, center] = ncc_at(pyr, seed, s, ang, search_sides, t)
% best normalized correlation of the template t with the photo (pyramid pyr), for an array of side s (pixels) tilted by ang (degrees),
% its center within search_sides x s of seed; center: where it is best (sub-pixel)
persistent spectra                                                          % Fourier transforms of the templates (by template and size)
tp = 64*t.px;   step = s/tp;                                                % pixels of the photo per pixel of the template
th = t.size(1);     tw = t.size(2);
search = search_sides*tp;
ph = ceil(max(th, tw)/2 + search) + 1;                                      % half size of the patch sampled around the seed
l = find([pyr.f] <= step, 1, 'last');   if isempty(l), l = 1; end
B = pyr(l).Y;   fl = pyr(l).f;  stl = step/fl;  sl = (seed - 0.5)/fl + 0.5;
hx = ceil((ph + 2)*stl*1.5) + 2;
x0 = max(1, round(sl(1)) - hx);     x1 = min(size(B, 2), round(sl(1)) + hx);
y0 = max(1, round(sl(2)) - hx);     y1 = min(size(B, 1), round(sl(2)) + hx);
best = -inf;    center = seed;
if x1 - x0 < 3 || y1 - y0 < 3, return, end
B = B(y0:y1, x0:x1);
if stl > 1.2, B = gauss_smooth(B, 0.5*stl); end
ca = cosd(ang);     sa = sind(ang);
[U, V] = meshgrid(-ph:ph, -ph:ph);
P = interp2(B, sl(1) - x0 + 1 + stl*(ca*U - sa*V), sl(2) - y0 + 1 + stl*(sa*U + ca*V), 'linear', NaN);
valid = ~isnan(P);
if nnz(valid) < 0.5*numel(P), return, end
P(~valid) = mean(P(valid));
[gx, gy] = gradients(P);    Pc = {abs(gx), abs(gy)};
n = 2*ph + 1;   Hf = 2^nextpow2(n);
key = sprintf('%s_%d', t.id, Hf);
if isempty(spectra), spectra = containers.Map(); end
if ~isKey(spectra, key)
    spectra(key) = {conj(fft2(t.W, Hf, Hf)), conj(fft2(t.Tz{1}, Hf, Hf)), conj(fft2(t.Tz{2}, Hf, Hf))};
end
FS = spectra(key);  FW = FS{1};
num = 0;    vr = 0;
for c = 1:2
    F = fft2(Pc{c}, Hf, Hf);
    num = num + real(ifft2(F.*FS{c + 1}));
    S1 = real(ifft2(F.*FW));    S2 = real(ifft2(fft2(Pc{c}.^2, Hf, Hf).*FW));
    vr = vr + max(S2 - S1.^2/t.Nw, 0);
end
nu = n - th + 1;    nv = n - tw + 1;
ncc = num(1:nu, 1:nv)./sqrt(max(vr(1:nu, 1:nv), 1e-12)*t.tn);
[cu, cv] = meshgrid((1:nv) + t.center(1) - 1 - (ph + 1), (1:nu) + t.center(2) - 1 - (ph + 1));    % center of the array from the seed
ncc(hypot(cu, cv) > search) = -inf;
[best, j] = max(ncc(:));
[iu, iv] = ind2sub(size(ncc), j);   du = 0;     dv = 0;                     % sub-pixel: parabola through the neighbors
if iu > 1 && iu < nu && all(isfinite(ncc([iu-1 iu+1], iv)))
    d = ncc(iu-1, iv) - 2*ncc(iu, iv) + ncc(iu+1, iv);    if d < 0, dv = 0.5*(ncc(iu-1, iv) - ncc(iu+1, iv))/d; end
end
if iv > 1 && iv < nv && all(isfinite(ncc(iu, [iv-1 iv+1])))
    d = ncc(iu, iv-1) - 2*ncc(iu, iv) + ncc(iu, iv+1);    if d < 0, du = 0.5*(ncc(iu, iv-1) - ncc(iu, iv+1))/d; end
end
uc = cu(j) + du;    vc = cv(j) + dv;
center = seed + step*[ca*uc - sa*vc, sa*uc + ca*vc];
end

function [side, how] = find_landmarks(I, circle)
% side of the photo where the "3D" text of the chip is, from the gold landmarks of the board outside the big circle ([x y radius], pixels
% of I): the connector (2 long columns of gold contacts, on the left of the well when the text is on the right), the triangle (top-right
% corner) and the "3D" text (along the right edge, below the center of the well). side = "" when none is visible; how: the landmarks used
side = "";  how = "";
[V, ~, sat, hue] = color_features(I);
[H, W] = size(V);
[q, r] = meshgrid(1:W, 1:H);
u = (q - circle(1))/circle(3);  v = (r - circle(2))/circle(3);              % pixels of the photo, in radii of the circle (v downwards)
gold = hue >= 35 & hue <= 65 & sat > 0.40 & V > 0.45 & hypot(u, v) > 1.55;
unit = 1/circle(3)^2;                                                       % area of a pixel, in squared radii
votes = strings(0, 1);  names = strings(0, 1);
opposite = ["left", "bottom", "top", "right"];
% 1. connector: gold mass along one side, 1.8 to 3 radii from the center (bands: right, top, bottom, left)
bands = {abs(v) <= 1.7 & u >= 1.8 & u <= 3, abs(u) <= 1.7 & v <= -1.8 & v >= -3, abs(u) <= 1.7 & v >= 1.8 & v <= 3, ...
    abs(v) <= 1.7 & u <= -1.8 & u >= -3};
mass = cellfun(@(b) nnz(gold & b), bands)*unit;
[m, i] = max(mass);     rest = mass;    rest(i) = [];
connector = false(H, W);
if m >= 0.10 && m >= 4*max(rest)
    votes(end+1) = opposite(i);     names(end+1) = "connector";
    connector = bands{i};
end
% 2. triangle: a compact gold region (at least 0.03 squared radii) in a corner, 0.7 to 3 radii from the center along both axes
[L, n] = label_mask(gold & ~connector);
for c = 1:n
    idx = find(L == c);
    if numel(idx)*unit < 0.03, continue, end
    uc = mean(u(idx));  vc = mean(v(idx));
    e = eig(cov([u(idx) v(idx)]));
    if abs(uc) >= 0.7 && abs(vc) >= 0.7 && abs(uc) <= 3 && abs(vc) <= 3 && max(e)/max(min(e), eps) <= 6
        quadrant = ["top", "right"; "left", "bottom"];                      % (vc > 0) + 1, (uc > 0) + 1
        votes(end+1) = quadrant((vc > 0) + 1, (uc > 0) + 1);    names(end+1) = "triangle"; %#ok<AGROW>
        gold(idx) = false;
        break
    end
end
% 3. "3D" text: small gold marks grouped along an edge of the board, 1.3 to 2.4 radii from the center
marks = dilate_mask(gold & ~connector & hypot(u, v) <= 2.4 & max(abs(u), abs(v)) >= 1.0, max(1, round(0.08*circle(3))));
[L, n] = label_mask(marks);
best = 0;   best_mass = 0.01;
for c = 1:n
    idx = find(L == c & gold);
    if numel(idx)*unit > best_mass, best_mass = numel(idx)*unit; best = c; end
end
if best > 0
    idx = find(L == best & gold);
    uc = mean(u(idx));  vc = mean(v(idx));
    [E, D] = eig(cov([u(idx) v(idx)]));     [lam, j] = sort(diag(D), 'descend');
    if lam(1) >= 6*lam(2)                                                   % a line of marks
        if abs(E(1, j(1))) < abs(E(2, j(1)))                                % vertical: right (below the center) or left (above)
            if uc > 0 && vc > 0, votes(end+1) = "right"; names(end+1) = "text"; end
            if uc < 0 && vc < 0, votes(end+1) = "left"; names(end+1) = "text"; end
        else                                                                % horizontal: bottom (left of the center) or top (right)
            if vc > 0 && uc < 0, votes(end+1) = "bottom"; names(end+1) = "text"; end
            if vc < 0 && uc > 0, votes(end+1) = "top"; names(end+1) = "text"; end
        end
    end
end
if isempty(votes), return, end
side = votes(1);                                                            % by priority: connector > triangle > text
how = strjoin(names, ", ");
if any(votes ~= side)
    how = how + " (they disagree: " + strjoin(names + " " + votes, ", ") + ")";
end
end

function p = turn_point(p, k, H, W)
% pixel (x, y) of a photo of H x W pixels, in the photo turned by rot90(photo, k)
switch mod(k, 4)
    case 1, p = [p(2), W + 1 - p(1)];                                       % 90 degrees counterclockwise
    case 2, p = [W + 1 - p(1), H + 1 - p(2)];
    case 3, p = [H + 1 - p(2), p(1)];                                       % 90 degrees clockwise
end
end

function a = angle_of(z)
a = atan2(imag(z), real(z));
end

function [u, v] = electrode_to_photo(x, y, corners, n)
% pixels of the photo (reference orientation) of the electrode coordinates (x, y); corners: bottom-left, bottom-right, top-right, top-left
% (edges of the array: x = 0.5 and n + 0.5, y = 0.5 and n + 0.5)
a = (x - 0.5)/n;    b = (y - 0.5)/n;                                        % 0 to 1 across the array
u = (1 - a).*(1 - b)*corners(1,1) + a.*(1 - b)*corners(2,1) + a.*b*corners(3,1) + (1 - a).*b*corners(4,1);
v = (1 - a).*(1 - b)*corners(1,2) + a.*(1 - b)*corners(2,2) + a.*b*corners(3,2) + (1 - a).*b*corners(4,2);
end

function [neck, halves] = try_split(idx, sz)
% splits the pixels idx in two (Gaussian mixture); neck = length of the border between the halves / diameter of the smaller half
neck = nan;     halves = {};
if numel(idx) < 50, return, end
[r, q] = ind2sub(sz, idx);
X = [q r];
mu = mean(X);   C = cov(X);
[E, D] = eig(C);    [lam, j] = max(diag(D));    e = E(:, j)';
m = [mu - sqrt(lam)*e; mu + sqrt(lam)*e];                                   % start: on both sides along the long axis
S = cat(3, C/2, C/2);   p = [0.5 0.5];
for it = 1:100
    lik = zeros(size(X, 1), 2);
    for c = 1:2
        d = X - m(c,:);     Si = inv(S(:,:,c));
        lik(:, c) = p(c)*exp(-0.5*sum((d*Si).*d, 2))/sqrt(det(S(:,:,c)));
    end
    resp = lik./max(sum(lik, 2), realmin);
    nk = sum(resp, 1);
    if any(nk < 10), return, end
    p = nk/sum(nk);
    for c = 1:2
        m(c,:) = resp(:,c)'*X/nk(c);
        d = X - m(c,:);
        S(:,:,c) = (d.*resp(:,c))'*d/nk(c) + eye(2)*0.5;
    end
end
[~, a] = max(resp, [], 2);
if any(accumarray(a, 1, [2 1]) < 0.1*numel(idx)), return, end
M = zeros(sz);      M(idx) = a;
% keep each half in one piece (pixels cut off go to the other half)
for c = 1:2
    [Lc, nc] = label_mask(M == c);
    if nc > 1
        cnt = accumarray(Lc(Lc > 0), 1);    [~, big] = max(cnt);
        M(Lc > 0 & Lc ~= big) = 3 - c;
    end
end
A = M == 1;     B = M == 2;
border = A & (circshift(B, 1, 1) | circshift(B, -1, 1) | circshift(B, 1, 2) | circshift(B, -1, 2));
diam = 2*sqrt(min(nnz(A), nnz(B))/pi);
neck = nnz(border)/diam;
halves = {find(A), find(B)};
end

function parts = split_k(idx, k, sz)
% splits the pixels idx in k parts (Gaussian mixture of their coordinates, started from k pixels far from each other); each part is kept
% in one piece (pieces cut off go to the part with the nearest center)
[r, q] = ind2sub(sz, idx);
X = [q r];      n = size(X, 1);
start = zeros(k, 1);
[~, start(1)] = max(sum((X - mean(X)).^2, 2));                              % farthest pixels from each other
dmin = sum((X - X(start(1),:)).^2, 2);
for c = 2:k
    [~, start(c)] = max(dmin);      dmin = min(dmin, sum((X - X(start(c),:)).^2, 2));
end
m = X(start, :);    C = repmat(cov(X)/k, 1, 1, k);   p = ones(1, k)/k;
for it = 1:100
    lik = zeros(n, k);
    for c = 1:k
        d = X - m(c,:);
        lik(:, c) = p(c)*exp(-0.5*sum((d/C(:,:,c)).*d, 2))/sqrt(det(C(:,:,c)));
    end
    resp = lik./max(sum(lik, 2), realmin);
    nk = max(sum(resp, 1), 1e-9);      p = nk/sum(nk);
    for c = 1:k
        m(c,:) = resp(:,c)'*X/nk(c);
        d = X - m(c,:);     C(:,:,c) = (d.*resp(:,c))'*d/nk(c) + eye(2)*0.5;
    end
end
[~, a] = max(resp, [], 2);
M = zeros(sz);      M(idx) = a;
for c = 1:k                                                                 % one piece per part
    [Lc, nc] = label_mask(M == c);
    if nc > 1
        cnt = accumarray(Lc(Lc > 0), 1);    [~, big] = max(cnt);
        for piece = setdiff(1:nc, big)
            pix = find(Lc == piece);    [pr, pq] = ind2sub(sz, pix);
            dist = sum((mean([pq pr], 1) - m).^2, 2);  dist(c) = inf;
            [~, to] = min(dist);    M(pix) = to;
        end
    end
end
parts = arrayfun(@(c) find(M == c), 1:k, 'UniformOutput', false);
parts = parts(~cellfun(@isempty, parts));
end

function d = min_distance(P, Q, n)
% for each point of P (rows), distance to the nearest point of Q, as an n x n matrix (P = electrodes in the order of ndgrid)
d = inf(size(P, 1), 1);
for s = 1:2000:size(Q, 1)
    q = Q(s:min(s + 1999, size(Q, 1)), :);
    d = min(d, min(hypot(P(:,1) - q(:,1)', P(:,2) - q(:,2)'), [], 2));
end
d = reshape(d, n, n);
end

function xy = mask_outline(mask, xv, yv)
% outline of the largest piece of the mask, [x y] (closed polygon), in the coordinates xv (columns) and yv (rows)
C = contourc(xv, yv, double(mask), [0.5 0.5]);
xy = zeros(0, 2);   i = 1;
while i < size(C, 2)
    m = C(2, i);
    seg = C(:, i + 1:i + m)';
    if size(seg, 1) > size(xy, 1), xy = seg; end
    i = i + m + 1;
end
end

function draw_view(result, colors, show_numbers, highlight)
% photo straightened on the electrode coordinates, with the array, the grid and the organoids (current axes); highlight > 0: only that
% organoid with its electrodes kept, the others thin
v = result.view;    n = size(result.organoids(1).crop, 1);
image(v.x, v.y, v.image);   set(gca, 'YDir', 'normal');     axis image;     hold on
for g = 0.5:8:n + 0.5
    plot([g g], [0.5 n + 0.5], '-', 'Color', [1 1 1 0.35], 'LineWidth', 0.5);
    plot([0.5 n + 0.5], [g g], '-', 'Color', [1 1 1 0.35], 'LineWidth', 0.5);
end
plot([0.5 n+0.5 n+0.5 0.5 0.5], [0.5 0.5 n+0.5 n+0.5 0.5], 'y-', 'LineWidth', 1.5);
for k = 1:numel(result.organoids)
    o = result.organoids(k);
    if highlight > 0 && k ~= highlight
        plot(o.outline(:,1), o.outline(:,2), '-', 'Color', [1 1 1 0.6], 'LineWidth', 0.8);
        continue
    end
    plot(o.outline(:,1), o.outline(:,2), '-', 'Color', colors(k,:), 'LineWidth', 1.5);
    draw_crop_outline(o.crop, colors(k,:));
    if show_numbers
        text(o.centroid(1), o.centroid(2), sprintf('%d', k), 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 14, ...
            'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 1 0.6]);
    end
end
set(gca, 'XTick', [1 8:8:n], 'YTick', [1 8:8:n], 'FontSize', 8);
xlabel('x (electrodes)');   ylabel('y (electrodes)');
end

function draw_crops(result, colors)
% electrodes kept for each organoid (current axes)
n = size(result.organoids(1).crop, 1);
map = zeros(n, n, 3) + 0.92;
for k = 1:numel(result.organoids)
    crop = result.organoids(k).crop';                                       % rows = y
    for ch = 1:3
        plane = map(:,:,ch);    plane(crop) = 0.55 + 0.45*colors(k, ch);  map(:,:,ch) = plane;
    end
end
image(1:n, 1:n, map);   set(gca, 'YDir', 'normal');     axis image;     hold on
for k = 1:numel(result.organoids)
    o = result.organoids(k);
    plot(o.outline(:,1), o.outline(:,2), '-', 'Color', colors(k,:), 'LineWidth', 1.5);
    text(o.centroid(1), o.centroid(2), sprintf('%d', k), 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'center');
end
xlim([0.5 n + 0.5]);    ylim([0.5 n + 0.5]);
set(gca, 'XTick', [1 8:8:n], 'YTick', [1 8:8:n], 'FontSize', 8);
xlabel('x (electrodes)');   ylabel('y (electrodes)');
end

function draw_crop_outline(crop, color)
% outline of the electrodes kept (squares of 1 x 1 around each electrode), dashed
n = size(crop, 1);
up = kron(double(crop'), ones(4));                                          % 4 sub-pixels per electrode, rows = y
xs = linspace(0.5 + 1/8, n + 0.5 - 1/8, 4*n);
C = contourc(xs, xs, up, [0.5 0.5]);
i = 1;
while i < size(C, 2)
    m = C(2, i);
    plot(C(1, i + 1:i + m), C(2, i + 1:i + m), '--', 'Color', color, 'LineWidth', 1);
    i = i + m + 1;
end
end
