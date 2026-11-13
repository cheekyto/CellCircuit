function passed = test_waveform_features(code_dir)

% Function-level test of "WaveformFeatures" and "loadKilosortTemplates" on synthetic templates (new, November 2026)
% Like "test_cNEDetection", it does not compare two versions of the code: it checks that
%   1. the measures of a sampled template (19.75 kHz, as the recordings) match the same measures computed on the exact waveform
%      (sum of Gaussians, evaluated every microsecond, with its exact slope)
%   2. a Gaussian trough gives the values known from theory (half-width, slopes, phase-plot area, spectrum); a waveform still rising
%      (or falling) at the edge of the template gets the edge flags and empty values where the peak is outside the template; small ripples
%      are not counted as peaks
%   3. the measures do not depend on the amplitude, a positive waveform is flipped and flagged as non-somatic, and bombcell type 3 is
%      flagged too; cell array input = matrix input
%   4. the area of the phase plot (polygon of the loop) equals integral(slope^2) / integral(waveform^2)
%   5. edge cases: padding that is not zero (error), a template of zeros or NaN, no units
%   6. "loadKilosortTemplates" returns the template of each unit in the order of spikeTimestamp, from "waveforms.mat" or from the text
%      of "metadata.csv", and stops if a unit is missing
% No reference recording is needed
%
% Example (from a terminal):
%   matlab -batch "addpath('<repo>/tests'); test_waveform_features('<repo>')"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dir:             folder with the code to be tested

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% passed:               =true if all the checks passed (each check is printed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(char(code_dir));
FS = 19753.775390625;                                                       % sampling rate of the recordings, in Hz
N_PAD = 21;     N_KEEP = 61;    ALIGN = 21;                                 % Kilosort 2.5 template: 21 zeros + 61 samples, trough at sample 21
t = ((1:N_KEEP) - ALIGN) / FS * 1000;                                       % time of each sample, in msec (0 = alignment sample)
passed = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Sampled templates vs the exact waveform

% shapes: [amplitude, center (msec), width (msec)] of each Gaussian; the trough is not on a sample (sub-sample offsets)
shapes = { [0.25 -0.40 0.15; -1 0.012 0.09; 0.30 0.45 0.15], ...           % narrow
           [0.15 -0.50 0.20; -1 0.020 0.16; 0.25 1.20 0.35], ...           % broad
           [0.05 -0.60 0.20; -1 -0.017 0.12; 0.45 0.70 0.25] };            % small first peak, large second peak
names = {'narrow', 'broad', 'large P2'};
W = zeros(numel(shapes), N_PAD + N_KEEP);
for s=1:numel(shapes)
    W(s,:) = [zeros(1,N_PAD) 0.37*gauss_sum(shapes{s}, t)];                 % arbitrary amplitude, as in the templates
end
F = WaveformFeatures(W, FS, [1 1 2]);
for s=1:numel(shapes)
    R = exact_measures(shapes{s}, t(1), t(end));
    tol_time = 0.006;                                                       % msec (about 1/10 of a sample)
    ok = abs(F.troughToPeak(s) - R.troughToPeak) < tol_time && abs(F.halfWidth(s) - R.halfWidth) < tol_time && ...
         abs(F.slopeTroughToPeak(s) - R.slopeTroughToPeak) < tol_time;
    passed = report(sprintf("%s: times (TTP %.3f vs %.3f, half-width %.3f vs %.3f, slope TTP %.3f vs %.3f msec)", names{s}, ...
             F.troughToPeak(s), R.troughToPeak, F.halfWidth(s), R.halfWidth, F.slopeTroughToPeak(s), R.slopeTroughToPeak), ok) && passed;
    ok = abs(F.peakAfterRatio(s) - R.peakAfterRatio) < 0.005 && abs(F.peakBeforeRatio(s) - R.peakBeforeRatio) < 0.005 && ...
         abs(F.asymmetry(s) - R.asymmetry) < 0.01;
    passed = report(sprintf("%s: ratios (P2 %.3f vs %.3f, P1 %.3f vs %.3f, asymmetry %.3f vs %.3f)", names{s}, F.peakAfterRatio(s), ...
             R.peakAfterRatio, F.peakBeforeRatio(s), R.peakBeforeRatio, F.asymmetry(s), R.asymmetry), ok) && passed;
    rel = @(a,b) abs(a-b)/abs(b);
    ok = rel(F.maxDescent(s), R.maxDescent) < 0.01 && rel(F.maxRecovery(s), R.maxRecovery) < 0.01 && ...
         rel(F.descentRecoveryRatio(s), R.descentRecoveryRatio) < 0.01 && rel(F.phasePlotArea(s), R.phasePlotArea) < 0.01;
    passed = report(sprintf("%s: slopes within 1%% (descent %.2f vs %.2f, recovery %.2f vs %.2f, area %.1f vs %.1f)", names{s}, ...
             F.maxDescent(s), R.maxDescent, F.maxRecovery(s), R.maxRecovery, F.phasePlotArea(s), R.phasePlotArea), ok) && passed;
end
ok = F.troughToPeak(1) < F.troughToPeak(2) && F.halfWidth(1) < F.halfWidth(2) && F.spectralCentroid(1) > F.spectralCentroid(2) && ...
     F.powerAbove1kHz(1) > F.powerAbove1kHz(2) && F.phasePlotArea(1) > F.phasePlotArea(2);
passed = report("narrow vs broad: shorter times, higher spectral centroid, more power above 1 kHz, larger phase-plot area", ok) && passed;
ok = isequal(F.nPeaks(:)', [2 1 1]) && isequal(F.nTroughs(:)', [1 1 1]) && ~any(F.edgePeakAfter) && all(F.polarity == -1) && ~any(F.nonSomatic);
passed = report("only peaks with a prominence of 20% counted (2; 1 and 1, whose first peak is smaller), one trough each, no edge flags, negative polarity", ok) && passed;
rp = gauss_sum(shapes{2}, t) + 0.02*sin(2*pi*t/0.1).*(t > 0.4);             % small ripples on the broad peak after the trough (seen in templates)
Rp = WaveformFeatures([zeros(1,N_PAD) rp], FS);
passed = report(sprintf("ripples on the peak after the trough are not counted as peaks (nPeaks %d, nTroughs %d)", Rp.nPeaks, Rp.nTroughs), ...
         Rp.nPeaks == 1 && Rp.nTroughs == 1) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Gaussian trough: values known from theory

sigma = 0.12;                                                               % msec
g = [zeros(1,N_PAD) -gauss_sum([1 0 sigma], t)];
G = WaveformFeatures(g, FS);
s_f = 1/(2*sqrt(2)*pi*sigma*1e-3);                                          % power spectrum of a Gaussian: half-normal in f, with this SD (Hz)
theory = [-exp(-0.5)/sigma, exp(-0.5)/sigma, 1, 1/(2*sigma^2), s_f*sqrt(2/pi)/1000];       % slopes, their ratio, area, centroid (kHz)
got = [G.maxDescent, G.maxRecovery, G.descentRecoveryRatio, G.phasePlotArea, G.spectralCentroid];
passed = report(sprintf("Gaussian trough: slopes, their ratio, phase-plot area and spectral centroid within 1%% of theory (max error %.2f%%)", ...
         100*max(abs(got - theory)./abs(theory))), all(abs(got - theory)./abs(theory) < 0.01)) && passed;
theory = [2*sigma*sqrt(2*log(2)), 2*sigma];                                 % half-width and slope trough-to-peak, in msec
got = [G.halfWidth, G.slopeTroughToPeak];
passed = report(sprintf("Gaussian trough: half-width %.4f (theory %.4f) and slope TTP %.4f (theory %.4f) within 6 microsec", got(1), theory(1), ...
         got(2), theory(2)), all(abs(got - theory) < 0.006)) && passed;
theory = erfc(1000/s_f/sqrt(2));
passed = report(sprintf("Gaussian trough: power above 1 kHz %.4f (theory %.4f) within 0.005", G.powerAbove1kHz, theory), ...
         abs(G.powerAbove1kHz - theory) < 0.005) && passed;
r = [zeros(1,N_PAD) gauss_sum([-1 0 0.12; 0.3 2.6 0.5], t)];               % the peak after the trough comes after the end of the template
Rr = WaveformFeatures(r, FS);
ok = Rr.edgePeakAfter && isnan(Rr.troughToPeak) && isnan(Rr.peakAfterRatio) && isnan(Rr.asymmetry) && ...
     ~isnan(Rr.halfWidth) && ~isnan(Rr.maxRecovery) && abs(Rr.peakBeforeRatio) < 0.01;
passed = report("waveform still rising at the end: edge flag set, TTP, P2 and asymmetry empty, the other measures kept", ok) && passed;
q = [zeros(1,N_PAD) gauss_sum([-1 0 0.12; 0.3 0.6 0.2; 0.4 -1.3 0.3], t)];   % the peak before the trough comes before the start
Rq = WaveformFeatures(q, FS);
passed = report(sprintf("waveform falling from the start: P1 edge flag set, P1/|T| = value at the start (%.3f)", Rq.peakBeforeRatio), ...
         Rq.edgePeakBefore && abs(Rq.peakBeforeRatio - Rq.waveform(1)) < 1e-3 && ~isnan(Rq.troughToPeak)) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Amplitude, polarity, unit type, cell input

y = W(1,:);
H = WaveformFeatures([y/0.37; 5*y; -y; y], FS, [NaN NaN NaN 3]);
M = cell2mat(cellfun(@(m) H.(m), H.measures, 'UniformOutput', false));
same = @(a,b) all(abs(a - b) < 1e-9 | (isnan(a) & isnan(b)));
passed = report("same measures for any amplitude", same(M(1,:), M(2,:)) && same(M(1,:), M(4,:))) && passed;
passed = report("positive waveform: flipped (same measures), polarity +1, non-somatic", same(M(1,:), M(3,:)) && H.polarity(3) == 1 && H.nonSomatic(3)) && passed;
passed = report("bombcell type 3: non-somatic; others not", isequal(H.nonSomatic(:)', [false false true true])) && passed;
Hc = WaveformFeatures(num2cell(W, 2), FS, [1 1 2]);
Mc = cell2mat(cellfun(@(m) Hc.(m), Hc.measures, 'UniformOutput', false));
Mm = cell2mat(cellfun(@(m) F.(m), F.measures, 'UniformOutput', false));
passed = report("cell array input = matrix input", isequaln(Mc, Mm)) && passed;
passed = report("scaled waveform stored (trough = -1) with its time axis", abs(min(H.waveform(2,:)) + 1) < 1e-12 && numel(H.time) == N_KEEP && H.time(ALIGN) == 0) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Area of the phase plot = integral of the squared slope / energy

for s=1:numel(shapes)
    te = linspace(t(1), t(end), 30001);
    [v, dv] = gauss_sum(shapes{s}, te);
    scale = -min(v);    v = v/scale;    dv = dv/scale;
    loop = abs(sum(v(1:end-1).*dv(2:end) - v(2:end).*dv(1:end-1)) + v(end)*dv(1) - v(1)*dv(end)) / 2;   % polygon of the closed loop
    area = loop / trapz(te, v.^2);
    passed = report(sprintf("%s: loop area of the phase plot / energy = %.2f, measured %.2f", names{s}, area, F.phasePlotArea(s)), ...
             abs(area - F.phasePlotArea(s))/area < 0.01) && passed;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Edge cases

bad = W(1,:);    bad(5) = 0.01;
try
    WaveformFeatures(bad, FS);
    ok = false;
catch err
    ok = contains(err.message, 'N_PADDING');
end
passed = report("padding that is not zero: error that explains N_PADDING", ok) && passed;
E = WaveformFeatures([zeros(1, N_PAD+N_KEEP); nan(1, N_PAD+N_KEEP); W(1,:)], FS);
Me = cell2mat(cellfun(@(m) E.(m), E.measures, 'UniformOutput', false));
passed = report("template of zeros or NaN: all measures NaN, the other units measured", all(isnan(Me(1:2,:)), 'all') && same(Me(3,:), Mm(1,:))) && passed;
Z = WaveformFeatures({}, FS);
passed = report("no units: empty measures, no error", isempty(Z.troughToPeak) && isempty(Z.polarity)) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. loadKilosortTemplates

folder = tempname;    mkdir(folder);
cleanup = onCleanup(@() rmdir(folder, 's'));
ids = [0 2 5 7 9];      types = [3 1 2 1 0];                                % rows of metadata.csv (new_cluster_id is not the row number)
rows = cell(1, numel(ids));
for k=1:numel(ids)
    rows{k} = single([zeros(1,N_PAD) -gauss_sum([1 0 0.05*k], t)]);
end
fid = fopen(fullfile(folder, 'metadata.csv'), 'w');
fprintf(fid, ',unitType,cluster_id,new_cluster_id,waveform_template,peak_channel\n');
for k=1:numel(ids)
    fprintf(fid, '%d,%d,%d,%d,"[%s]",%d\n', k-1, types(k), ids(k), ids(k), sprintf(' %.8g', rows{k}), 10*k);
end
fclose(fid);
waveforms = rows;
save(fullfile(folder, 'waveforms.mat'), 'waveforms');
st = struct('cluster', {7, 2, 9}, 't_step', 1000/FS);                       % units kept after the QC, in another order
[T1, u1, fs1] = loadKilosortTemplates(folder, st);
ok = isequal(T1, double(vertcat(rows{[4 2 5]}))) && isequal(u1(:)', [1 1 0]) && abs(fs1 - FS) < 1e-6;
passed = report("loader: templates and unit types in the order of spikeTimestamp, sampling rate from t_step", ok) && passed;
delete(fullfile(folder, 'waveforms.mat'));
[T2, u2] = loadKilosortTemplates(folder, st);
passed = report("loader without waveforms.mat: same templates from the text of metadata.csv", max(abs(T2 - T1), [], 'all') < 1e-6 && isequal(u2, u1)) && passed;
try
    loadKilosortTemplates(folder, struct('cluster', {7, 4}, 't_step', 1000/FS));
    ok = false;
catch err
    ok = contains(err.message, 'not in metadata.csv');
end
passed = report("loader: a unit missing from metadata.csv stops with an error", ok) && passed;
[T3, u3] = loadKilosortTemplates(folder, struct('cluster', {}, 't_step', {}));
passed = report("loader: no units -> empty outputs", isempty(T3) && isempty(u3)) && passed;

if passed
    fprintf("\nALL CHECKS PASSED\n");
else
    fprintf("\nSOME CHECKS FAILED\n");
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v, dv] = gauss_sum(shape, t)

% sum of Gaussians (one row of "shape" per Gaussian: amplitude, center, SD) and its exact slope

v = zeros(size(t));     dv = zeros(size(t));
for k=1:size(shape,1)
    g = shape(k,1) * exp(-(t - shape(k,2)).^2 / (2*shape(k,3)^2));
    v = v + g;
    dv = dv - g .* (t - shape(k,2)) / shape(k,3)^2;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = exact_measures(shape, t0, t1)

% the measures of "WaveformFeatures", computed directly on the exact waveform (every 0.1 microsec, exact slope)

te = t0:1e-4:t1;
[v, dv] = gauss_sum(shape, te);
scale = -min(v);    v = v/scale;    dv = dv/scale;
[~, iT] = min(v);
[P1, ~] = max(v(1:iT));
[P2, k] = max(v(iT:end));    iP2 = iT + k - 1;
R.troughToPeak = te(iP2) - te(iT);
R.peakAfterRatio = P2;
R.peakBeforeRatio = P1;
R.asymmetry = (max(P2,0) - max(P1,0)) / (max(P2,0) + max(P1,0));
below = find(v <= -0.5);
R.halfWidth = te(below(end)) - te(below(1));                                % (one trough: the points below half depth are contiguous)
kD = find(te >= te(iT) - 0.5, 1):iT;
[R.maxDescent, k] = min(dv(kD));    iD = kD(k);
[R.maxRecovery, k] = max(dv(iT:iP2));    iR = iT + k - 1;
R.descentRecoveryRatio = abs(R.maxDescent) / R.maxRecovery;
R.slopeTroughToPeak = te(iR) - te(iD);
R.phasePlotArea = trapz(te, dv.^2) / trapz(te, v.^2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = report(name, ok)

if ok
    fprintf("  PASS  %s\n", name);
else
    fprintf("  FAIL  %s\n", name);
end

end
