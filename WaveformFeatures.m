function WaveformFeature = WaveformFeatures(waveforms, samplingRate, unitType)

% function to measure the shape of the spike waveform of each unit (new, November 2026)
% It works on the Kilosort 2.5 template of each unit on its peak electrode (file "waveforms.mat" of the recording folder; see
% "loadKilosortTemplates"). The templates are whitened and have unit norm, so their amplitude is arbitrary: only the shape is measured
% (times, ratios, and slopes of the waveform scaled to trough = -1). The units are compared with each other; thresholds from the literature,
% set on raw waveforms (e.g., "narrow if trough-to-peak <= 0.4 ms"), do not apply to these templates
% Choice of the measures (literature review of November 2026): trough-to-peak, half-width, peak ratios and asymmetry are the standard
% measures (Barthó et al. 2004; Jia et al. 2019; CellExplorer, Petersen et al. 2021; SpikeInterface template metrics; DeePhys, Hornauer et
% al. 2024). The slopes and the area of the phase plot (waveform vs its slope) adapt the intracellular phase plot: the extracellular spike
% is close to the negative slope of the membrane potential (Bean 2007; Henze et al. 2000), so the trough marks the fastest rise of the
% action potential and P2 its fastest fall. The area of the phase plot divided by the energy of the waveform equals (2*pi*RMS frequency)^2,
% so it also measures the bandwidth of the spectrum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   waveforms:      template of each unit on its peak electrode, as written by Kilosort 2.5 (82 samples, the first N_PADDING are zeros):
%                   cell array of vectors (as in "waveforms.mat") or matrix with one row per unit
%   samplingRate:   sampling rate of the recording, in Hz
%   unitType:       (optional) type of each unit given by bombcell (=0, noise; =1, good; =2, MUA; =3, non-somatic); [] or omitted if not known

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Landmarks of the waveform scaled to trough = -1: T (trough, the lowest point), P1 (highest point before T), P2 (highest point after T)
% The 12 measures (one value per unit, NaN if it cannot be measured):
% WaveformFeature.troughToPeak(i)          time from T to P2, in msec
% WaveformFeature.halfWidth(i)             width of the trough at half its depth (-0.5), in msec
% WaveformFeature.peakAfterRatio(i)        P2 / |T|
% WaveformFeature.peakBeforeRatio(i)       P1 / |T|
% WaveformFeature.asymmetry(i)             (P2 - P1) / (P2 + P1), with a peak below zero counted as 0 (CellExplorer "AB ratio")
% WaveformFeature.maxDescent(i)            steepest slope on the way into the trough (most negative slope in the DESCENT_WINDOW msec before T),
%                                          in 1/msec (units of the scaled waveform per msec)
% WaveformFeature.maxRecovery(i)           steepest slope after the trough (most positive slope between T and P2), in 1/msec
% WaveformFeature.descentRecoveryRatio(i)  |maxDescent| / maxRecovery
% WaveformFeature.slopeTroughToPeak(i)     time from maxDescent to maxRecovery, in msec (CellExplorer "derivative trough-to-peak")
% WaveformFeature.phasePlotArea(i)         area inside the phase plot (waveform vs its slope) divided by the energy of the waveform:
%                                          integral(slope^2) / integral(waveform^2) over the whole template, in 1/msec^2
% WaveformFeature.spectralCentroid(i)      power-weighted mean frequency of the waveform, in kHz
% WaveformFeature.powerAbove1kHz(i)        fraction of the power of the waveform above SPLIT_FREQUENCY (1 kHz)
% Information about each unit:
% WaveformFeature.polarity(i)              =-1, the largest deflection is negative (trough); =+1, it is positive: the waveform is flipped before
%                                          measuring (as CellExplorer does), so the measures describe its main (positive) peak
% WaveformFeature.unitType(i)              type given by bombcell (NaN if not provided)
% WaveformFeature.nonSomatic(i)            =true for the non-somatic units of bombcell (type 3) and the units with a positive largest deflection
%                                          (typical of axons); they are measured, but kept out of the waveform clustering
% WaveformFeature.nPeaks(i)                number of peaks with a prominence of at least PEAK_FRACTION of the trough depth (quality check, as
%                                          bombcell: it calls "noise" the units with more than 2 peaks or more than 1 trough). The prominence is
%                                          the height of a peak above the higher of the two lowest points that separate it from a higher point
%                                          (or from the end of the template) on each side, as in "findpeaks"; small ripples of the template do
%                                          not count
% WaveformFeature.nTroughs(i)              number of troughs with the same prominence (quality check)
% WaveformFeature.edgePeakBefore(i)        =true if P1 is the first sample of the template: there is no peak before the trough within the window,
%                                          and peakBeforeRatio is the value at the start of the template (usually close to 0)
% WaveformFeature.edgePeakAfter(i)         =true if P2 is the last sample of the template: the real peak comes later, so troughToPeak,
%                                          peakAfterRatio and asymmetry are NaN
% WaveformFeature.waveform(i,:)            template without the padding, scaled to trough = -1 (flipped if polarity = +1)
% WaveformFeature.time                     time of each sample of "waveform", in msec, from the sample ALIGN_SAMPLE (where Kilosort 2.5
%                                          aligns the trough, +/- 1 sample)
% WaveformFeature.measures                 names of the 12 measures (fields above), in this order
% WaveformFeature.parameters               parameters used

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_PADDING = 21;                     % number of zeros at the start of each template, added by the export of Kilosort 2.5; they are removed before
                                    % measuring (=0 for templates without padding)
ALIGN_SAMPLE = 21;                  % sample (after removing the padding) where Kilosort 2.5 aligns the trough; only used for the time axis
UPSAMPLING = 10;                    % factor of the cubic-spline interpolation before measuring (x10: 5 microsec steps at 19.75 kHz); without
                                    % it, the slopes are underestimated by 10-30% and the times jump in steps of one sample (50 microsec)
DESCENT_WINDOW = 0.5;               % window before the trough where the steepest descent is searched, in msec (as CellExplorer)
PEAK_FRACTION = 0.2;                % peaks / troughs counted by nPeaks / nTroughs: prominence of at least this fraction of the trough depth
                                    % (as bombcell, "minThreshDetectPeaksTroughs")
SPLIT_FREQUENCY = 1000;             % frequency that splits the spectrum for "powerAbove1kHz", in Hz
N_FFT = 1024;                       % number of points of the Fourier transform (the template is padded with zeros: smoother spectrum,
                                    % same information)
apply_parameters('WaveformFeatures');                                       % values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

MEASURES = {'troughToPeak','halfWidth','peakAfterRatio','peakBeforeRatio','asymmetry','maxDescent','maxRecovery', ...
            'descentRecoveryRatio','slopeTroughToPeak','phasePlotArea','spectralCentroid','powerAbove1kHz'};
K_MSEC_SEC = 1/1000;                % factor to rescale from msec to sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check the inputs and remove the padding

if iscell(waveforms) && isempty(waveforms)                                  % no units (e.g., a recording without units after the QC)
    W = zeros(0, N_PADDING + 1);
elseif iscell(waveforms)                                                    % one vector per unit (as in "waveforms.mat")
    lengths = cellfun(@numel, waveforms);
    if any(lengths ~= lengths(1))
        error('WaveformFeatures: the templates do not have the same number of samples');
    end
    W = zeros(numel(waveforms), lengths(1));
    for i=1:numel(waveforms)
        W(i,:) = double(waveforms{i}(:))';
    end
else
    W = double(waveforms);
end
n_units = size(W,1);
if nargin < 3 || isempty(unitType)
    unitType = nan(n_units,1);
end
unitType = double(unitType(:));
if numel(unitType) ~= n_units
    error('WaveformFeatures: %g unit types for %g templates', numel(unitType), n_units);
end
if n_units == 0 && size(W,2) <= N_PADDING
    W = zeros(0, N_PADDING + 1);
end
if N_PADDING > 0
    pad = W(:, 1:min(N_PADDING, size(W,2)));     pad(isnan(pad)) = 0;      % (a unit without template, all NaN, is not checked)
    if size(W,2) <= N_PADDING || any(pad ~= 0, 'all')
        error(['WaveformFeatures: the first %g samples of the templates are not all zeros (padding of Kilosort 2.5). ' ...
               'For templates without padding, set N_PADDING = 0'], N_PADDING);
    end
    W = W(:, N_PADDING+1:end);
end
n_samples = size(W,2);
dt = 1/samplingRate/K_MSEC_SEC;                                             % time step, in msec
t = (0:n_samples-1)*dt;                                                     % time of each sample, in msec
tu = linspace(t(1), t(end), (n_samples-1)*UPSAMPLING + 1);                  % time of the interpolated samples, in msec
f = (0:N_FFT/2)*samplingRate/N_FFT;                                         % frequencies of the one-sided spectrum, in Hz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Measure each unit

for m=1:numel(MEASURES)
    WaveformFeature.(MEASURES{m}) = nan(n_units,1);
end
WaveformFeature.polarity = nan(n_units,1);
WaveformFeature.unitType = unitType;
WaveformFeature.nonSomatic = unitType == 3;
WaveformFeature.nPeaks = nan(n_units,1);
WaveformFeature.nTroughs = nan(n_units,1);
WaveformFeature.edgePeakBefore = false(n_units,1);
WaveformFeature.edgePeakAfter = false(n_units,1);
WaveformFeature.waveform = nan(n_units, n_samples);

for i=1:n_units
    y = W(i,:);
    if any(~isfinite(y)) || all(y == 0)                                    % no template: every measure stays NaN
        continue
    end

    % polarity and scaling (trough = -1)
    [y_min, ~] = min(y);    [y_max, ~] = max(y);
    if y_max > -y_min                                                       % the largest deflection is positive: flip the waveform
        WaveformFeature.polarity(i) = 1;
        WaveformFeature.nonSomatic(i) = true;
        y = -y;
    else
        WaveformFeature.polarity(i) = -1;
    end
    y = y / (-min(y));                                                      % trough = -1 (the baseline of the template is a true zero)
    WaveformFeature.waveform(i,:) = y;

    % interpolated waveform and its slope
    yu = interp1(t, y, tu, 'spline');                                       % cubic spline (not-a-knot), x UPSAMPLING
    yu = yu / (-min(yu));                                                   % the interpolated trough can be slightly below the sampled one
    dyu = gradient(yu, tu);                                                 % slope, in 1/msec
    n_u = numel(yu);

    % landmarks: T, P1, P2
    [~, iT] = min(yu);
    [P1, iP1] = max(yu(1:iT));
    [P2, k] = max(yu(iT:end));    iP2 = iT + k - 1;
    WaveformFeature.edgePeakBefore(i) = iP1 == 1;
    WaveformFeature.edgePeakAfter(i) = iP2 == n_u;
    WaveformFeature.peakBeforeRatio(i) = P1;                                % |T| = 1
    if ~WaveformFeature.edgePeakAfter(i)
        WaveformFeature.troughToPeak(i) = tu(iP2) - tu(iT);
        WaveformFeature.peakAfterRatio(i) = P2;
        a = max(P1, 0);    b = max(P2, 0);                                  % a peak below zero counts as no peak
        if a + b > 0
            WaveformFeature.asymmetry(i) = (b - a) / (b + a);
        end
    end

    % half-width: crossings of -0.5 on each side of the trough (linear interpolation between interpolated samples)
    kL = find(yu(1:iT) > -0.5, 1, 'last');                                  % last point above half depth before the trough
    kR = iT - 1 + find(yu(iT:end) > -0.5, 1, 'first');                     % first point above half depth after the trough
    if ~isempty(kL) && ~isempty(kR)
        tL = tu(kL) + (-0.5 - yu(kL)) * (tu(kL+1) - tu(kL)) / (yu(kL+1) - yu(kL));
        tR = tu(kR-1) + (-0.5 - yu(kR-1)) * (tu(kR) - tu(kR-1)) / (yu(kR) - yu(kR-1));
        WaveformFeature.halfWidth(i) = tR - tL;
    end

    % slopes: steepest descent before the trough, steepest recovery between the trough and P2 (or the end of the template)
    kD = find(tu >= tu(iT) - DESCENT_WINDOW, 1, 'first'):iT;
    [descent, k] = min(dyu(kD));    iD = kD(k);
    [recovery, k] = max(dyu(iT:iP2));    iR = iT + k - 1;
    WaveformFeature.maxDescent(i) = descent;
    WaveformFeature.maxRecovery(i) = recovery;
    if recovery > 0
        WaveformFeature.descentRecoveryRatio(i) = abs(descent) / recovery;
    end
    WaveformFeature.slopeTroughToPeak(i) = tu(iR) - tu(iD);

    % area of the phase plot (= integral of the squared slope), divided by the energy of the waveform
    WaveformFeature.phasePlotArea(i) = trapz(tu, dyu.^2) / trapz(tu, yu.^2);

    % spectrum of the template (sampled waveform; the zeros before and after it are its true baseline)
    X = fft(y, N_FFT);
    P = abs(X(1:N_FFT/2+1)).^2;
    P(2:end-1) = 2*P(2:end-1);                                              % one-sided power spectrum
    WaveformFeature.spectralCentroid(i) = sum(f .* P) / sum(P) * K_MSEC_SEC;   % in kHz
    WaveformFeature.powerAbove1kHz(i) = sum(P(f > SPLIT_FREQUENCY)) / sum(P);

    % number of peaks and troughs with enough prominence (interpolated waveform; the trough depth is 1)
    WaveformFeature.nPeaks(i) = count_prominent(yu, PEAK_FRACTION);
    WaveformFeature.nTroughs(i) = count_prominent(-yu, PEAK_FRACTION);
end

WaveformFeature.time = ((1:n_samples) - ALIGN_SAMPLE) * dt;                % in msec
WaveformFeature.measures = MEASURES;
WaveformFeature.parameters = struct('N_PADDING', N_PADDING, 'ALIGN_SAMPLE', ALIGN_SAMPLE, 'UPSAMPLING', UPSAMPLING, ...
    'DESCENT_WINDOW', DESCENT_WINDOW, 'PEAK_FRACTION', PEAK_FRACTION, 'SPLIT_FREQUENCY', SPLIT_FREQUENCY, 'N_FFT', N_FFT, ...
    'samplingRate', samplingRate);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = count_prominent(y, min_prominence)

% number of local maxima of y with a prominence of at least min_prominence (same definition as "findpeaks", without its toolbox):
% the prominence of a peak is its height above the higher of the two lowest points between the peak and the nearest higher point on
% each side (or the end of y)

k = find(y(2:end-1) > y(1:end-2) & y(2:end-1) >= y(3:end)) + 1;            % local maxima (the first point of a flat top)
n = 0;
for j=k(:)'
    left = find(y(1:j-1) > y(j), 1, 'last');                                % nearest higher point on the left
    if isempty(left), left = 1; end
    right = j + find(y(j+1:end) > y(j), 1, 'first');                        % nearest higher point on the right
    if isempty(right), right = numel(y); end
    if y(j) - max(min(y(left:j)), min(y(j:right))) >= min_prominence
        n = n + 1;
    end
end

end
