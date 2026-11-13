function passed = test_firing_pattern_features(code_dir)

% Function-level test of "FiringPatternFeatures" on synthetic spike trains (new, November 2026)
% Like "test_waveform_features", it does not compare two versions of the code: it checks that
%   1. exact values on deterministic trains: LV and CV of alternating intervals, a periodic train, and the burst index of spike pairs
%   2. Poisson trains give LV, CV ~ 1 and burst index ~ 0; trains of short bursts give burst index > 0 and LV > 1; gamma trains are regular
%   3. the spikes inside the network bursts (widened by GUARD) are left out, no interval or pair of the autocorrelogram spans a burst,
%      the bursts given as the structure NetworkBurstDetected or as a matrix give the same result, and a NetworkBurstDetected with
%      nBurst = NaN (candidates rejected by the final check) or 0 means no network bursts (revised November 2026)
%   4. relativeRate = log10 of the rate outside the bursts minus the median of the reference units; NaN with too few reference units
%   5. MIN_SPIKES: the features are NaN below it; edge cases: no units, no bursts, bursts that cover the whole window
% No reference recording is needed
%
% Example (from a terminal):
%   matlab -batch "addpath('<repo>/tests'); test_firing_pattern_features('<repo>')"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dir:             folder with the code to be tested

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% passed:               =true if all the checks passed (each check is printed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(char(code_dir));
passed = true;
rng(20261111);
T0 = 300;   T1 = 900;                                                       % analysis window, in sec (as the baselines)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Exact values on deterministic trains

isi = repmat([10 20], 1, 100);                                              % alternating 10 and 20 msec
t = T0*1000 + 1000 + cumsum(isi)';
F = FiringPatternFeatures(train({t}), [T0 T1], []);
lv_theory = 3 * ((10-20)/30)^2;                                             % = 1/3 for every pair
cv_theory = std(isi(2:end)) / mean(isi(2:end));                             % 200 spikes: the 199 intervals after the first spike
passed = report(sprintf("alternating 10/20 msec: LV %.6f (theory 1/3), CV %.6f (theory %.6f)", F.lv, F.cv, cv_theory), ...
    abs(F.lv - lv_theory) < 1e-12 && abs(F.cv - cv_theory) < 1e-12) && passed;

t = (T0*1000 + 500 + (0:199)*100)';                                         % periodic, 10 Hz
F = FiringPatternFeatures(train({t}), [T0 T1], []);
passed = report(sprintf("periodic train: LV %.2g and CV %.2g are 0", F.lv, F.cv), F.lv < 1e-12 && F.cv < 1e-12) && passed;

n_pairs = 100;                                                              % pairs of spikes 4 msec apart, one pair per sec: the
starts = T0*1000 + 1000*(1:n_pairs)';                                       % autocorrelogram has n_pairs counts in the 4-5 msec bin and
t = sort([starts; starts + 4]);                                             % none at 200-300 msec
F = FiringPatternFeatures(train({t}), [T0 T1], []);
bi_theory = log10((n_pairs/2 + 0.5) / 0.5);
passed = report(sprintf("pairs 4 msec apart: burst index %.4f (theory %.4f)", F.burstIndex, bi_theory), abs(F.burstIndex - bi_theory) < 1e-12) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Statistical trains

poisson = sort(T0*1000 + (T1-T0)*1000*rand(3000, 1));                       % 5 Hz Poisson (used in section 3 too)
poisson20 = sort(T0*1000 + (T1-T0)*1000*rand(12000, 1));                    % 20 Hz Poisson: ~240 counts per msec of the autocorrelogram
gamma_isi = gamrnd(20, 200/20, 2900, 1);                                    % gamma, shape 20, mean 200 msec (5 Hz): regular
gamma_t = T0*1000 + cumsum(gamma_isi);   gamma_t = gamma_t(gamma_t < T1*1000);
b_start = sort(T0*1000 + (T1-T0)*1000*rand(800, 1));                        % bursts of 3 spikes 3-5 msec apart at random times
bursty = sort([b_start; b_start + 3 + 2*rand(800,1); b_start + 7 + 2*rand(800,1)]);
F = FiringPatternFeatures(train({poisson20, gamma_t, bursty}), [T0 T1], []);
passed = report(sprintf("Poisson (20 Hz): LV %.3f, CV %.3f within 0.08 of 1, burst index %.3f within 0.1 of 0", F.lv(1), F.cv(1), F.burstIndex(1)), ...
    abs(F.lv(1) - 1) < 0.08 && abs(F.cv(1) - 1) < 0.08 && abs(F.burstIndex(1)) < 0.1) && passed;
passed = report(sprintf("gamma (shape 20): LV %.3f and CV %.3f below 0.3", F.lv(2), F.cv(2)), F.lv(2) < 0.3 && F.cv(2) < 0.3) && passed;
passed = report(sprintf("bursts of 3 spikes: burst index %.2f > 1, LV %.2f > 1.3", F.burstIndex(3), F.lv(3)), F.burstIndex(3) > 1 && F.lv(3) > 1.3) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Network bursts left out

bursts = [400 401.5; 550.2 552; 700 702.5];                                 % [start last], in sec
guard = 0.5;
inside = sort(cell2mat(arrayfun(@(j) bursts(j,1)*1000 + (bursts(j,2)-bursts(j,1))*1000*rand(300,1), (1:3)', 'UniformOutput', false)));
t = sort([poisson; inside]);
D.IndivBurst = struct('start', num2cell(bursts(:,1))', 'last', num2cell(bursts(:,2))');
Fs = FiringPatternFeatures(train({t}), [T0 T1], D);
Fm = FiringPatternFeatures(train({t}), [T0 T1], bursts);
keep = true(size(poisson));                                                 % the Poisson spikes outside the widened bursts
for j=1:3
    keep = keep & ~(poisson >= (bursts(j,1)-guard)*1000 & poisson < (bursts(j,2)+guard)*1000);
end
periods = [T0 bursts(1,1)-guard; bursts(1,2)+guard bursts(2,1)-guard; bursts(2,2)+guard bursts(3,1)-guard; bursts(3,2)+guard T1]*1000;
[lv_ref, cv_ref, bi_ref] = by_period(poisson, periods);
time_out = (T1-T0) - sum(bursts(:,2) - bursts(:,1) + 2*guard);
passed = report(sprintf("spikes in the widened bursts left out: %d of %d outside (expected %d), time outside %.1f sec (expected %.1f)", ...
    Fs.nSpikesOutside, numel(t), sum(keep), Fs.timeOutside, time_out), Fs.nSpikesOutside == sum(keep) && Fs.nSpikes == numel(t) && ...
    abs(Fs.timeOutside - time_out) < 1e-9 && abs(Fs.rateOutside - sum(keep)/time_out) < 1e-12) && passed;
passed = report(sprintf("no interval or pair spans a burst: LV %.6f, CV %.6f, burst index %.6f equal the values computed period by period", ...
    Fs.lv, Fs.cv, Fs.burstIndex), abs(Fs.lv - lv_ref) < 1e-12 && abs(Fs.cv - cv_ref) < 1e-12 && abs(Fs.burstIndex - bi_ref) < 1e-12) && passed;
passed = report("bursts as NetworkBurstDetected = bursts as a matrix", isequaln(Fs, Fm)) && passed;
D.nBurst = 3;
passed = report("NetworkBurstDetected with nBurst = 3: its bursts are left out", isequaln(FiringPatternFeatures(train({t}), [T0 T1], D), Fm)) && passed;
Fw = FiringPatternFeatures(train({t}), [T0 T1], []);
D.nBurst = NaN;
Fn = FiringPatternFeatures(train({t}), [T0 T1], D);
D.nBurst = 0;
Fz = FiringPatternFeatures(train({t}), [T0 T1], D);
passed = report(sprintf("nBurst = NaN (candidates rejected by the final check) or 0: no network bursts, the whole window is used (%.0f sec)", ...
    Fn.timeOutside), isequaln(Fn, Fw) && isequaln(Fz, Fw) && Fn.timeOutside == T1-T0) && passed;
Fu = FiringPatternFeatures(train({t}), [T0 T1], flipud(bursts));
passed = report("the order of the bursts does not matter", isequaln(Fu, Fm)) && passed;
F0 = FiringPatternFeatures(train({t}), [T0 T1], []);
passed = report(sprintf("without bursts, the whole window is used (%d spikes, %.0f sec)", F0.nSpikesOutside, F0.timeOutside), ...
    F0.nSpikesOutside == numel(t) && F0.timeOutside == T1-T0 && isequal(F0.periods, [T0 T1]*1000)) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Relative rate

rates = [0.5 1 2 4 8 16];
tr = arrayfun(@(r) sort(T0*1000 + (T1-T0)*1000*rand(round(r*(T1-T0)), 1)), rates, 'UniformOutput', false);
F = FiringPatternFeatures(train(tr), [T0 T1], []);
expected = log10(F.rateOutside) - median(log10(F.rateOutside));
passed = report(sprintf("relativeRate = log10(rate) - median over the units (%s)", mat2str(F.relativeRate', 3)), ...
    max(abs(F.relativeRate - expected)) < 1e-12) && passed;
ref = logical([1 1 1 1 1 0]);
F = FiringPatternFeatures(train(tr), [T0 T1], [], ref);
expected = log10(F.rateOutside) - median(log10(F.rateOutside(ref)));
passed = report("reference units: the median of the reference units only", max(abs(F.relativeRate - expected)) < 1e-12) && passed;
F = FiringPatternFeatures(train(tr), [T0 T1], [], logical([1 1 1 1 0 0]));
passed = report("fewer than 5 reference units: relativeRate NaN, the other features measured", all(isnan(F.relativeRate)) && ...
    all(~isnan(F.lv))) && passed;
ok = false;
try
    FiringPatternFeatures(train(tr), [T0 T1], [], true(1, 3));
catch
    ok = true;
end
passed = report("'reference' of the wrong size: error", ok) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Minimum number of spikes and edge cases

t49 = sort(T0*1000 + (T1-T0)*1000*rand(49, 1));     t50 = sort(T0*1000 + (T1-T0)*1000*rand(50, 1));
F = FiringPatternFeatures(train([{t49, t50}, tr(1:4)]), [T0 T1], []);
passed = report("49 spikes outside: the 4 features NaN, enoughSpikes false; 50 spikes: measured", ...
    all(isnan([F.relativeRate(1) F.lv(1) F.cv(1) F.burstIndex(1)])) && ~F.enoughSpikes(1) && ...
    all(~isnan([F.relativeRate(2) F.lv(2) F.cv(2) F.burstIndex(2)])) && F.enoughSpikes(2)) && passed;
F = FiringPatternFeatures(train({}), [T0 T1], []);
passed = report("no units: empty outputs", isempty(F.lv) && isempty(F.relativeRate) && isempty(F.nSpikes)) && passed;
F = FiringPatternFeatures(train({poisson}), [T0 T1], [T0 T1]);
passed = report("a burst that covers the whole window: no time outside, rate NaN, features NaN", F.timeOutside == 0 && ...
    isnan(F.rateOutside) && isnan(F.lv) && ~F.enoughSpikes) && passed;
D0.IndivBurst = struct();
F = FiringPatternFeatures(train({poisson}), [T0 T1], D0);
passed = report("NetworkBurstDetected without bursts: the whole window is used", F.timeOutside == T1-T0) && passed;
t = [poisson; (T0 - 10)*1000; (T1 + 10)*1000];
F = FiringPatternFeatures(train({t}), [T0 T1], []);
passed = report("spikes outside the analysis window are ignored", F.nSpikes == numel(poisson)) && passed;

if passed
    fprintf("test_firing_pattern_features: all checks passed\n");
else
    fprintf("test_firing_pattern_features: SOME CHECKS FAILED\n");
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function st = train(times)
% structure array as PerCondition(cond).spikeTrain (only the field used)
st = struct('timeStamp', times);
if isempty(times)
    st = struct('timeStamp', {});
end
end

function [lv, cv, bi] = by_period(t, periods)
% LV, CV and burst index computed separately in each period, then pooled (independent of the code tested)
I1 = [];  I2 = [];  isi = [];  short = 0;  long = 0;
for p=1:size(periods, 1)
    tp = t(t >= periods(p,1) & t < periods(p,2));
    d = diff(tp);
    isi = [isi; d];                                                         %#ok<AGROW>
    I1 = [I1; d(1:end-1)];  I2 = [I2; d(2:end)];                            %#ok<AGROW>
    lags = tp' - tp;                                                        % all pairs (later minus earlier)
    lags = lags(lags > 0);
    short = short + sum(lags >= 3 & lags < 5);
    long = long + sum(lags >= 200 & lags <= 300);
end
lv = 3 * mean(((I1 - I2) ./ (I1 + I2)).^2);
cv = std(isi) / mean(isi);
bi = log10((short/2 + 0.5) / (long/100 + 0.5));
end

function ok = report(name, ok)
if ok
    fprintf("  PASS  %s\n", name);
else
    fprintf("  FAIL  %s\n", name);
end
end
