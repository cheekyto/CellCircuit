function passed = test_cNECompare(code_dir)

% Function-level test of "cNECompare" and of the tracking of "cNEDetection" (TEMPLATE) on synthetic spike trains (November 2026). It checks
%   1. tracking the cNEs of a window on that same window gives the activations of the detection exactly; tracking uses the activation
%      settings of the template; a template with other units is refused
%   2. a "drug" that changes the ensembles (reference 0-600 s, test 600-1200 s): an ensemble kept (found again, same activation rate and
%      member rates), an ensemble whose members keep their rates but no longer fire together (activations fall, member rates do not; lost
%      by the separate detection), an ensemble whose members are silenced (silent members, no activations), a new ensemble (new cNE);
%      chance level: the kept ensemble far above it, the uncoordinated one at it in the test window (excess ratio near 0)
%   3. no change (same process in both windows, with large network events where the ensembles fire in a fixed order): every cNE found
%      again, same rates, same roles, the large events of the test window in the order of the reference (template match 1, agreement
%      above 0.85, p < 0.01); with the order reversed in the test window: template match -1, agreement below 0.15
%   4. a silenced test window (no spikes), no cNEs in the reference, reuse of a given detection, units that differ; the figure is drawn
%
% Example (from a terminal):
%   matlab -batch "addpath('<repo>/tests'); test_cNECompare('<repo>')"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dir:             folder with the code to be tested

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% passed:               =true if all the checks passed (each check is printed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(char(code_dir));
passed = true;
N = 40;     T_SEC = 1200;   REF = [0 600];  TEST = [600 1200];
G = {1:6, 10:15, 25:31, 33:38};                                             % planted ensembles (G{4} only in the test window)
FAST = struct('ANALYSIS', struct('N_PERM', 500));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. A "drug" that changes the ensembles

rng(1);
trains = background(N, T_SEC);
trains = plant(trains, G{1}, events(0.3, [0 T_SEC]), 0.6);                  % kept
trains = plant(trains, G{2}, events(0.3, REF), 0.6);                        % coordinated in the reference only ...
for i = G{2}                                                                % ... same rate in the test window, independent spikes
    n_extra = round(0.6*0.3*diff(TEST));
    trains{i} = [trains{i}; TEST(1)*1000 + rand(n_extra,1)*diff(TEST)*1000];
end
trains = plant(trains, G{3}, events(0.3, REF), 0.6);                        % silenced in the test window
for i = G{3}
    trains{i} = trains{i}(trains{i} < TEST(1)*1000);
end
trains = plant(trains, G{4}, events(0.3, TEST), 0.6);                       % new in the test window
trains = cellfun(@sort, trains, 'UniformOutput', false);

tic;
[C, fig] = cNECompare(struct('spikeTimestamp', {trains}, 'tRange', REF, 'name', "baseline"), ...
    struct('spikeTimestamp', {trains}, 'tRange', TEST, 'name', "drug"), FAST);
close(fig);
fprintf("  (comparison: %.0f sec)\n", toc);
T = C.cNEs;     S = C.summary;
k = zeros(1,3);
for g = 1:3
    hit = find(cellfun(@(m) isequal(sort(m(:))', G{g}), C.reference.cNE.members), 1);
    if ~isempty(hit)
        k(g) = hit;
    end
end
passed = report(sprintf("reference: the 3 planted ensembles detected (%d cNEs)", S.nCNE), S.nCNE == 3 && all(k > 0)) && passed;
if all(k > 0)
    passed = report(sprintf("kept ensemble: found again (similarity %.2f), activation ratio %.2f, member rate ratio %.2f", ...
        T.similarity(k(1)), T.rateRatio(k(1)), T.memberRateRatio(k(1))), T.match(k(1)) > 0 && T.similarity(k(1)) > 0.9 && ...
        abs(log(T.rateRatio(k(1)))) < log(1.4) && abs(log(T.memberRateRatio(k(1)))) < log(1.2)) && passed;
    passed = report(sprintf("uncoordinated ensemble: activation ratio %.2f, member rate ratio %.2f, lost", T.rateRatio(k(2)), ...
        T.memberRateRatio(k(2))), T.rateRatio(k(2)) < 0.2 && abs(log(T.memberRateRatio(k(2)))) < log(1.2) && T.match(k(2)) == 0) && passed;
    passed = report(sprintf("silenced ensemble: %d silent members, no activations, member rate 0, lost", T.silentMembers(k(3))), ...
        T.silentMembers(k(3)) == numel(G{3}) && T.rateTest(k(3)) == 0 && T.memberRateTest(k(3)) == 0 && T.match(k(3)) == 0) && passed;
    passed = report(sprintf(['chance level: kept ensemble %.0f%% / %.0f%% of its activations (excess ratio %.2f); uncoordinated ensemble in the ' ...
        'test window %.0f%% (excess ratio %.2f); %d reference cNEs at chance'], T.chanceShareReference(k(1)), T.chanceShareTest(k(1)), ...
        T.excessRatio(k(1)), T.chanceShareTest(k(2)), T.excessRatio(k(2)), S.nAtChance), all([T.chanceShareReference(k(1)) ...
        T.chanceShareTest(k(1)) T.chanceShareReference(k(2))] < 25) && abs(log(T.excessRatio(k(1)))) < log(1.4) && ...
        T.chanceShareTest(k(2)) > 40 && T.chanceShareTest(k(2)) < 250 && T.excessRatio(k(2)) < 0.1 && S.nAtChance == 0) && passed;
end
new_members = cellfun(@(m) sort(m(:))', C.detected.cNE.members(C.newCNEs), 'UniformOutput', false);
passed = report(sprintf("summary: %d found, %d lost, %d new (members %s)", S.nFound, S.nLost, S.nNew, strjoin(cellfun(@mat2str, ...
    new_members, 'UniformOutput', false), ", ")), S.nFound == 1 && S.nLost == 2 && S.nNew == 1 && isequal(new_members{1}, G{4}) && ...
    abs(S.foundShare - 100/3) < 1e-9) && passed;
passed = report("no network events: roles and order not compared (NaN)", isnan(S.rolesKept) && isnan(S.referenceAgreement) && ...
    S.noBurstingReference && S.noBurstingTest) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Tracking on the same window = detection

C0 = C.reference.cNE;
C1 = cNEDetection(trains, REF, struct('TEMPLATE', C0, 'BIN_SIZE', 10));
same = isequal(C1.events, C0.events) && isequal(C1.activityThreshold, C0.activityThreshold) && max(abs(C1.activity(:) - C0.activity(:))) ...
    < 1e-9 && isequal(C1.members, C0.members) && isequal(C1.weights, C0.weights) && isequal(C1.binTime, C0.binTime);
passed = report("tracking on the window of the detection = its activations (bin size of the template kept)", same) && passed;
try
    cNEDetection(trains(1:end-1), TEST, struct('TEMPLATE', C0));
    refused = false;
catch
    refused = true;
end
passed = report("tracking refuses a template with other units", refused) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. No change, with large network events (the ensembles fire in the order 1 -> 2 -> 3); then the order reversed in the test window

rng(2);
s = (5:10:T_SEC - 5)';                                                      % large events of 1 s every 10 s, tails of 1 s
NB = network(s);
trains = background(N, T_SEC);
for g = 1:3
    trains = plant(trains, G{g}, events(0.05, [0 T_SEC]), 0.6);             % a few events outside the network events
    trains = plant(trains, G{g}, s + 0.05 + 0.2*(g - 1), 0.8);              % in each large event, in order
end
trains = cellfun(@sort, trains, 'UniformOutput', false);
C = cNECompare(struct('spikeTimestamp', {trains}, 'tRange', REF, 'NetworkBurstDetected', NB), ...
    struct('spikeTimestamp', {trains}, 'tRange', TEST, 'NetworkBurstDetected', NB), FAST);
S = C.summary;  T = C.cNEs;
passed = report(sprintf("no change: %d of %d cNEs found again, %d new; activation ratios %s; member rate ratios %s", S.nFound, S.nCNE, ...
    S.nNew, mat2str(T.rateRatio', 2), mat2str(T.memberRateRatio', 2)), S.nCNE == 3 && S.nFound == 3 && S.nNew == 0 && ...
    all(abs(log(T.rateRatio)) < log(1.3)) && all(abs(log(T.memberRateRatio)) < log(1.1))) && passed;
passed = report(sprintf("no change: roles kept %d, lost %d, gained %d (%s)", S.rolesKept, S.rolesLost, S.rolesGained, ...
    strjoin(T.rolesReference, " | ")), S.rolesKept >= 3 && S.rolesLost == 0 && S.rolesGained == 0) && passed;
passed = report(sprintf("no change: large events in the reference order, agreement %.2f (p %.4f), template match %.2f", ...
    S.referenceAgreement, S.referenceAgreementP, S.referenceMatch), S.referenceAgreement > 0.85 && S.referenceAgreementP < 0.01 && ...
    S.referenceMatch == 1) && passed;

for g = 1:3                                                                 % order reversed in the test window (3 -> 2 -> 1)
    for i = G{g}
        late = trains{i} >= TEST(1)*1000;
        t = trains{i}(late)/1000;
        q = discretize(t, reshape([s, s + 1]', [], 1));                     % spikes inside a large event (odd intervals)
        inside = ~isnan(q) & mod(q, 2) == 1;
        t(inside) = t(inside) + 0.2*(3 - g) - 0.2*(g - 1);                   % 0.05 + 0.2*(3 - g) instead of 0.05 + 0.2*(g - 1)
        trains{i} = sort([trains{i}(~late); t*1000]);
    end
end
C_same = C;
C = cNECompare(struct('spikeTimestamp', {trains}, 'tRange', REF, 'NetworkBurstDetected', NB, 'cNE', C_same.reference.cNE), ...
    struct('spikeTimestamp', {trains}, 'tRange', TEST, 'NetworkBurstDetected', NB), FAST);
S = C.summary;
passed = report(sprintf("order reversed: %d cNEs found again, agreement %.2f (p %.4f), template match %.2f", S.nFound, ...
    S.referenceAgreement, S.referenceAgreementP, S.referenceMatch), S.nFound == 3 && S.referenceAgreement < 0.15 && ...
    S.referenceAgreementP > 0.99 && S.referenceMatch == -1) && passed;
passed = report("reuse of a given detection: same reference report", isequaln(C.reference.report.summary, ...
    C_same.reference.report.summary) && isequaln(C.reference.report.cNEs, C_same.reference.report.cNEs)) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Silenced test window, no reference cNEs, reuse of a detection, units that differ

silenced = cellfun(@(t) t(t < TEST(1)*1000), trains, 'UniformOutput', false);
[C2, fig] = cNECompare(struct('spikeTimestamp', {silenced}, 'tRange', REF, 'NetworkBurstDetected', NB, 'cNE', C.reference.cNE), ...
    struct('spikeTimestamp', {silenced}, 'tRange', TEST), FAST);
close(fig);
S = C2.summary;
passed = report(sprintf("silenced test window: no activations (active %.0f%%), no cNEs detected, all %d silent members", S.activeShare, ...
    sum(C2.cNEs.silentMembers)), S.activeShare == 0 && S.nCNETest == 0 && S.nFound == 0 && S.nLost == 3 && ...
    all(C2.cNEs.memberRateTest == 0) && sum(C2.cNEs.silentMembers) == sum(C.cNEs.nMembers)) && passed;

rng(3);
flat = cellfun(@sort, background(N, T_SEC), 'UniformOutput', false);
[C3, fig] = cNECompare(struct('spikeTimestamp', {flat}, 'tRange', REF), struct('spikeTimestamp', {flat}, 'tRange', TEST), FAST);
close(fig);
passed = report(sprintf("no reference cNEs: empty comparison (%d cNEs, found share %g)", C3.summary.nCNE, C3.summary.foundShare), ...
    C3.summary.nCNE == 0 && height(C3.cNEs) == 0 && isnan(C3.summary.foundShare)) && passed;

st = struct('timeStamp', trains, 'cluster', num2cell((1:N)'));
st2 = st;   st2(5).cluster = 99;
try
    cNECompare(struct('spikeTimestamp', st, 'tRange', REF), struct('spikeTimestamp', st2, 'tRange', TEST), FAST);
    refused = false;
catch
    refused = true;
end
passed = report("units that differ between the conditions are refused", refused) && passed;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = report(name, ok)

if ok
    fprintf("  PASS  %s\n", name);
else
    fprintf("  FAIL  %s\n", name);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trains = background(n_units, T_sec)
% independent Poisson units (1-4 Hz), spike times in msec
rates = 1 + 3*rand(n_units,1);
trains = cell(n_units,1);
for i=1:n_units
    trains{i} = rand(poissrnd(rates(i)*T_sec),1)*T_sec*1000;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = events(rate, range)
% Poisson event times (sec) in the range
t = sort(range(1) + 0.01 + rand(poissrnd(rate*diff(range)),1)*(diff(range) - 0.02));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trains = plant(trains, group, t, p)
% at each event time t (sec), each member fires with probability p, within +/- 2 msec
for i=group
    part = rand(numel(t),1) < p;
    trains{i} = [trains{i}; t(part)*1000 + (2*rand(sum(part),1) - 1)*2];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NB = network(s)
% output of "NetworkBurstDetection" (option 5) with large events [s, s + 1] and tails of 1 s
E = struct('start', num2cell(s), 'last', num2cell(s + 1), 'peak', num2cell(s + 0.3), 'level', num2cell(ones(numel(s),1)));
NB = struct('nBurst', numel(E), 'IndivBurst', E, 'closedZones', [s, s + 1, s + 2]);
end
