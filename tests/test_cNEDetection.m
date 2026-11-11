function passed = test_cNEDetection(code_dir)

% Function-level test of "cNEDetection" on synthetic spike trains (independent Poisson units + planted ensembles)
% (October 2026) Unlike the other tests, it does not compare two versions of the code: it checks that
%   1. the planted ensembles are recovered (number of cNEs and members, also a unit shared by two ensembles)
%   2. the planted events are detected
%   3. every output matches a dense, brute-force implementation of the same method (zscore, corr, circshift, prctile; same random numbers)
%   4. same seed = same result; structure array or cell array input = same result
%   5. edge cases: silent unit, spikes outside the window, no cNEs
%   6. chance level of the activations: below 25% of the activations of planted ensembles, no effect on the cNEs and their events
% (November 2026) also the robust definition (default parameters): units below MIN_RATE left out, bound from shifted data and the
% member-count activation rule (both against the dense implementation), each resample of the consensus = a single run on its blocks,
% stability of the planted ensembles, at least MIN_MEMBERS members
% (November 2026) member-count rule with the interval-jitter null (default): against a dense implementation, and on synthetic slow
% co-modulation (units that share seconds-long high-rate states without precise timing: chance-level activations with "jitter" where
% the rate is constant within the jitter window, many with "shift")
%
% Example (from a terminal):
%   matlab -batch "addpath('<repo>/tests'); test_cNEDetection('<repo>')"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dir:             folder with the code to be tested

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% passed:               =true if all the checks passed (each check is printed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(char(code_dir));
T_SEC = 600;                                                                % length of the synthetic recording, in sec
GROUPS = {1:6, 6:11, 25:31};                                                % planted ensembles (unit 6 belongs to two of them)
SEE_ET_AL = struct('MIN_RATE', 0, 'N_CNE_RULE', "MP", 'N_RESAMPLE', 0, 'MIN_MEMBERS', 0, 'MEMBER_COUNT_RULE', false, 'BIN_SIZE', 10);
passed = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1-2. Planted ensembles and events (default parameters: consensus of resamples)

[trains, events] = synthetic_trains(40, T_SEC, GROUPS, 0.3, 0.6, 2, 1);    % 40 units; events at 0.3 Hz, 60% participation, +/-2 ms jitter
trains{35} = [];                                                            % silent unit
trains{36} = [trains{36}; -500; T_SEC*1000 + 500];                          % spikes outside the window
trains{37} = sort(rand(60,1)*T_SEC*1000);                                   % unit firing at 0.1 Hz (below MIN_RATE)
tic; cNE = cNEDetection(trains, [0 T_SEC]); el = toc;
fprintf("  (default parameters: %.0f sec)\n", el);

found = cellfun(@(m) sort(m(:))', cNE.members, 'UniformOutput', false);
recovered = cNE.nCNE == numel(GROUPS) && all(cellfun(@(g) any(cellfun(@(m) isequal(m, g), found)), GROUPS));
passed = report("planted ensembles recovered (members " + strjoin(cellfun(@mat2str, found, 'UniformOutput', false), ", ") + ")", recovered) && passed;
passed = report("planted ensembles found in every resample (stability " + mat2str(cNE.stability) + ")", all(cNE.stability == 1)) && passed;

if recovered
    for g=1:numel(GROUPS)
        k = find(cellfun(@(m) isequal(m, GROUPS{g}), found));
        planted_bin = floor(events{g}/cNE.parameters.BIN_SIZE) + 1;
        hit = mean(arrayfun(@(b) any(abs(cNE.events{k} - b) <= 1), planted_bin));         % planted events detected (+/- 1 bin)
        precision = mean(arrayfun(@(b) any(abs(planted_bin - b) <= 1), cNE.events{k}));   % detected events that were planted
        passed = report(sprintf("ensemble %d: %.0f%% of the planted events detected, %.0f%% of the detected events planted", ...
                        g, 100*hit, 100*precision), hit > 0.6 && precision > 0.5) && passed;
    end
end
passed = report("silent unit left out (weight 0)", ~cNE.active(35) && all(cNE.weights(35,:) == 0)) && passed;
passed = report("unit below MIN_RATE left out (weight 0)", ~cNE.active(37) && all(cNE.weights(37,:) == 0) && cNE.active(36)) && passed;
passed = report("events = bins that pass both activation rules", all(cellfun(@(e, a, c) isequal(e, intersect(a, c)), cNE.events, ...
                cNE.activityEvents, cNE.memberCountEvents))) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Each resample of the consensus = a single run on its blocks (joined)

r = 1;
pick = cNE.resamples(r).blocks;     block = cNE.parameters.RESAMPLE_BLOCK*1000;                % in msec
used = find(cNE.active);
joined = cell(numel(used),1);
for i=1:numel(used)
    t = trains{used(i)};    y = [];
    for j=1:numel(pick)
        y = [y; t(t >= (pick(j)-1)*block & t < pick(j)*block) - (pick(j)-1)*block + (j-1)*block];   %#ok<AGROW>
    end
    joined{i} = y;
end
single = cNEDetection(joined, [0 numel(pick)*block/1000], struct('N_RESAMPLE', 0, 'MIN_RATE', 0, 'MIN_MEMBERS', 0, ...
                      'MEMBER_COUNT_RULE', false, 'SEED', 1000*cNE.parameters.SEED + r));
same_members = isequal(cellfun(@(m) used(m), single.members, 'UniformOutput', false), cNE.resamples(r).members);
passed = report(sprintf("resample %d (blocks %s) = single run on its blocks", r, mat2str(pick)), single.nCNE == cNE.resamples(r).nCNE && ...
                abs(single.lambdaMax - cNE.resamples(r).lambdaMax) < 1e-10 && max(abs(single.weights - cNE.resamples(r).weights(used,:)), [], 'all') < 1e-10 ...
                && same_members) && passed;
passed = report("every cNE has at least MIN_MEMBERS members", all(cellfun(@numel, cNE.members) >= cNE.parameters.MIN_MEMBERS)) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Same results as the dense implementation (single runs)

% the steps of See et al.
opts = SEE_ET_AL;   opts.N_SHUFFLE = 20;
cNE = cNEDetection(trains, [0 T_SEC], opts);
ref = dense_reference(trains, [0 T_SEC], 10, 20, 0, 1.5, 99.9, struct('rule', "MP"));
passed = compare_dense(cNE, ref, "See et al.") && passed;

% robust definition, single run: bound from shifted data, units below MIN_RATE left out, member-count rule (both nulls; the rule is the
% last step, so the rest of the run is the same with both)
opts2 = struct('N_RESAMPLE', 0, 'MIN_MEMBERS', 0, 'N_SHUFFLE', 20, 'N_SHIFT_EIG', 20, 'N_SHIFT_COUNT', 50, 'COUNT_NULL', "shift");
cNE2 = cNEDetection(trains, [0 T_SEC], opts2);
ref2 = dense_reference(trains, [0 T_SEC], 15, 20, 0, 1.5, 99.9, struct('rule', "shift", 'n_shift_eig', 20, 'eig_percentile', 99, ...
                       'min_rate', 0.5, 'n_shift_count', 50, 'count_percentile', 99.9, 'min_firing', 2, 'jitter_bins', 333));
passed = compare_dense(cNE2, ref2, "shifted bound, MIN_RATE") && passed;
passed = report("shifted bound, MIN_RATE: member-count rule, circular shifts (members needed and events)", isequal(cNE2.membersNeeded, ref2.membersNeeded) && ...
                isequal(cNE2.memberCountEvents, ref2.memberCountEvents) && isequal(cNE2.events, ref2.events)) && passed;
cNE3 = cNEDetection(trains, [0 T_SEC], setfield(opts2, 'COUNT_NULL', "jitter"));                                       %#ok<SFLD>
same_run = isequaln(rmfield(cNE3, {'membersNeeded', 'memberCountEvents', 'events', 'eventTimes', 'chanceEvents', 'parameters'}), ...
                    rmfield(cNE2, {'membersNeeded', 'memberCountEvents', 'events', 'eventTimes', 'chanceEvents', 'parameters'}));
passed = report(sprintf("shifted bound, MIN_RATE: member-count rule, 5-s interval jitter (events; %d of %d bins of the circular-shift rule kept)", ...
                numel(vertcat(cNE3.events{:})), numel(vertcat(cNE2.events{:}))), same_run && all(isnan(cNE3.membersNeeded)) && ...
                isequal(cNE3.memberCountEvents, ref2.jitterCountEvents) && isequal(cNE3.events, ref2.jitterEvents)) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Slow co-modulation: units that share seconds-long high-rate states, without precise timing
% NOTE: such a group also raises the shifted-data bound for every component (circular shifts keep each unit's long high-rate periods,
% which then overlap by chance): with a planted precise ensemble of 7 units (eigenvalue 1.23) added to this recording, the bound was
% 1.81-1.88 and only this group was found (bound 1.06 and the ensemble found without the group). The ensemble is left out here

[slow, up] = comodulated_trains(T_SEC, 1:8, 7);
opts4 = struct('N_RESAMPLE', 0, 'N_SHUFFLE', 20, 'N_SHIFT_EIG', 20, 'N_SHIFT_COUNT', 50);
jit = cNEDetection(slow, [0 T_SEC], opts4);
shf = cNEDetection(slow, [0 T_SEC], setfield(opts4, 'COUNT_NULL', "shift"));                                           %#ok<SFLD>
n_state = cellfun(@(m) numel(intersect(m, 1:8)), jit.members);
[~, ks] = max([n_state 0]);
found_state = ~isempty(n_state) && n_state(ks) >= 6;
passed = report(sprintf("slow co-modulation: the units sharing the high-rate states form a cNE (%d of 8 units)", max([n_state 0])), found_state) && passed;
if found_state
    t = jit.binTime(:);     up_bin = false(size(t));
    for j=1:size(up,1)
        up_bin = up_bin | (t >= up(j,1) & t < up(j,2));
    end
    w = ceil((1:numel(t))'/round(jit.parameters.JITTER_WINDOW*1000/jit.parameters.BIN_SIZE));   % jitter window of each bin
    inside = accumarray(w, up_bin, [], @all);   inside = find(up_bin & inside(w));                  % state bins whose window is inside the state
    f_shift = mean(ismember(inside, shf.events{ks}));   f_jit = mean(ismember(inside, jit.events{ks}));
    f_all = mean(ismember(find(up_bin), jit.events{ks}));
    passed = report(sprintf(['slow co-modulation: active in %.1f%% of the high-rate bins with circular shifts and %.2f%% with interval jitter ' ...
                    '(bins whose jitter window is inside a state, chance level 0.1%%; %.2f%% of all the high-rate bins)'], ...
                    100*f_shift, 100*f_jit, 100*f_all), f_shift > 0.01 && f_jit < 0.005 && f_jit < f_shift/10) && passed;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4-5. Reproducibility, input formats, no cNEs

spikeTimestamp = struct('timeStamp', trains);
passed = report("same seed = same result; structure array = cell array", isequaln(cNE, cNEDetection(spikeTimestamp, [0 T_SEC], opts))) && passed;
other = cNEDetection(trains, [0 T_SEC], setfield(opts, 'SEED', 1));                                                   %#ok<SFLD>
member_sets = @(c) sort(cellfun(@(m) mat2str(sort(m(:))'), c.members, 'UniformOutput', false));   % (the order of the ICs is arbitrary)
passed = report("another seed: same members", isequal(member_sets(cNE), member_sets(other))) && passed;
fast = struct('N_SHUFFLE', 20);
consensus = cNEDetection(trains, [0 T_SEC], fast);
passed = report("consensus: same seed = same result; structure array = cell array", isequaln(consensus, cNEDetection(spikeTimestamp, [0 T_SEC], fast))) && passed;
no_chance = cNEDetection(trains, [0 T_SEC], setfield(fast, 'N_CHANCE', 0));                                           %#ok<SFLD>
share = consensus.chanceEvents./cellfun(@numel, consensus.events);
passed = report(sprintf("chance level: %s%% of the activations of the planted ensembles; the cNEs and events do not change; N_CHANCE = 0: NaN", ...
    mat2str(round(100*share))), all(share < 0.25) && isequal(no_chance.events, consensus.events) && isequal(no_chance.weights, consensus.weights) ...
    && all(isnan(no_chance.chanceEvents))) && passed;

% independent units (no ensembles). NOTE: with the steps of See et al., the Marcenko-Pastur bound is exceeded by chance in some independent
% recordings (a marginal eigenvalue, ~1 in 5 synthetic recordings like this one); the shifted-data bound and the consensus remove them
independent = synthetic_trains(60, T_SEC, {}, 0, 0, 0, 101);
none = cNEDetection(independent, [0 T_SEC], fast);
passed = report("independent units: no cNEs", none.nCNE == 0 && isempty(none.members) && size(none.activity,1) == 0) && passed;
none = cNEDetection(independent, [0 T_SEC], setfield(opts, 'N_SHUFFLE', 20));                                          %#ok<SFLD>
passed = report("independent units, steps of See et al.: no cNEs", none.nCNE == 0 && isempty(none.members) && size(none.activity,1) == 0) && passed;

if passed
    fprintf("\ncNEDetection: all checks passed\n");
else
    fprintf("\ncNEDetection: SOME CHECKS FAILED\n");
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function passed = compare_dense(cNE, ref, name)

tol = 1e-8;
passed = report(name + ": eigenvalues, bound and number of cNEs", cNE.nCNE == ref.nCNE && max(abs(cNE.eigenvalues - ref.eigenvalues)) < tol ...
                && abs(cNE.lambdaMax - ref.lambdaMax) < tol);
passed = report(name + ": units analyzed", isequal(cNE.active, ref.active)) && passed;
passed = report(name + ": IC weights", max(abs(cNE.weights(cNE.active,:) - ref.weights), [], 'all') < tol) && passed;
passed = report(name + ": members", isequal(cNE.members, ref.members)) && passed;
passed = report(name + ": cNE activity", max(abs(cNE.activity' - ref.activity), [], 'all') < tol*max(abs(ref.activity(:)))) && passed;
passed = report(name + ": membership thresholds", max(abs(cNE.memberThreshold - ref.memberThreshold), [], 'all') < tol) && passed;
passed = report(name + ": activity thresholds", max(abs(cNE.activityThreshold - ref.activityThreshold)) < tol*max(ref.activityThreshold)) && passed;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trains, events] = synthetic_trains(n_units, T_sec, groups, event_rate, p_participate, jitter, seed)

% independent Poisson units (1-4 Hz) plus planted ensembles: at each event of an ensemble (Poisson, event_rate Hz) each member fires
% with probability p_participate, within +/- jitter msec of the event. Spike times in msec

rng(seed);
rates = 1 + 3*rand(n_units,1);
trains = cell(n_units,1);
for i=1:n_units
    trains{i} = rand(poissrnd(rates(i)*T_sec),1)*T_sec*1000;
end
events = cell(numel(groups),1);
for g=1:numel(groups)
    events{g} = sort(rand(poissrnd(event_rate*T_sec),1)*(T_sec*1000 - 20) + 10);
    for i=groups{g}
        part = rand(numel(events{g}),1) < p_participate;
        trains{i} = [trains{i}; events{g}(part) + (2*rand(sum(part),1) - 1)*jitter];
    end
end
trains = cellfun(@sort, trains, 'UniformOutput', false);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trains, up] = comodulated_trains(T_sec, state_units, seed)

% 30 independent Poisson units (1-4 Hz). The state units fire at 1 Hz, and at 30 Hz (Poisson, independently) during shared 15-s
% high-rate states that start every 40 sec, longer than the jitter windows (up: [start end] of each state, in msec). Spike times in msec

rng(seed);
trains = cell(30,1);
for i=1:30
    trains{i} = rand(poissrnd((1 + 3*rand)*T_sec),1)*T_sec*1000;
end
up = (5:40:T_sec-16)'*1000 + [0 15000];
for i=state_units
    t = rand(poissrnd(T_sec),1)*T_sec*1000;
    for j=1:size(up,1)
        t = [t; up(j,1) + rand(poissrnd(30*15),1)*15000];                   %#ok<AGROW>
    end
    trains{i} = t;
end
trains = cellfun(@sort, trains, 'UniformOutput', false);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ref = dense_reference(trains, tRange, bin_size, n_shuffle, seed, member_sd, percentile, opt)

% the same method (single run) with dense matrices and the standard functions, and the same random numbers as "cNEDetection" (default
% fastICA settings). opt.rule: "MP" or "shift" (opt.n_shift_eig, opt.eig_percentile); opt.min_rate (Hz, default 0); member-count rule if
% opt.n_shift_count is given (opt.count_percentile, opt.min_firing), and also with the interval-jitter null if opt.jitter_bins is given
% (exact distribution of the count by convolution, bin by bin)

edges = tRange(1)*1000 : bin_size : tRange(2)*1000;
n_bin = numel(edges) - 1;
counts = zeros(n_bin, numel(trains));
for i=1:numel(trains)
    counts(:,i) = histcounts(trains{i}, edges);
end
min_rate = 0;
if isfield(opt, 'min_rate'), min_rate = opt.min_rate; end
ref.active = (std(counts) > 0 & sum(counts,1)/(n_bin*bin_size/1000) >= min_rate)';
A = counts(:, ref.active);
Z = zscore(A);
n_active = size(Z,2);
[E, D] = eig(corr(Z));
[lambda, order] = sort(diag(D), 'descend');     E = E(:,order);
ref.eigenvalues = lambda;

stream = RandStream("twister", "Seed", seed);
if opt.rule == "shift"
    top = zeros(opt.n_shift_eig, 1);
    for it=1:opt.n_shift_eig
        shift = randi(stream, n_bin, n_active, 1) - 1;
        Z_shift = Z;
        for i=1:n_active
            Z_shift(:,i) = circshift(Z(:,i), shift(i));
        end
        top(it) = max(eig(corr(Z_shift)));
    end
    ref.lambdaMax = prctile(top, opt.eig_percentile);
else
    ref.lambdaMax = (1 + sqrt(n_active/n_bin))^2;
end
ref.nCNE = sum(lambda > ref.lambdaMax);
K = ref.nCNE;

ref.weights = dense_ica(Z, E(:,1:K), lambda(1:K), stream);
ref.activity = dense_activity(Z, ref.weights);
null_weights = [];      null_activity = zeros(n_bin*n_shuffle, K);
for it=1:n_shuffle
    shift = randi(stream, n_bin, n_active, 1) - 1;
    Z_shuffle = Z;
    for i=1:n_active
        Z_shuffle(:,i) = circshift(Z(:,i), shift(i));
    end
    [E_shuffle, D_shuffle] = eig(corr(Z_shuffle));
    [lambda_shuffle, order] = sort(diag(D_shuffle), 'descend');
    V = dense_ica(Z_shuffle, E_shuffle(:,order(1:K)), lambda_shuffle(1:K), stream);
    null_weights = [null_weights; V(:)];                                    %#ok<AGROW>
    null_activity((it-1)*n_bin + (1:n_bin), :) = dense_activity(Z_shuffle, ref.weights);
end
ref.memberThreshold = mean(null_weights) + member_sd*std(null_weights)*[-1 1];
ref.activityThreshold = prctile(null_activity, percentile, 1)';

idx = find(ref.active);
ref.members = cell(1,K);    member_col = cell(1,K);
for k=1:K
    [w, order] = sort(ref.weights(:,k), 'descend');
    member_col{k} = order(w > ref.memberThreshold(2));
    ref.members{k} = idx(member_col{k});
end

if isfield(opt, 'n_shift_count')                                            % member-count rule (after the null distributions)
    ref.membersNeeded = nan(1,K);   ref.memberCountEvents = cell(1,K);  ref.events = cell(1,K);
    for k=1:K
        F = A(:, member_col{k}) > 0;                                        % members with at least one spike in each bin
        n_mem = size(F,2);
        h = zeros(1, n_mem + 1);
        for s=1:opt.n_shift_count
            shift = randi(stream, n_bin, n_mem, 1);
            F_shift = F;
            for j=1:n_mem
                F_shift(:,j) = circshift(F(:,j), shift(j));
            end
            h = h + histcounts(sum(F_shift, 2), -0.5:1:n_mem+0.5);
        end
        chance = find(cumsum(h)/sum(h) >= opt.count_percentile/100, 1) - 1;
        ref.membersNeeded(k) = max(chance + 1, opt.min_firing);
        ref.memberCountEvents{k} = find(sum(F, 2) >= ref.membersNeeded(k));
        ref.events{k} = intersect(find(ref.activity(:,k) > ref.activityThreshold(k)), ref.memberCountEvents{k});
    end
end
if isfield(opt, 'jitter_bins')                                              % member-count rule, interval-jitter null
    ref.jitterCountEvents = cell(1,K);  ref.jitterEvents = cell(1,K);
    m = opt.jitter_bins;
    for k=1:K
        C = A(:, member_col{k});    n_firing = sum(C > 0, 2);
        pass = false(n_bin, 1);
        for b=find(n_firing >= opt.min_firing)'
            first = floor((b-1)/m)*m + 1;   last = min(first + m - 1, n_bin);              % jitter window of the bin b
            p = 1 - (1 - 1/(last - first + 1)).^sum(C(first:last,:), 1);   % P(each member fires in the bin)
            d = 1;
            for j=1:numel(p)
                d = conv(d, [1 - p(j), p(j)]);                              % distribution of the number of members firing
            end
            pass(b) = sum(d(n_firing(b)+1:end)) < 1 - opt.count_percentile/100;
        end
        ref.jitterCountEvents{k} = find(pass);
        ref.jitterEvents{k} = intersect(find(ref.activity(:,k) > ref.activityThreshold(k)), ref.jitterCountEvents{k});
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = dense_activity(Z, V)

R = zeros(size(Z,1), size(V,2));
for k=1:size(V,2)
    P = V(:,k)*V(:,k)';
    P(1:size(P,1)+1:end) = 0;                                               % zero diagonal
    R(:,k) = sum((Z*P).*Z, 2);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = dense_ica(Z, P, lambda, stream)

% symmetric fastICA (logcosh) on the whitened projection onto the principal components P; unit-length weights, largest one positive

X = (Z*P)/diag(sqrt(lambda));
[n_sample, K] = size(X);
U = decorrelate(randn(stream, K, K));
for it=1:1000
    G = tanh(X*U');
    U_new = decorrelate((G'*X)/n_sample - diag(mean(1 - G.^2, 1))*U);
    change = max(abs(abs(diag(U_new*U')) - 1));
    U = U_new;
    if change < 1e-6
        break
    end
end
V = P*diag(1./sqrt(lambda))*U';
for k=1:K
    V(:,k) = V(:,k)/norm(V(:,k));
    [~, imax] = max(abs(V(:,k)));
    V(:,k) = V(:,k)*sign(V(imax,k));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = decorrelate(U)

[E, D] = eig((U*U' + (U*U')')/2);
U = E*diag(1./sqrt(diag(D)))*E'*U;

end
