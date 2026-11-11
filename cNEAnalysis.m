function [report, fig] = cNEAnalysis(cNE, NetworkBurstDetected, spikeTimestamp, inputs)

% function to describe the network activity of a recording-condition with its coordinated neuronal ensembles (cNEs, "cNEDetection"):
% when each cNE is active relative to the network events, in which order the cNEs are active within the large events, and which cNEs
% the small and the secondary events recruit
% (November 2026) The cNEs do not detect network events; they extend their analysis. Network events of "NetworkBurstDetection" (option 5):
% large (level-1) events, small (level-0) events, secondary events (reverberations, inside the large events, ".Secondary") and the tail of
% each large event (from its end to the end of its closed zone, ".closedZones")
%   1. periods: each bin of the cNE detection belongs to one period, by precedence: secondary (start to end of a secondary event) > core
%      (start to end of a large event) > small (start to end of a small event) > tail > onset (the ONSET_WINDOW sec before the start of a
%      large event) > quiet baseline (the rest). Per cNE: share of its activations and activations per minute in each period
%   2. time of a cNE in an event: TIME_PERCENTILE percentile of the times of its active bins, from the start of the event (with fewer than
%      MIN_BINS_TIME active bins: the midpoint of the first one and the median). Window of a large event: from ONSET_WINDOW sec before its
%      start (only after the end of the previous tail) to the end of its tail (only before the next onset, and never shorter than the
%      event); of a small event: from ONSET_WINDOW sec before its start (outside the large events, their onsets and their tails) to its
%      end (outside the large events); of a secondary event: from its start to its end. A cNE is "active in an event" when it has at least
%      one active bin in that window
%   3. template: median time of each cNE over the large events; rank and position of each cNE in that order (0 = first, 1 = last)
%   4. order consistency of the large events, S = sum_ij n_ij |2 p_ij - 1| / sum_ij n_ij over the pairs of cNEs (n_ij: events where both
%      are active; p_ij: share of them where i comes first, ties count 1/2). Null: the times are permuted among the active cNEs of each
%      event (N_PERM permutations). Reported scaled to chance, (S - S0)/(1 - S0), S0 = mean of the null (0 = chance, 1 = the same order
%      in every event), with p = share of the null at least as large as S. (November 2026) Timing of the order: for each pair of cNEs
%      active together in at least 3 large events, the lag between them in each of these events (later minus earlier cNE of the
%      template); its median (lag) and its median absolute deviation (jitter: how precise the timing is); medians over the pairs
%   5. per-event match: Kendall tau (tau-b) between the order of the active cNEs of a large event and the template of the other large
%      events (events with at least 3 active cNEs that have a time there); median over the events
%   6. small and secondary events vs the template: mean position of the cNEs they recruit (over the events and their active cNEs; 0 = the
%      leading cNEs, 1 = the late ones), and agreement A = share of the pairs of active cNEs in the template order (ties count 1/2; 0.5 =
%      chance, 1 = the same order), with the same null as S; also the consistency of their own order (as in 4; small events: November
%      2026, with its own generator of permutations, so that the other p values do not change). (November 2026) Share of the cNEs that an
%      event recruits (active cNEs / all the cNEs), mean over the events of the type that recruit at least one cNE (NaN: none does)
%   7. member participation: share of the members of a cNE that fire (at least one spike) in each of its active bins (median); it falls
%      with the number of members, so it is meant to compare a cNE with itself (e.g., across conditions)
%   8. roles (a cNE can have several): core = largest share of its activations in the core; leader = active in the onset of at least
%      LEADER_ONSETS % of the large events; tail = at least TAIL_SHARE % of its activations in the tails and active in at least TAIL_COVER %
%      of the tails; reverberation = active in at least REVERBERATION_COVER % of the secondary events; small-event = active in at least
%      SMALL_COVER % of the small events; baseline = largest share of its activations in the quiet baseline. The reverberation and
%      small-event roles need at least MIN_EVENTS_ROLE events of that type
%   9. (November 2026) chance level: activations expected by chance (cNE.chanceEvents of "cNEDetection", interval-jitter surrogates), as a
%      rate and as a share of the activations; excess rate = activation rate - chance rate. A cNE is flagged "at chance" when the chance
%      level is at least MAX_CHANCE_SHARE % of its activations (nothing is left out: the flag tells which activation counts are mostly
%      what the rules let through by chance)
% Minimums: the order measures (4, 5, ranks and positions, and the comparison of 6) need at least MIN_CNE_ORDER cNEs and MIN_LARGE_ORDER
% large events; the comparisons of 6 also need at least MIN_EVENTS_TYPE events with an active cNE (position) or with 2 active cNEs
% (agreement, own order). (November 2026) With 2 cNEs (the minimum), the order is which cNE leads: the consistency (4) tells whether the
% same one leads in every event; the per-event match (5) stays blank (it needs 3 active cNEs in an event)
% Events the cNEs do not describe: when fewer than MIN_COVERAGE % of the large events have an active cNE, the roles, the times of the
% cNEs and the measures of 3-6 are left blank (NaN); summary.notDescribed = true. The periods and the share of the events of each type
% where each cNE is active are still reported
% No large event in the window (summary.noBursting = true): only the cNEs, their activation rates and member participation
% Only the events entirely inside the window of the cNE detection are counted (the others still set the period of the bins they cover)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   cNE:                    output of "cNEDetection" for the condition
%   NetworkBurstDetected:   output of "NetworkBurstDetection" (option 5) for the same condition; [] or nBurst = 0: no network events
%   spikeTimestamp:         (optional) spike trains given to "cNEDetection" (structure array with .timeStamp, in msec, or cell array),
%                           for the member participation; [] = not computed (NaN)
%   inputs:                 (optional) structure that overwrites any of the PARAMETERS below, e.g., struct('N_PERM', 500)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report.summary            one row per condition (structure of scalars; shares and coverages in %, rates per minute):
%   .nUnits, .nCNE              units analyzed by "cNEDetection" (cNE.active) and number of cNEs
%   .medianMembers, .unitsInCNE median number of members, and % of the analyzed units that are members of at least one cNE
%   .activationRate             median over the cNEs of their activations per minute (whole window)
%   .excessRate, .nAtChance     (9) median over the cNEs of their excess rate (per minute), and number of cNEs at chance
%   .memberParticipation        median over the cNEs of their member participation (7)
%   .nLarge, .nSmall, .nSecondary   number of events of each type (small events that overlap a large event are left out: .nSmallLeftOut)
%   .largeWithCNE               % of the large events with at least one active cNE
%   .notDescribed, .noBursting  flags (see above)
%   .timeOnset ... .timeQuiet   % of the window in each period
%   .shareOnset ... .shareQuiet median over the cNEs of the % of their activations in each period (NaN: period absent)
%   .nCore, .nLeader, .nTail, .nReverberation, .nSmallEvent, .nBaseline   number of cNEs with each role
%   .leaderLead                 median over the leaders of their median lead (start of the event minus their first onset activation), in sec
%   .orderConsistency, .orderP  order consistency of the large events scaled to chance (4), and its p
%   .orderMatch                 median per-event match (5)
%   .orderSpan, .orderSpanFraction  time from the first to the last cNE of the template, in sec and as a fraction of the median duration
%                               of the large events
%   .orderLag, .orderJitter     (November 2026) timing of the order (4): median lag between the cNEs of a pair and its median absolute
%                               deviation over the large events (jitter), in sec; medians over the pairs (NaN: no pair active together
%                               in at least 3 large events)
%   .smallWithCNE, .smallPosition, .smallAgreement, .smallAgreementP           small events (6): % with an active cNE, position, A, p
%   .smallRecruited, .secondaryRecruited    (November 2026) small and secondary events (6): mean % of the cNEs that an event recruits,
%                               over the events with an active cNE
%   .smallConsistency, .smallConsistencyP   (November 2026) small events (6): consistency of their own order, and its p
%   .secondaryWithCNE, .secondaryPosition, .secondaryAgreement, .secondaryAgreementP, .secondaryConsistency, .secondaryConsistencyP
%   .referenceAgreement, .referenceAgreementP   (with REFERENCE_TIME) agreement A of the large events with the reference order (as in 6,
%                               with the reference template; at least MIN_CNE_ORDER cNEs and MIN_EVENTS_TYPE events with 2 active cNEs that
%                               have a reference time), and its p; NaN without REFERENCE_TIME
%   .referenceMatch             (with REFERENCE_TIME) Kendall tau-b between the template of this condition and the reference template
%                               (with the order measures of this condition, and at least 3 cNEs with both times)
% report.cNEs               table, one row per cNE:
%   cNE, nMembers, members (index into spikeTimestamp), memberCluster (.cluster of spikeTimestamp, if any; to track the cNE), stability,
%   rate (activations per minute), chanceRate, excessRate (per minute), chanceShare (% of its activations), atChance (9),
%   memberParticipation (median; also participationInLarge and participationOutside the large events),
%   shareOnset ... shareQuiet (% of its activations), rateOnset ... rateQuiet (per minute of each period), largeCover, onsetCover,
%   tailCover, smallCover, secondaryCover (% of the events of each type where it is active; onsets and tails: active in that period),
%   onsetLead (median lead in the onsets, in sec), time (template, in sec from the start of the large events), rank, position, roles
% report.events              per event (kept for the figures and later analyses):
%   .large: start, last, peak, tailEnd, window ([start end] of the window of 2), time(j,k) (time of the cNE k in the event j, NaN = not
%           active), nBins(j,k) (active bins), match (5); .peth(l,k), .lags(l): share of the large events where the cNE k is active at the
%           lag l from their start (bins outside the windows left out); .peakTime, .coreEnd, .windowEnd (medians, from the start)
%   .small: start, last, inTail (inside the tail of a large event), time, nBins
%   .secondary: start, last, time, nBins
% report.binPeriod(b)       period of the bin b of the cNE detection (1-6, in the order of report.periods)
% report.periods            names of the periods: "onset" "core" "secondary" "tail" "small" "quiet"
% report.parameters         parameters used
% fig                       (only if requested) figure with the shares of the periods, the activity around the large events and the
%                           participation in each type of event

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ONSET_WINDOW = 0.5;              % onset of the large events: this time before their start, in sec (also the lead-in of the small events)
TIME_PERCENTILE = 25;            % time of a cNE in an event: this percentile of the times of its active bins
MIN_BINS_TIME = 4;               % with fewer active bins: the midpoint of the first one and the median
N_PERM = 2000;                   % permutations of the nulls of the order consistency and of the agreement
SEED = 1;                        % seed of the permutations (same seed + same inputs = same result)
MIN_COVERAGE = 50;               % the cNEs describe the network events if at least this % of the large events have an active cNE
MIN_CNE_ORDER = 2;               % order measures: minimum number of cNEs (November 2026: 2, was 3)
MIN_LARGE_ORDER = 5;             % order measures: minimum number of large events
MIN_EVENTS_TYPE = 3;             % small / secondary events vs the template: minimum number of events with an active cNE (position) or with
                                 % 2 active cNEs (agreement, own order)
MIN_EVENTS_ROLE = 3;             % reverberation and small-event roles: minimum number of events of that type
LEADER_ONSETS = 25;              % leader: active in the onset of at least this % of the large events
TAIL_SHARE = 30;                 % tail: at least this % of its activations in the tails ...
TAIL_COVER = 50;                 %       ... and active in at least this % of the tails
REVERBERATION_COVER = 50;        % reverberation: active in at least this % of the secondary events
SMALL_COVER = 50;                % small-event: active in at least this % of the small events
REFERENCE_TIME = [];             % (November 2026, "cNECompare") template of the same cNEs in a reference condition (report.cNEs.time of
                                 % that condition, in sec; cNEs tracked with "cNEDetection"): the large events of this condition are
                                 % compared with that order ([] = not compared)
MAX_CHANCE_SHARE = 50;           % (November 2026) a cNE is "at chance" when the activations expected by chance are at least this % of its
                                 % activations
apply_parameters('cNEAnalysis');                                           % values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")
if nargin > 3 && ~isempty(inputs)
    override_inputs(inputs);
end
if nargin < 3
    spikeTimestamp = [];
end
parameters = struct('ONSET_WINDOW', ONSET_WINDOW, 'TIME_PERCENTILE', TIME_PERCENTILE, 'MIN_BINS_TIME', MIN_BINS_TIME, 'N_PERM', N_PERM, ...
    'SEED', SEED, 'MIN_COVERAGE', MIN_COVERAGE, 'MIN_CNE_ORDER', MIN_CNE_ORDER, 'MIN_LARGE_ORDER', MIN_LARGE_ORDER, ...
    'MIN_EVENTS_TYPE', MIN_EVENTS_TYPE, 'MIN_EVENTS_ROLE', MIN_EVENTS_ROLE, 'LEADER_ONSETS', LEADER_ONSETS, 'TAIL_SHARE', TAIL_SHARE, ...
    'TAIL_COVER', TAIL_COVER, 'REVERBERATION_COVER', REVERBERATION_COVER, 'SMALL_COVER', SMALL_COVER, 'REFERENCE_TIME', REFERENCE_TIME, ...
    'MAX_CHANCE_SHARE', MAX_CHANCE_SHARE);

K_MSEC_SEC = 1/1000;                                                        % factor to rescale from msec to sec
PERIODS = ["onset" "core" "secondary" "tail" "small" "quiet"];
W = ONSET_WINDOW;
stream = RandStream('mt19937ar', 'Seed', SEED);                             % permutations (the global generator is not touched)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bins and activations of the cNEs

tb = double(cNE.binTime(:))/1000;                                           % midpoint of each bin, in sec
nb = numel(tb);
bin_size = cNE.parameters.BIN_SIZE*K_MSEC_SEC;                              % in sec
bin_min = bin_size/60;                                                      % in min
if isfield(cNE.parameters, 'tRange')
    win = cNE.parameters.tRange;                                            % window of the cNE detection, in sec
else
    win = [tb(1) tb(end)] + [-1 1]*bin_size/2;
end
nk = cNE.nCNE;
act = false(nb, nk);                                                        % active bins of each cNE
for k=1:nk
    act(cNE.events{k}, k) = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Network events and periods

[Lall, pk_all, S0all, S2all, Z] = network_events(NetworkBurstDetected);     % [start last] of the large, small and secondary events, closed zones
[~, o] = sort(Lall(:,1));   Lall = Lall(o,:);   pk_all = pk_all(o);
nA = size(Lall,1);
ov = Lall(:,2) >= Z(:,1)' & Lall(:,1) <= Z(:,2)';                           % closed zones that overlap each large event
tail_end = Lall(:,2);
for q=1:nA
    tail_end(q) = max([Z(ov(q,:),3); Lall(q,2)]);                           % end of its tail (= its end without a closed zone)
end
inside = @(S) S(:,1) >= win(1) & S(:,2) <= win(2);                         % events entirely inside the window
iL = find(inside(Lall));    L = Lall(iL,:);     nL = numel(iL);

in_core = in_intervals(tb, Lall);
label = 6*ones(nb,1);                                                       % quiet baseline
label(in_intervals(tb, [Lall(:,1) - W, Lall(:,1)]) & ~in_core) = 1;        % onset
label(in_intervals(tb, [Lall(:,2), tail_end]) & ~in_core) = 4;             % tail
label(in_intervals(tb, S0all)) = 5;                                         % small events
label(in_core) = 2;                                                         % core
label(in_intervals(tb, S2all)) = 3;                                         % secondary events
period_min = accumarray(label, 1, [6 1])'*bin_min;                          % time in each period, in min

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time of each cNE in each event

w0 = Lall(:,1) - W;     w1 = tail_end;                                      % windows of the large events
for q=2:nA
    w0(q) = max(w0(q), tail_end(q-1));                                      % the onset only after the previous tail
end
for q=1:nA-1
    w1(q) = min(w1(q), Lall(q+1,1) - W);                                    % up to the next onset
end
w1 = max(w1, Lall(:,2));                                                    % (never shorter than the event)
XL = nan(nL, nk);   nbinL = zeros(nL, nk);
for j=1:nL
    q = iL(j);
    [XL(j,:), nbinL(j,:)] = event_times(act, tb, find(tb >= w0(q) & tb <= w1(q)), Lall(q,1), TIME_PERCENTILE, MIN_BINS_TIME);
end

lags = zeros(1,0);  peth = zeros(0, nk);                                    % activity around the large events
if nL > 0
    lags = -W : bin_size : prctile(w1(iL) - L(:,1), 75);
    nl = numel(lags);   PE = nan(nL, nl, nk);
    for j=1:nL
        q = iL(j);
        sel = tb >= w0(q) & tb <= w1(q);
        idx = round((lags + Lall(q,1) - tb(1))/bin_size) + 1;     ok = idx >= 1 & idx <= nb;
        in_w = false(1, nl);    in_w(ok) = sel(idx(ok));
        for k=1:nk
            v = nan(1, nl);     v(ok) = act(idx(ok), k);    v(~in_w) = NaN;
            PE(j,:,k) = v;
        end
    end
    peth = reshape(mean(PE, 1, 'omitnan'), nl, nk);
end

keep0 = true(size(S0all,1), 1);                                             % small events that overlap a large event are left out
for q=1:size(S0all,1)
    keep0(q) = ~any(S0all(q,2) >= Lall(:,1) & S0all(q,1) <= Lall(:,2));
end
i0 = find(keep0 & inside(S0all));   n0 = numel(i0);
in_big = in_intervals(tb, [Lall(:,1) - W, tail_end]);
X0 = nan(n0, nk);   nbin0 = zeros(n0, nk);  tail0 = false(n0,1);
for j=1:n0
    s = S0all(i0(j),1);     e = S0all(i0(j),2);
    sel = find((tb >= s - W & tb < s & ~in_big) | (tb >= s & tb <= e & ~in_core));
    [X0(j,:), nbin0(j,:)] = event_times(act, tb, sel, s, TIME_PERCENTILE, MIN_BINS_TIME);
    tail0(j) = any(s <= tail_end & e >= Lall(:,2));
end
i2 = find(inside(S2all));   n2 = numel(i2);
X2 = nan(n2, nk);   nbin2 = zeros(n2, nk);
for j=1:n2
    s = S2all(i2(j),1);     e = S2all(i2(j),2);
    [X2(j,:), nbin2(j,:)] = event_times(act, tb, find(tb >= s & tb <= e), s, TIME_PERCENTILE, MIN_BINS_TIME);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Validity of the measures

large_with = 100*mean(any(~isnan(XL), 2));                                  % NaN without large events
no_bursting = nL == 0;
not_described = ~no_bursting && large_with < MIN_COVERAGE;
events_ok = ~no_bursting && ~not_described;                                 % roles and measures of the network events
order_ok = events_ok && nk >= MIN_CNE_ORDER && nL >= MIN_LARGE_ORDER;       % order measures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Member participation

[part, part_in, part_out] = member_participation(cNE, spikeTimestamp, in_core);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Each cNE: periods, events, roles

share = nan(nk, 6);     rate_p = nan(nk, 6);    rate = nan(nk,1);
onset_cover = nan(nk,1);    onset_lead = nan(nk,1);     tail_cover = nan(nk,1);
chance = nan(nk,1);     chance_share = nan(nk,1);                           % (9) chance level (NaN: not computed by "cNEDetection")
if isfield(cNE, 'chanceEvents') && numel(cNE.chanceEvents) == nk
    chance = cNE.chanceEvents(:);
end
roles = strings(nk,1);
for k=1:nk
    e = cNE.events{k}(:);   lab = label(e);     t = tb(e);  n = numel(e);
    cnt = accumarray(lab, 1, [6 1])';
    share(k,:) = 100*cnt/max(n,1);      rate_p(k,:) = cnt./period_min;
    share(k, period_min == 0) = NaN;    rate_p(k, period_min == 0) = NaN;   % period absent
    rate(k) = n/(nb*bin_min);
    if n > 0
        chance_share(k) = 100*chance(k)/n;
    end
    first = nan(nL,1);  hit_tail = false(nL,1);
    for j=1:nL
        q = iL(j);
        ta = t(lab == 1 & t >= Lall(q,1) - W & t < Lall(q,1));
        if ~isempty(ta)
            first(j) = Lall(q,1) - min(ta);                                 % lead of its first onset activation
        end
        hit_tail(j) = any(lab == 4 & t > Lall(q,2) & t <= tail_end(q));
    end
    onset_cover(k) = 100*mean(~isnan(first));   onset_lead(k) = median(first, 'omitnan');   tail_cover(k) = 100*mean(hit_tail);
end
chance_rate = chance/(nb*bin_min);  excess_rate = rate - chance_rate;     % (9) chance level, per minute
at_chance = chance_share >= MAX_CHANCE_SHARE;                               % (false when not computed, or without activations)
large_cover = 100*mean(~isnan(XL), 1)';     small_cover = 100*mean(~isnan(X0), 1)';    secondary_cover = 100*mean(~isnan(X2), 1)';
if nk == 0
    [large_cover, small_cover, secondary_cover] = deal(zeros(0,1));
end
if no_bursting
    [large_cover, small_cover, secondary_cover] = deal(nan(nk,1));
end

is_role = false(nk, 6);                                                     % core, leader, tail, reverberation, small-event, baseline
ROLE_NAMES = ["core" "leader" "tail" "reverberation" "small-event" "baseline"];
if events_ok
    for k=1:nk
        if isempty(cNE.events{k})
            continue
        end
        [~, mx] = max(share(k,:));                                          % (absent periods are NaN and ignored)
        is_role(k,:) = [mx == 2, onset_cover(k) >= LEADER_ONSETS, share(k,4) >= TAIL_SHARE && tail_cover(k) >= TAIL_COVER, ...
            n2 >= MIN_EVENTS_ROLE && secondary_cover(k) >= REVERBERATION_COVER, n0 >= MIN_EVENTS_ROLE && small_cover(k) >= SMALL_COVER, mx == 6];
        roles(k) = strjoin(ROLE_NAMES(is_role(k,:)), ", ");
    end
end

tmpl = median(XL, 1, 'omitnan');                                            % template (NaN: never active in the large events)
if nL == 0
    tmpl = nan(1, nk);
end
[rk, pos] = template_rank(tmpl);
if ~events_ok
    tmpl(:) = NaN;
end
if ~order_ok
    rk(:) = NaN;    pos(:) = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Order of the large events, and the small and secondary events vs that order

[order_S, order_p, order_match, span, span_fraction, order_lag, order_jitter] = deal(NaN);
tau = nan(nL,1);
if order_ok
    [S, order_p, S0] = consistency(XL, N_PERM, stream);
    order_S = scale_to_chance(S, S0);
    tau = loo_tau(XL);
    order_match = median(tau, 'omitnan');
    span = max(tmpl) - min(tmpl);
    span_fraction = span/median(L(:,2) - L(:,1));
    [order_lag, order_jitter] = pair_lags(XL, tmpl);
end
sm = event_type(X0, tmpl, pos, events_ok, order_ok, false, nk >= MIN_CNE_ORDER, MIN_EVENTS_TYPE, N_PERM, stream);
se = event_type(X2, tmpl, pos, events_ok, order_ok, true, nk >= MIN_CNE_ORDER, MIN_EVENTS_TYPE, N_PERM, stream);
sm_own = event_type(X0, tmpl, pos, events_ok, false, true, nk >= MIN_CNE_ORDER, MIN_EVENTS_TYPE, N_PERM, ...
    RandStream('mt19937ar', 'Seed', SEED + 1));     % (November 2026) own order of the small events only; own generator: other p unchanged
[ref_A, ref_pA, ref_tau] = deal(NaN);                                       % large events vs the order of a reference condition
if ~isempty(REFERENCE_TIME)
    ref_time = reshape(double(REFERENCE_TIME), 1, []);
    if numel(ref_time) ~= nk
        error("cNEAnalysis: REFERENCE_TIME has %d values, there are %d cNEs", numel(ref_time), nk);
    end
    if events_ok && nk >= MIN_CNE_ORDER && sum(sum(~isnan(XL) & ~isnan(ref_time), 2) >= 2) >= MIN_EVENTS_TYPE
        [ref_A, ref_pA] = agreement(XL, ref_time, N_PERM, stream);
    end
    both = ~isnan(tmpl) & ~isnan(ref_time);
    if order_ok && sum(both) >= 3
        ref_tau = kendall_b(tmpl(both), ref_time(both));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output

n_members = cellfun(@numel, cNE.members(:));
all_members = cellfun(@(m) m(:), cNE.members(:), 'UniformOutput', false);
summary = struct();
summary.nUnits = sum(cNE.active);
summary.nCNE = nk;
summary.medianMembers = median(n_members);
summary.unitsInCNE = 100*numel(unique(vertcat(all_members{:})))/max(summary.nUnits, 1);
summary.activationRate = median(rate);
summary.excessRate = median(excess_rate);
summary.nAtChance = sum(at_chance);
summary.memberParticipation = median(part);
summary.nLarge = nL;
summary.nSmall = n0;
summary.nSmallLeftOut = sum(~keep0 & inside(S0all));
summary.nSecondary = n2;
summary.largeWithCNE = large_with;
summary.notDescribed = not_described;
summary.noBursting = no_bursting;
for p=1:6
    summary.(char("time" + capital(PERIODS(p)))) = 100*period_min(p)/sum(period_min);
end
for p=1:6
    summary.(char("share" + capital(PERIODS(p)))) = median(share(:,p), 'omitnan');
end
role_fields = ["nCore" "nLeader" "nTail" "nReverberation" "nSmallEvent" "nBaseline"];
for r=1:6
    summary.(char(role_fields(r))) = sum(is_role(:,r));
    if ~events_ok
        summary.(char(role_fields(r))) = NaN;
    end
end
summary.leaderLead = median(onset_lead(is_role(:,2)));
if ~any(is_role(:,2))
    summary.leaderLead = NaN;
end
summary.orderConsistency = order_S;
summary.orderP = order_p;
summary.orderMatch = order_match;
summary.orderSpan = span;
summary.orderSpanFraction = span_fraction;
summary.orderLag = order_lag;
summary.orderJitter = order_jitter;
summary.smallWithCNE = sm.with;
summary.smallRecruited = sm.recruited;
summary.smallPosition = sm.position;
summary.smallAgreement = sm.A;
summary.smallAgreementP = sm.pA;
summary.smallConsistency = sm_own.S;
summary.smallConsistencyP = sm_own.pS;
summary.secondaryWithCNE = se.with;
summary.secondaryRecruited = se.recruited;
summary.secondaryPosition = se.position;
summary.secondaryAgreement = se.A;
summary.secondaryAgreementP = se.pA;
summary.secondaryConsistency = se.S;
summary.secondaryConsistencyP = se.pS;
summary.referenceAgreement = ref_A;
summary.referenceAgreementP = ref_pA;
summary.referenceMatch = ref_tau;

cluster = cell(nk,1);
for k=1:nk
    cluster{k} = [];
    if isstruct(spikeTimestamp) && isfield(spikeTimestamp, 'cluster')
        cluster{k} = [spikeTimestamp(cNE.members{k}).cluster];
    end
end
stability = nan(nk,1);
if isfield(cNE, 'stability') && numel(cNE.stability) == nk
    stability = cNE.stability(:);
end
cols = {(1:nk)', n_members, cellfun(@(m) m(:)', cNE.members(:), 'UniformOutput', false), cluster, stability, rate, chance_rate, ...
    excess_rate, chance_share, at_chance, part, part_in, part_out};
names = {'cNE', 'nMembers', 'members', 'memberCluster', 'stability', 'rate', 'chanceRate', 'excessRate', 'chanceShare', 'atChance', ...
    'memberParticipation', 'participationInLarge', 'participationOutside'};
for p=1:6
    cols{end+1} = share(:,p);   names{end+1} = char("share" + capital(PERIODS(p))); %#ok<AGROW>
end
for p=1:6
    cols{end+1} = rate_p(:,p);  names{end+1} = char("rate" + capital(PERIODS(p))); %#ok<AGROW>
end
cols = [cols, {large_cover, onset_cover, onset_lead, tail_cover, small_cover, secondary_cover, tmpl(:), rk(:), pos(:), roles}];
names = [names, {'largeCover', 'onsetCover', 'onsetLead', 'tailCover', 'smallCover', 'secondaryCover', 'time', 'rank', 'position', 'roles'}];
cNEs = table(cols{:}, 'VariableNames', names);

events = struct();
events.large = struct('start', L(:,1), 'last', L(:,2), 'peak', pk_all(iL), 'tailEnd', tail_end(iL), 'window', [w0(iL) w1(iL)], ...
    'time', XL, 'nBins', nbinL, 'match', tau, 'peth', peth, 'lags', lags, 'peakTime', median(pk_all(iL) - L(:,1)), ...
    'coreEnd', median(L(:,2) - L(:,1)), 'windowEnd', median(w1(iL) - L(:,1)));
events.small = struct('start', S0all(i0,1), 'last', S0all(i0,2), 'inTail', tail0, 'time', X0, 'nBins', nbin0);
events.secondary = struct('start', S2all(i2,1), 'last', S2all(i2,2), 'time', X2, 'nBins', nbin2);

report = struct('summary', summary, 'cNEs', cNEs, 'events', events, 'binPeriod', uint8(label), 'periods', PERIODS, 'parameters', parameters);
if nargout > 1
    fig = draw_report(report);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions

function [L, pk, S0, S2, Z] = network_events(D)
% [start last] of the large (level-1), small (level-0) and secondary events, peak of the large events, closed zones [start end tail_end]
% (events without a level are large)
L = zeros(0,2);     pk = zeros(0,1);    S0 = zeros(0,2);    S2 = zeros(0,2);    Z = zeros(0,3);
if isempty(D) || ~isfield(D, 'nBurst') || isempty(D.nBurst) || D.nBurst == 0 || ~isfield(D, 'IndivBurst') || isempty(D.IndivBurst)
    return
end
E = D.IndivBurst;
if isfield(E, 'level')
    lv = [E.level];
else
    lv = ones(1, numel(E));
end
L = reshape([[E(lv == 1).start]' [E(lv == 1).last]'], [], 2);    pk = reshape([E(lv == 1).peak], [], 1);    % (also when all the events are small)
S0 = reshape([[E(lv == 0).start]' [E(lv == 0).last]'], [], 2);
if isfield(D, 'Secondary') && isfield(D.Secondary, 'nBurst') && ~isempty(D.Secondary.nBurst) && D.Secondary.nBurst > 0
    S2 = [[D.Secondary.IndivBurst.start]' [D.Secondary.IndivBurst.last]'];
end
if isfield(D, 'closedZones') && ~isempty(D.closedZones)
    Z = D.closedZones;
end
end

function m = in_intervals(tb, S)
% bins whose midpoint is inside any of the intervals [S(q,1) S(q,2)] (closed); tb sorted
m = false(size(tb));
for q=1:size(S,1)
    a = find(tb >= S(q,1), 1, 'first');     b = find(tb <= S(q,2), 1, 'last');
    if ~isempty(a) && ~isempty(b) && a <= b
        m(a:b) = true;
    end
end
end

function [x, nbin] = event_times(act, tb, sel, t0, pc, min_bins)
% time of each cNE in the bins sel, from t0: percentile pc of the times of its active bins (fewer than min_bins: midpoint of the first
% and the median)
nk = size(act,2);   x = nan(1,nk);  nbin = zeros(1,nk);
for k=1:nk
    ta = tb(sel(act(sel,k))) - t0;
    nbin(k) = numel(ta);
    if isempty(ta)
        continue
    end
    if numel(ta) >= min_bins
        x(k) = prctile(ta, pc);
    else
        x(k) = (min(ta) + median(ta))/2;
    end
end
x = round(x, 9);                                                            % (to the nsec: equal times are ties, without rounding errors)
end

function [rk, pos] = template_rank(tmpl)
% rank of each cNE in the template (ties: mean rank) and position (0 = first, 1 = last); NaN for the cNEs without a time
rk = nan(size(tmpl));   pos = nan(size(tmpl));
ok = find(~isnan(tmpl));    v = tmpl(ok);
for i=1:numel(ok)
    rk(ok(i)) = sum(v < v(i)) + (sum(v == v(i)) + 1)/2;
end
if numel(ok) >= 2
    pos(ok) = (rk(ok) - 1)/(numel(ok) - 1);
end
end

function [part, part_in, part_out] = member_participation(cNE, spikeTimestamp, in_large)
% median share of the members of each cNE that fire in its active bins (all, inside and outside the large events)
nk = cNE.nCNE;  part = nan(nk,1);   part_in = nan(nk,1);    part_out = nan(nk,1);
if isempty(spikeTimestamp) || nk == 0
    return
end
if iscell(spikeTimestamp)
    trains = spikeTimestamp(:);
else
    trains = {spikeTimestamp.timeStamp}';
end
if numel(trains) ~= cNE.nUnits
    error("cNEAnalysis: spikeTimestamp has %d units, the cNE detection had %d", numel(trains), cNE.nUnits);
end
bs = cNE.parameters.BIN_SIZE;
edges = [double(cNE.binTime(:)) - bs/2; double(cNE.binTime(end)) + bs/2];  % bin edges, in msec
units = unique(cell2mat(cellfun(@(m) m(:), cNE.members(:), 'UniformOutput', false)));
F = false(numel(units), numel(cNE.binTime));                                % bins where each member fires
for u=1:numel(units)
    b = discretize(double(trains{units(u)}(:)), edges);
    F(u, b(~isnan(b))) = true;
end
for k=1:nk
    [~, rows] = ismember(cNE.members{k}, units);
    e = cNE.events{k}(:);
    if isempty(e)
        continue
    end
    f = mean(F(rows, e), 1)';
    part(k) = median(f);
    part_in(k) = median(f(in_large(e)));
    part_out(k) = median(f(~in_large(e)));
end
end

function out = event_type(X, tmpl, pos, events_ok, order_ok, own, enough_cne, min_events, n_perm, stream)
% small or secondary events: % with an active cNE, % of the cNEs recruited, mean position of the recruited cNEs, agreement with the
% template, own order
out = struct('with', NaN, 'recruited', NaN, 'position', NaN, 'A', NaN, 'pA', NaN, 'S', NaN, 'pS', NaN);
if ~events_ok || isempty(X)
    return
end
na = sum(~isnan(X), 2);
out.with = 100*mean(na >= 1);
if any(na >= 1)                                                             % (November 2026) % of the cNEs an event recruits, mean
    out.recruited = 100*mean(na(na >= 1))/size(X,2);                        % over the events with an active cNE
end
if order_ok && sum(na >= 1) >= min_events
    P = repmat(pos, size(X,1), 1);
    out.position = mean(P(~isnan(X) & ~isnan(P)));
end
if order_ok && sum(na >= 2) >= min_events
    [out.A, out.pA] = agreement(X, tmpl, n_perm, stream);
end
if own && enough_cne && sum(na >= 2) >= min_events
    [S, out.pS, S0] = consistency(X, n_perm, stream);
    out.S = scale_to_chance(S, S0);
end
end

function v = scale_to_chance(S, S0)
v = (S - S0)/(1 - S0);
if isnan(S) || S0 >= 1 - 1e-12
    v = NaN;
end
end

function [S, p, S0] = consistency(X, n_perm, stream)
% order consistency of the events (rows) and its null: times permuted among the active cNEs of each event (S and p are NaN when the null
% is always 1, e.g., when every pair of cNEs is active together in a single event)
S = consistency_stat(X);    p = NaN;    S0 = NaN;
if isnan(S)
    return
end
null = zeros(n_perm,1);
for it=1:n_perm
    null(it) = consistency_stat(permute_rows(X, stream));
end
p = (sum(null >= S - 1e-12) + 1)/(n_perm + 1);
S0 = mean(null);
if S0 >= 1 - 1e-12
    p = NaN;
end
end

function S = consistency_stat(X)
k = size(X,2);  num = 0;    den = 0;
for a=1:k-1
    for b=a+1:k
        both = ~isnan(X(:,a)) & ~isnan(X(:,b));     n = sum(both);
        if n == 0
            continue
        end
        pr = mean((X(both,a) < X(both,b)) + 0.5*(X(both,a) == X(both,b)));
        num = num + n*abs(2*pr - 1);    den = den + n;
    end
end
S = num/den;
if den == 0
    S = NaN;
end
end

function Y = permute_rows(X, stream)
Y = X;
for q=1:size(X,1)
    a = find(~isnan(X(q,:)));
    Y(q,a) = X(q, a(randperm(stream, numel(a))));
end
end

function [A, p] = agreement(X, tmpl, n_perm, stream)
% share of the pairs of active cNEs in the template order (pooled over the events) and its null
A = agreement_stat(X, tmpl);    p = NaN;
if isnan(A)
    return
end
cnt = 0;
for it=1:n_perm
    cnt = cnt + (agreement_stat(permute_rows(X, stream), tmpl) >= A - 1e-12);
end
p = (cnt + 1)/(n_perm + 1);
end

function A = agreement_stat(X, tmpl)
num = 0;    n = 0;
for q=1:size(X,1)
    a = find(~isnan(X(q,:)) & ~isnan(tmpl));
    for u=1:numel(a)-1
        for v=u+1:numel(a)
            i = a(u);   j = a(v);
            if tmpl(i) == tmpl(j)
                continue
            end
            num = num + (sign(X(q,i) - X(q,j)) == sign(tmpl(i) - tmpl(j))) + 0.5*(X(q,i) == X(q,j));
            n = n + 1;
        end
    end
end
A = num/n;
if n == 0
    A = NaN;
end
end

function [lag, jitter] = pair_lags(X, tmpl)
% timing of the order: for each pair of cNEs active together in at least 3 large events, the lag between them in each of these events
% (later minus earlier cNE of the template); median lag and median absolute deviation (jitter), medians over the pairs
[~, o] = sort(tmpl);    o = o(~isnan(tmpl(o)));     G = [];     J = [];
for a=1:numel(o)-1
    for b=a+1:numel(o)
        both = ~isnan(X(:,o(a))) & ~isnan(X(:,o(b)));
        if sum(both) >= 3
            d = X(both,o(b)) - X(both,o(a));
            G(end+1) = median(d);                                           %#ok<AGROW>
            J(end+1) = median(abs(d - median(d)));                          %#ok<AGROW>
        end
    end
end
[lag, jitter] = deal(NaN);
if ~isempty(G)
    lag = median(G);    jitter = median(J);
end
end

function tau = loo_tau(X)
% per event: Kendall tau between its order and the template of the other events (events with >= 3 active cNEs with a template time)
nE = size(X,1);     tau = nan(nE,1);
for q=1:nE
    a = find(~isnan(X(q,:)));
    if numel(a) < 3
        continue
    end
    tm = median(X(setdiff(1:nE, q), a), 1, 'omitnan');     ok = ~isnan(tm);
    if sum(ok) < 3
        continue
    end
    tau(q) = kendall_b(X(q, a(ok)), tm(ok));
end
end

function t = kendall_b(x, y)
% Kendall tau-b (ties adjusted)
n = numel(x);   s = 0;  tx = 0;     ty = 0;     n0 = n*(n-1)/2;
for i=1:n-1
    dx = sign(x(i+1:n) - x(i));     dy = sign(y(i+1:n) - y(i));
    s = s + sum(dx.*dy);    tx = tx + sum(dx == 0);     ty = ty + sum(dy == 0);
end
t = s/sqrt((n0 - tx)*(n0 - ty));
if (n0 - tx) == 0 || (n0 - ty) == 0
    t = NaN;
end
end

function s = capital(s)
s = upper(extractBefore(s, 2)) + extractAfter(s, 1);
end

function fig = draw_report(report)
% shares of the periods (left), activity around the large events (middle) and participation in each type of event (right); cNEs in
% the order of the large events (the ones without a time at the end)
T = report.cNEs;    S = report.summary;     E = report.events;  nk = height(T);
blues = interp1(linspace(0, 1, 8), [255 255 255; 205 226 251; 158 197 244; 109 167 236; 57 135 229; 37 106 191; 24 79 149; 13 54 107]/255, ...
    linspace(0, 1, 256));
t = T.time;     t(isnan(t)) = Inf;  [~, ord] = sort(t);
labels = arrayfun(@(k) sprintf('cNE %d (%d), %.1f/min%s', T.cNE(k), T.nMembers(k), T.rate(k), chance_text(T.chanceShare(k))), ord, ...
    'UniformOutput', false);
fig = figure('Color', 'w', 'Position', [50 50 1600 120 + 40*max(nk, 4)]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
head = sprintf('%d cNEs of %d units; %d large, %d small, %d secondary events', nk, S.nUnits, S.nLarge, S.nSmall, S.nSecondary);
if S.noBursting
    head = [head '; no network bursting'];
else
    head = [head sprintf('; %.0f%% of the large events with an active cNE', S.largeWithCNE)];
end
if S.notDescribed
    head = [head ' (the cNEs do not describe these events)'];
end
title(tl, head, 'FontSize', 11, 'Interpreter', 'none');
if nk == 0
    return
end

ax = nexttile(tl);      % shares of the periods
V = [T.shareOnset T.shareCore T.shareSecondary T.shareTail T.shareSmall T.shareQuiet];  V = V(ord,:);
imagesc(ax, V, [0 100]);    colormap(ax, blues);    hold(ax, 'on');
for i=1:nk
    for p=1:6
        if isnan(V(i,p))
            txt = '-';
        else
            txt = sprintf('%.0f', V(i,p));
        end
        text(ax, p, i, txt, 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.1 0.1 0.1] + 0.9*(V(i,p) > 60));
    end
end
ax.XTick = 1:6;     ax.XTickLabel = cellstr(report.periods);    ax.YTick = 1:nk;
ax.YTickLabel = arrayfun(@(i) sprintf('%s  %s', labels{i}, T.roles(ord(i))), (1:nk)', 'UniformOutput', false);
ax.TickLabelInterpreter = 'none';   ax.FontSize = 8.5;  ax.TickLength = [0 0];
title(ax, {'% of the activations of each cNE in each period', 'cNE (members), activations per minute (% at chance), roles'}, 'FontSize', 9.5, ...
    'FontWeight', 'normal');

ax = nexttile(tl);      % activity around the large events
if S.nLarge > 0 && nk > 0
    M = E.large.peth(:, ord)';  mx = max(M, [], 2);     M = M./max(mx, 0.25);   % rows below 0.25 are not stretched
    imagesc(ax, E.large.lags, 1:nk, M, [0 1]);  colormap(ax, blues);    hold(ax, 'on');
    yl = [0.5 nk + 0.5];
    plot(ax, [0 0], yl, '-', 'Color', [0.85 0.3 0.1]);
    plot(ax, [1 1]*E.large.peakTime, yl, ':', 'Color', [0.85 0.3 0.1], 'LineWidth', 1.2);
    plot(ax, [1 1]*E.large.coreEnd, yl, '--', 'Color', [0.85 0.3 0.1]);
    ax.YTick = 1:nk;    ax.YTickLabel = arrayfun(@(i) sprintf('cNE %d [%.2f]', T.cNE(ord(i)), mx(i)), (1:nk)', 'UniformOutput', false);
    ax.TickLabelInterpreter = 'none';   ax.FontSize = 8.5;  ax.TickDir = 'out';
    xlabel(ax, 'time from the start of the large events (s)');
    st = sprintf('order consistency %.2f (p %.4f), per-event match %.2f, lag %.0f ms (jitter %.0f ms)', S.orderConsistency, S.orderP, ...
        S.orderMatch, 1000*S.orderLag, 1000*S.orderJitter);
    if isnan(S.orderConsistency)
        st = 'order: the cNEs are rarely active in the same large events';
        if S.notDescribed
            st = 'order: not computed (the cNEs do not describe the large events)';
        elseif S.nCNE < report.parameters.MIN_CNE_ORDER || S.nLarge < report.parameters.MIN_LARGE_ORDER
            st = 'order: below the minimums (cNEs, large events)';
        end
    end
    title(ax, {'share of the large events where each cNE is active [maximum]', st}, 'FontSize', 9.5, 'FontWeight', 'normal');
else
    axis(ax, 'off');    title(ax, 'no large events', 'FontSize', 9.5, 'FontWeight', 'normal');
end

ax = nexttile(tl);      % participation in each type of event
Y = [T.largeCover T.smallCover T.secondaryCover];   Y = Y(ord,:);
if nk > 0 && any(~isnan(Y(:)))
    x = 1:nk;
    if nk == 1
        Y = [Y; nan(1,3)];  x = 1:2;
    end
    b = bar(ax, x, Y, 'grouped', 'EdgeColor', 'none');
    colors = [42 120 214; 235 104 52; 27 175 122]/255;
    for i=1:3
        b(i).FaceColor = colors(i,:);
    end
    ax.XTick = 1:nk;    ax.XTickLabel = arrayfun(@(k) sprintf('cNE %d', k), T.cNE(ord), 'UniformOutput', false);
    ax.XLim = [0.4 nk + 0.6];   ax.YLim = [0 105];
    ax.TickLabelInterpreter = 'none';   ax.FontSize = 8.5;  ax.YGrid = 'on';    ax.Box = 'off';
    ylabel(ax, '% of the events where the cNE is active');
    legend(ax, b, ["large" "small" "secondary"], 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
else
    axis(ax, 'off');    title(ax, 'no network events', 'FontSize', 9.5, 'FontWeight', 'normal');
end
end

function txt = chance_text(share)
% chance share for the labels of the figure
txt = '';
if ~isnan(share)
    txt = sprintf(' (%.0f%%)', share);
end
end
