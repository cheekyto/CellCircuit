function passed = test_cNEAnalysis(code_dir)

% Function-level test of "cNEAnalysis" on synthetic cNE activations and network events (November 2026). It checks that
%   1. every bin gets the right period (against a direct labelling of the bins; precedence secondary > core > small > tail > onset > quiet),
%      and the events outside the window, or small events that overlap a large event, are not counted
%   2. shares and rates of the periods of each cNE (against a direct count), onset coverage and lead
%   3. a planted order (the same in every large event) gives an order consistency of 1 (p < 0.01) and a per-event match of 1, the right
%      ranks and positions; a random order gives a consistency near 0 (p > 0.05); lag and jitter between the cNEs (direct count, and
%      planted lags with a known median and median absolute deviation)
%   4. small and secondary events: share with an active cNE, position of the recruited cNEs, agreement with the template, own order;
%      the large events against a reference order (REFERENCE_TIME): the same order, reversed, none
%   5. roles (core, leader, tail, reverberation, small-event, baseline) and member participation (planted spikes)
%   6. events the cNEs do not describe (blank roles and order), 1 cNE (no order), 2 cNEs (which one leads), no network bursting, only
%      small events; chance level of the activations (shares, flag, excess rate)
%   7. same seed = same result; cell array of spike trains = structure array; the figure is drawn
%
% Example (from a terminal):
%   matlab -batch "addpath('<repo>/tests'); test_cNEAnalysis('<repo>')"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dir:             folder with the code to be tested

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% passed:               =true if all the checks passed (each check is printed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(char(code_dir));
passed = true;
BIN = 15;   T_SEC = 600;    W = 0.5;
tb = ((1:T_SEC*1000/BIN)' - 0.5)*BIN/1000;                                  % midpoints of the bins, in sec
at = @(t) reshape(round(t/(BIN/1000) + 0.5), [], 1);                       % bins that contain the times t (sec), as a column

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Synthetic network events

s = 30 + 55*(0:9)';                                                         % 10 large events of 1 s, tails of 2 s
L = [s, s + 1];     Z = [s, s + 1, s + 3];
L_out = [-10 -9; 599.5 600.5];                                              % outside the window (partly or entirely)
S0 = [s(1:4) + 25, s(1:4) + 25.3; s(5) + 2.2, s(5) + 2.4];                  % 4 small events in the baseline, 1 in a tail
S0_overlap = [s(6) + 0.4, s(6) + 1.2];                                      % overlaps a large event: left out
S2 = [s(1:4) + 0.62, s(1:4) + 0.8];                                         % secondary events (inside the large events)
NB = network(L, L_out, S0, S0_overlap, S2, [Z; L_out, L_out(:,2) + 1.5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Synthetic cNEs (bins where each one is active), in the order 1 -> 2 -> 3 -> 5 -> 4 within the large events

ev = cell(5,1);
ev{1} = [at(s(1:6) - 0.2); at(s + 0.05); at(s(1:3) + 24.9)];               % leader (60% of the onsets), core; before 3 small events
ev{2} = at(s + [0.2 0.25 0.3 0.35]);                                        % core
ev{3} = [at(s + [0.5 0.55 0.6]); at(s(1:4) + 0.65)];                        % core and the secondary events
ev{4} = at(s + [0.9 1.5 1.6 2.0]);                                          % tail
ev{5} = [at(s(1:4) + 0.7); at(s(1:3) + 25.1); at(s + 40)];                  % secondary events, 3 small events, quiet baseline
members = {1:3, 4:7, 8:10, 11:12, 13:16};
cNE = cne_struct(tb, ev, members, 20, BIN, T_SEC);

trains = cell(20,1);                                                        % spikes: cNE 2, 2 of its 4 members in each activation; cNE 5,
for b = ev{2}(:)'                                                           % all 4 inside the large events, 1 outside
    trains{4}(end+1,1) = tb(b)*1000;    trains{5}(end+1,1) = tb(b)*1000;
end
in_large5 = any(tb(ev{5}) >= L(:,1)' & tb(ev{5}) <= L(:,2)', 2);
for b = ev{5}(:)'
    for u = 13:(13 + 3*any(tb(b) >= L(:,1) & tb(b) <= L(:,2)))
        trains{u}(end+1,1) = tb(b)*1000;
    end
end
spikeTimestamp = struct('timeStamp', trains, 'cluster', num2cell(100 + (1:20)'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1-2. Periods

[R, fig] = cNEAnalysis(cNE, NB, spikeTimestamp);
close(fig);
allL = [L; L_out];  allZ = [Z; L_out, L_out(:,2) + 1.5];
inS = @(I) any(tb >= I(:,1)' & tb <= I(:,2)', 2);
ref = 6*ones(numel(tb),1);
ref(inS([allL(:,1) - W, allL(:,1)]) & ~inS(allL)) = 1;
ref(inS(allZ(:,[2 3])) & ~inS(allL)) = 4;
ref(inS([S0; S0_overlap])) = 5;     ref(inS(allL)) = 2;     ref(inS(S2)) = 3;
passed = report("period of every bin (direct labelling)", isequal(double(R.binPeriod), ref)) && passed;
passed = report("events counted: 10 large, 5 small (1 left out), 4 secondary", R.summary.nLarge == 10 && R.summary.nSmall == 5 && ...
    R.summary.nSmallLeftOut == 1 && R.summary.nSecondary == 4) && passed;
passed = report("small event in a tail flagged", isequal(R.events.small.inTail, [false; false; false; false; true])) && passed;

tm = accumarray(ref, 1, [6 1])'*BIN/1000/60;
ok = true;
for k=1:5
    c = accumarray(ref(ev{k}), 1, [6 1])';
    sh = [R.cNEs.shareOnset(k) R.cNEs.shareCore(k) R.cNEs.shareSecondary(k) R.cNEs.shareTail(k) R.cNEs.shareSmall(k) R.cNEs.shareQuiet(k)];
    rt = [R.cNEs.rateOnset(k) R.cNEs.rateCore(k) R.cNEs.rateSecondary(k) R.cNEs.rateTail(k) R.cNEs.rateSmall(k) R.cNEs.rateQuiet(k)];
    ok = ok && max(abs(sh - 100*c/numel(ev{k}))) < 1e-9 && max(abs(rt - c./tm)) < 1e-9 && abs(R.cNEs.rate(k) - numel(ev{k})/10) < 1e-9;
end
passed = report("shares and rates of the periods (direct count)", ok) && passed;
lead = median(s(1:6) - tb(at(s(1:6) - 0.2)));
passed = report(sprintf("onsets: cNE 1 in 60%% of them, lead %.4f s", R.cNEs.onsetLead(1)), R.cNEs.onsetCover(1) == 60 && ...
    abs(R.cNEs.onsetLead(1) - lead) < 1e-9 && all(R.cNEs.onsetCover(2:5) == 0)) && passed;
passed = report("coverage: tails (cNE 4), large events, secondary events", R.cNEs.tailCover(4) == 100 && all(R.cNEs.largeCover(1:4) == 100) && ...
    R.cNEs.largeCover(5) == 40 && R.summary.largeWithCNE == 100 && isequal(R.cNEs.secondaryCover', [0 0 100 0 100])) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3-4. Order

S = R.summary;
passed = report(sprintf("planted order: consistency %.2f (p %.4f), per-event match %.2f", S.orderConsistency, S.orderP, S.orderMatch), ...
    abs(S.orderConsistency - 1) < 1e-12 && S.orderP < 0.01 && S.orderMatch == 1) && passed;
passed = report("ranks and positions (1 -> 2 -> 3 -> 5 -> 4)", isequal(R.cNEs.rank', [1 2 3 5 4]) && max(abs(R.cNEs.position' - [0 0.25 0.5 1 0.75])) < 1e-12) && passed;
span = R.cNEs.time(4) - R.cNEs.time(1);
passed = report("span of the order", abs(S.orderSpan - span) < 1e-12 && abs(S.orderSpanFraction - span/1) < 1e-12) && passed;
X = R.events.large.time;    [~, o] = sort(R.cNEs.time');    G = [];     J = [];
for pr = nchoosek(o, 2)'                                                    % pairs (earlier, later) in the template order
    d = X(:,pr(2)) - X(:,pr(1));    d = d(~isnan(d));
    if numel(d) >= 3
        G(end+1) = median(d);   J(end+1) = median(abs(d - median(d)));     %#ok<AGROW>
    end
end
passed = report(sprintf("timing of the order: lag %.3f s, jitter %.3f s (direct count over %d pairs)", S.orderLag, S.orderJitter, ...
    numel(G)), numel(G) == 10 && abs(S.orderLag - median(G)) < 1e-12 && abs(S.orderJitter - median(J)) < 1e-12) && passed;
k = [-4 -2 0 0 2 4 6 -6 1 -1]';                                             % planted lags of 20 + k bins: median 0.3 s, jitter 2 bins
Rj = cNEAnalysis(cne_struct(tb, {at(s + 0.05); at(s + 0.05) + 20 + k}, {1:3, 4:6}, 20, BIN, T_SEC), NB);
passed = report(sprintf("planted lags: lag %.3f s (0.300), jitter %.3f s (0.030)", Rj.summary.orderLag, Rj.summary.orderJitter), ...
    abs(Rj.summary.orderLag - 0.3) < 1e-9 && abs(Rj.summary.orderJitter - 0.03) < 1e-9) && passed;
passed = report(sprintf("small events: %.0f%% with a cNE, position %.3f, agreement %.2f (p %.3f), own order %.2f", S.smallWithCNE, ...
    S.smallPosition, S.smallAgreement, S.smallAgreementP, S.smallConsistency), S.smallWithCNE == 60 && ...
    abs(S.smallPosition - 0.375) < 1e-12 && S.smallAgreement == 1 && abs(S.smallConsistency - 1) < 1e-12) && passed;
ev_mix = ev;    ev_mix{5}(ev_mix{5} == at(s(3) + 25.1)) = at(s(3) + 24.8);   % cNE 5 before cNE 1 in 1 of the 3 small events: p 2/3
Rm = cNEAnalysis(cne_struct(tb, ev_mix, members, 20, BIN, T_SEC), NB);
passed = report(sprintf("small events, mixed own order: consistency %.3f (about -1/3, below chance; null mean from permutations)", ...
    Rm.summary.smallConsistency), abs(Rm.summary.smallConsistency + 1/3) < 0.05) && passed;
passed = report(sprintf("secondary events: %.0f%% with a cNE, position %.3f, agreement %.2f, own order %.2f", S.secondaryWithCNE, ...
    S.secondaryPosition, S.secondaryAgreement, S.secondaryConsistency), S.secondaryWithCNE == 100 && abs(S.secondaryPosition - 0.625) < 1e-12 && ...
    S.secondaryAgreement == 1 && abs(S.secondaryConsistency - 1) < 1e-12) && passed;
rec = @(X) 100*mean(nonzeros(sum(~isnan(X), 2)))/size(X,2);                % direct count: mean % of the cNEs an event recruits,
                                                                            % over the events with an active cNE
passed = report(sprintf("cNEs recruited per event with a cNE: small %.0f%% (2 of 5 in 3 events; 2 events without), secondary %.0f%%", ...
    S.smallRecruited, S.secondaryRecruited), abs(S.smallRecruited - 40) < 1e-9 && abs(S.secondaryRecruited - 40) < 1e-9 && ...
    abs(S.smallRecruited - rec(R.events.small.time)) < 1e-9 && abs(S.secondaryRecruited - rec(R.events.secondary.time)) < 1e-9) && passed;
ev_rec = ev;    ev_rec{2} = [ev{2}; at(s(1) + 25.2)];                       % cNE 2 also in the 1st small event: 3, 2, 2 of 5 cNEs
Rr = cNEAnalysis(cne_struct(tb, ev_rec, members, 20, BIN, T_SEC), NB);
ev_none = ev;   ev_none{1} = setdiff(ev{1}, at(s(1:3) + 24.9));    ev_none{5} = setdiff(ev{5}, at(s(1:3) + 25.1));   % no cNE in the small events
Rz = cNEAnalysis(cne_struct(tb, ev_none, members, 20, BIN, T_SEC), NB);
passed = report(sprintf("cNEs recruited: 3, 2, 2 of 5 -> %.1f%% (mean; median 40%%, mean over all events 28%%); no small event with a cNE: %.0f%% with, recruited %g", ...
    Rr.summary.smallRecruited, Rz.summary.smallWithCNE, Rz.summary.smallRecruited), abs(Rr.summary.smallRecruited - 140/3) < 1e-9 && ...
    Rz.summary.smallWithCNE == 0 && isnan(Rz.summary.smallRecruited)) && passed;
Ra = cNEAnalysis(cNE, NB, [], struct('REFERENCE_TIME', R.cNEs.time'));      % against a reference order: the same, then reversed
Rb = cNEAnalysis(cNE, NB, [], struct('REFERENCE_TIME', -R.cNEs.time'));
passed = report(sprintf("reference order: same %.2f (p %.4f, match %.0f), reversed %.2f (match %.0f); none: NaN", ...
    Ra.summary.referenceAgreement, Ra.summary.referenceAgreementP, Ra.summary.referenceMatch, Rb.summary.referenceAgreement, ...
    Rb.summary.referenceMatch), Ra.summary.referenceAgreement == 1 && Ra.summary.referenceAgreementP < 0.01 && Ra.summary.referenceMatch == 1 ...
    && Rb.summary.referenceAgreement == 0 && Rb.summary.referenceMatch == -1 && isnan(S.referenceAgreement) && isnan(S.referenceMatch)) && passed;

[Lr, evr] = random_order(80);                                              % 4 cNEs, one activation each at a random time in each event
Rr = cNEAnalysis(cne_struct(tb, evr, {1:3, 4:6, 7:9, 10:12}, 12, BIN, T_SEC), network(Lr, zeros(0,2), zeros(0,2), zeros(0,2), zeros(0,2), zeros(0,3)));
passed = report(sprintf("random order: consistency %.2f (p %.2f), per-event match %.2f", Rr.summary.orderConsistency, Rr.summary.orderP, ...
    Rr.summary.orderMatch), abs(Rr.summary.orderConsistency) < 0.3 && Rr.summary.orderP > 0.05) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Roles and member participation

expected = ["core, leader, small-event"; "core"; "core, reverberation"; "tail"; "reverberation, small-event, baseline"];
passed = report("roles: " + strjoin(R.cNEs.roles, " | "), isequal(R.cNEs.roles, expected)) && passed;
passed = report("number of cNEs per role", isequal([S.nCore S.nLeader S.nTail S.nReverberation S.nSmallEvent S.nBaseline], [3 1 1 2 2 1]) && ...
    abs(S.leaderLead - lead) < 1e-9) && passed;
passed = report("member participation (planted spikes)", R.cNEs.memberParticipation(2) == 0.5 && R.cNEs.participationInLarge(5) == 1 && ...
    R.cNEs.participationOutside(5) == 0.25 && R.cNEs.memberParticipation(5) == median(0.25 + 0.75*in_large5) && ...
    isnan(R.cNEs.participationOutside(2))) && passed;
passed = report("members and their cluster ids", isequal(R.cNEs.members{2}, 4:7) && isequal(R.cNEs.memberCluster{2}, 104:107)) && passed;
passed = report("summary of the cNEs", S.nCNE == 5 && S.nUnits == 20 && S.medianMembers == 3 && S.unitsInCNE == 80 && ...
    S.activationRate == median(cellfun(@numel, ev)/10)) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Events the cNEs do not describe, minimums, no network bursting

few = cellfun(@(e) e(tb(e) < s(4) - 1 | tb(e) > s(10) + 5), ev, 'UniformOutput', false);  % active in the first 3 large events only
Rf = cNEAnalysis(cne_struct(tb, few, members, 20, BIN, T_SEC), NB);
passed = report(sprintf("cNEs in %.0f%% of the large events: not described, roles and order blank", Rf.summary.largeWithCNE), ...
    Rf.summary.notDescribed && Rf.summary.largeWithCNE == 30 && isnan(Rf.summary.nCore) && all(Rf.cNEs.roles == "") && ...
    isnan(Rf.summary.orderConsistency) && all(isnan(Rf.cNEs.time)) && isnan(Rf.summary.smallWithCNE) && isnan(Rf.summary.smallRecruited) && ~isnan(Rf.summary.shareCore) && ...
    Rf.cNEs.largeCover(2) == 30) && passed;
R1 = cNEAnalysis(cne_struct(tb, ev(1), members(1), 20, BIN, T_SEC), NB);
passed = report("1 cNE: order measures blank, roles kept", isnan(R1.summary.orderConsistency) && isnan(R1.summary.smallPosition) && ...
    isnan(R1.summary.orderLag) && isnan(R1.summary.orderJitter) && isnan(R1.summary.smallConsistency) && ...
    isnan(R1.cNEs.rank) && ~isnan(R1.cNEs.time) && contains(R1.cNEs.roles(1), "leader")) && passed;
R2 = cNEAnalysis(cne_struct(tb, ev(1:2), members(1:2), 20, BIN, T_SEC), NB);
passed = report(sprintf("2 cNEs: which one leads, consistency %.2f (p %.4f), per-event match blank", R2.summary.orderConsistency, ...
    R2.summary.orderP), abs(R2.summary.orderConsistency - 1) < 1e-12 && R2.summary.orderP < 0.01 && isnan(R2.summary.orderMatch) && ...
    isequal(R2.cNEs.rank', [1 2]) && isequal(R2.cNEs.position', [0 1]) && R2.cNEs.roles(2) == "core") && passed;
R5 = cNEAnalysis(cNE, network(L(1:4,:), zeros(0,2), zeros(0,2), zeros(0,2), zeros(0,2), Z(1:4,:)));
passed = report("4 large events: order measures blank", isnan(R5.summary.orderConsistency) && R5.summary.nLarge == 4 && ~isnan(R5.summary.nCore)) && passed;
[Rn, fig] = cNEAnalysis(cNE, struct('nBurst', 0, 'closedZones', Z), spikeTimestamp);
close(fig);
passed = report("no network bursting: cNEs, rates, participation only", Rn.summary.noBursting && ~Rn.summary.notDescribed && ...
    Rn.summary.nLarge == 0 && Rn.summary.shareQuiet == 100 && isnan(Rn.summary.nLeader) && isnan(Rn.summary.orderConsistency) && ...
    all(isnan(Rn.cNEs.largeCover)) && Rn.summary.activationRate == S.activationRate && ~isnan(Rn.summary.memberParticipation)) && passed;
Rs0 = cNEAnalysis(cNE, network(zeros(0,2), zeros(0,2), S0, zeros(0,2), zeros(0,2), zeros(0,3)));     % only small events
passed = report(sprintf("only small events: no bursting, %d small events", Rs0.summary.nSmall), Rs0.summary.noBursting && ...
    Rs0.summary.nLarge == 0 && Rs0.summary.nSmall == 5 && isnan(Rs0.summary.smallConsistency)) && passed;
nev = cellfun(@numel, cNE.events);                                         % chance level: 0, 20, 50, 60 and 100% of the activations
cNEc = cNE;     cNEc.chanceEvents = nev.*[0 0.2 0.5 0.6 1];
Rc = cNEAnalysis(cNEc, NB);
mins = T_SEC/60;
passed = report(sprintf("chance level: shares %s%%, flagged %s, excess = rate - chance rate; without it: NaN, none flagged", ...
    mat2str(round(Rc.cNEs.chanceShare')), mat2str(Rc.cNEs.atChance')), max(abs(Rc.cNEs.chanceShare' - [0 20 50 60 100])) < 1e-9 && ...
    isequal(Rc.cNEs.atChance', [false false true true true]) && Rc.summary.nAtChance == 3 && ...
    max(abs(Rc.cNEs.excessRate' - nev.*[1 0.8 0.5 0.4 0]/mins)) < 1e-9 && all(isnan(R.cNEs.chanceShare)) && ~any(R.cNEs.atChance)) && passed;
R0 = cNEAnalysis(cne_struct(tb, {}, {}, 20, BIN, T_SEC), NB);
passed = report("no cNEs", R0.summary.nCNE == 0 && height(R0.cNEs) == 0 && R0.summary.unitsInCNE == 0) && passed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. Reproducibility

passed = report("same seed = same result; cell array of spikes = structure array", isequaln(R, cNEAnalysis(cNE, NB, spikeTimestamp)) && ...
    isequaln(R.cNEs.memberParticipation, cNEAnalysis(cNE, NB, trains).cNEs.memberParticipation)) && passed;
Rs = cNEAnalysis(cNE, NB, [], struct('SEED', 2));
passed = report("another seed: same measures except the p values", Rs.summary.orderConsistency == S.orderConsistency && ...
    isequaln(Rs.cNEs(:, setdiff(Rs.cNEs.Properties.VariableNames, {'memberParticipation', 'participationInLarge', 'participationOutside', 'memberCluster'})), ...
    R.cNEs(:, setdiff(R.cNEs.Properties.VariableNames, {'memberParticipation', 'participationInLarge', 'participationOutside', 'memberCluster'})))) && passed;

if passed
    fprintf("\ncNEAnalysis: all checks passed\n");
else
    fprintf("\ncNEAnalysis: SOME CHECKS FAILED\n");
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
function cNE = cne_struct(tb, ev, members, n_units, bin, t_sec)
% fields of the output of "cNEDetection" used by "cNEAnalysis"

nk = numel(ev);
cNE = struct('nCNE', nk, 'nUnits', n_units, 'binTime', tb*1000, 'active', true(n_units,1), 'stability', ones(1,nk));
cNE.events = cellfun(@(e) unique(e(:))', ev(:)', 'UniformOutput', false);
cNE.members = members(:)';
cNE.parameters = struct('BIN_SIZE', bin, 'tRange', [0 t_sec]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NB = network(L, L_out, S0, S0_overlap, S2, Z)
% output of "NetworkBurstDetection" (option 5) with these events

big = [L; L_out];   small = [S0; S0_overlap];
E = struct('start', num2cell([big(:,1); small(:,1)]), 'last', num2cell([big(:,2); small(:,2)]), ...
    'peak', num2cell([big(:,1) + 0.1; small(:,1) + 0.1]), 'level', num2cell([ones(size(big,1),1); zeros(size(small,1),1)]));
NB = struct('nBurst', numel(E), 'IndivBurst', E, 'closedZones', Z);
NB.Secondary = struct('nBurst', size(S2,1), 'IndivBurst', struct('start', num2cell(S2(:,1)), 'last', num2cell(S2(:,2))));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L, ev] = random_order(n)
% n large events of 1 s; 4 cNEs, each active in one bin at a random time of each event

stream = RandStream('mt19937ar', 'Seed', 7);
L = 20 + 7*(0:n-1)';    L = [L, L + 1];
ev = cell(4,1);
for k=1:4
    ev{k} = round((L(:,1) + 0.8*rand(stream, n, 1))/0.015 + 0.5);
end

end
