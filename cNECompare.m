function [comparison, fig] = cNECompare(reference, test, inputs)

% function to compare the coordinated neuronal ensembles (cNEs) of two conditions of the same recording, e.g., the baseline and a drug,
% or the two halves of the baseline (the change that time alone brings, to compare a drug with)
% (November 2026) Two complementary views:
%   1. tracking: the cNEs of the reference condition ("cNEDetection") are followed in the test condition: same weights and members, with
%      their activation computed on the test condition's own data ("cNEDetection" with TEMPLATE). Both conditions are described with
%      "cNEAnalysis" (each with its own network events; the test condition also against the order of the reference). Per cNE: the change
%      of its activation rate, of the firing rate of its members (if the activations fall more than the members' rates, the members fire
%      less together, not only less), of its member participation, of the share of the large events where it is active, of its roles and
%      of its position in the order
%   2. separate detection: the cNEs of the test condition are detected on their own ("cNEDetection") and matched one to one with the
%      reference cNEs by the similarity of their weights (cosine; greedy, from the most similar pair; a match needs at least
%      MATCH_SIMILARITY, the cut of the consensus of "cNEDetection"): reference cNEs found again, lost, and new cNEs
% Ratios are test/reference (NaN when the reference value is 0)
% Chance level (November 2026): the activations of a cNE include the ones the rules let through by chance ("cNEDetection", chanceEvents);
% the excess rate (activations - chance, per minute) is compared too, and the cNEs whose reference activations are mostly at chance
% ("cNEAnalysis", atChance) are flagged and left out of the summary of the excess

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   reference, test:    structures, one per condition, with the fields
%       spikeTimestamp          spike trains (as for "cNEDetection"); the same units, in the same order, in both conditions (checked with
%                               .cluster, when present)
%       tRange                  [start end] of the window of the condition, in sec
%       NetworkBurstDetected    (optional) network events of the condition ("NetworkBurstDetection", option 5); [] = none
%       name                    (optional) name of the condition (figure)
%       cNE                     (optional) output of "cNEDetection" for this window, to reuse it (otherwise it is detected here)
%   inputs:             (optional) structure that overwrites any of the PARAMETERS below, e.g., struct('ANALYSIS', struct('N_PERM', 500))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% comparison.summary        structure of scalars (shares in %):
%   .nCNE, .nCNETest            cNEs of the reference; cNEs detected in the test condition on their own
%   .nFound, .nLost, .nNew      reference cNEs matched by a test cNE, not matched; test cNEs not matched (2)
%   .foundShare                 % of the reference cNEs found again
%   .medianSimilarity           median over the reference cNEs of the similarity of the most similar test cNE (NaN without test cNEs)
%   .activeShare                % of the reference cNEs active at least once in the test condition (1)
%   .rateRatio                  median over the reference cNEs of the ratio of their activation rates
%   .nAtChance, .excessRatio    reference cNEs at chance; median ratio of the excess rates over the other reference cNEs (the excess of the
%                               test condition counted as 0 when negative)
%   .memberRateRatio            median over the reference cNEs of the ratio of the mean firing rate of their members
%   .participationChange        median change of the member participation (test - reference)
%   .largeCoverChange           median change of the share of the large events where the cNE is active (percentage points; NaN if a
%                               condition has no large events)
%   .rolesKept, .rolesLost, .rolesGained    number of roles (cNE and role) kept, lost and gained (NaN when the roles of a condition are not
%                               defined: no large events, or the cNEs do not describe them)
%   .referenceAgreement, .referenceAgreementP, .referenceMatch     the large events of the test condition against the order of the
%                               reference (see "cNEAnalysis", REFERENCE_TIME; computed only when the reference has the order measures)
%   .nLargeReference, .nLargeTest, .notDescribedReference, .notDescribedTest, .noBurstingReference, .noBurstingTest
% comparison.cNEs           table, one row per reference cNE:
%   cNE, nMembers, memberCluster, silentMembers (members without spikes in the test window), rateReference, rateTest, rateRatio
%   (activations per minute), chanceShareReference, chanceShareTest (% of the activations expected by chance), atChance (reference
%   activations mostly at chance), excessReference, excessTest, excessRatio (excess rates, per minute), memberRateReference,
%   memberRateTest, memberRateRatio (mean firing rate of the members, in Hz),
%   participationReference, participationTest, largeCoverReference, largeCoverTest, positionReference, positionTest, rolesReference,
%   rolesTest, match (test cNE matched, 0 = none), similarity (to that cNE; NaN = none), bestSimilarity (to the most similar test cNE)
% comparison.similarity     weight similarity of each reference cNE (rows) with each test cNE (columns)
% comparison.newCNEs        test cNEs not matched (index into comparison.detected.cNE)
% comparison.reference      .name, .tRange, .cNE, .report ("cNEAnalysis")
% comparison.tracked        .cNE (reference cNEs tracked in the test window), .report ("cNEAnalysis", with the reference order)
% comparison.detected       .cNE (separate detection in the test window), .report ("cNEAnalysis")
% comparison.test           .name, .tRange
% comparison.parameters     parameters used
% fig                       (only if requested) figure: activation rate and member firing rate of each cNE in both conditions, and the
%                           weight similarity with the cNEs detected in the test condition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MATCH_SIMILARITY = 0.5;          % separate detection: a test cNE matches a reference cNE with at least this weight similarity (cosine)
DETECTION = struct();            % inputs of "cNEDetection" (both conditions; e.g., struct('SEED', 1))
ANALYSIS = struct();             % inputs of "cNEAnalysis" (both conditions; e.g., struct('N_PERM', 500))
apply_parameters('cNECompare');                                             % values of the parameter set of "run_pipeline", if any
if nargin > 2 && ~isempty(inputs)
    override_inputs(inputs);
end
parameters = struct('MATCH_SIMILARITY', MATCH_SIMILARITY, 'DETECTION', DETECTION, 'ANALYSIS', ANALYSIS);

reference = complete(reference, "reference");
test = complete(test, "test");
trains_ref = as_trains(reference.spikeTimestamp);   trains_test = as_trains(test.spikeTimestamp);
if numel(trains_ref) ~= numel(trains_test)
    error("cNECompare: the reference has %d units, the test condition %d (the same units are needed)", numel(trains_ref), numel(trains_test));
end
if isstruct(reference.spikeTimestamp) && isstruct(test.spikeTimestamp) && isfield(reference.spikeTimestamp, 'cluster') && ...
        isfield(test.spikeTimestamp, 'cluster') && ~isequal([reference.spikeTimestamp.cluster], [test.spikeTimestamp.cluster])
    error("cNECompare: the units of the two conditions are not the same (.cluster differs)");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reference cNEs, tracked cNEs and separate detection

C0 = reference.cNE;
if isempty(C0)
    C0 = cNEDetection(reference.spikeTimestamp, reference.tRange, DETECTION);
end
R0 = cNEAnalysis(C0, reference.NetworkBurstDetected, reference.spikeTimestamp, ANALYSIS);

C1 = cNEDetection(test.spikeTimestamp, test.tRange, setfield(DETECTION, 'TEMPLATE', C0));                       %#ok<SFLD>
an = ANALYSIS;
if any(~isnan(R0.cNEs.rank))                                                % the reference has the order measures
    an.REFERENCE_TIME = R0.cNEs.time';
end
R1 = cNEAnalysis(C1, test.NetworkBurstDetected, test.spikeTimestamp, an);

C2 = test.cNE;
if isempty(C2)
    C2 = cNEDetection(test.spikeTimestamp, test.tRange, DETECTION);
end
R2 = cNEAnalysis(C2, test.NetworkBurstDetected, test.spikeTimestamp, ANALYSIS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matching of the separate detection (2)

nk = C0.nCNE;   nt = C2.nCNE;
sim = C0.weights'*C2.weights;                                               % cosine (unit-length weights), nk x nt
match = zeros(nk,1);    sim_match = nan(nk,1);
avail = sim;
for it=1:min(nk, nt)
    [v, idx] = max(avail(:));
    if isempty(v) || v < MATCH_SIMILARITY
        break
    end
    [i, j] = ind2sub(size(avail), idx);
    match(i) = j;   sim_match(i) = sim(i,j);
    avail(i,:) = -Inf;  avail(:,j) = -Inf;
end
best = nan(nk,1);
if nt > 0
    best = max(sim, [], 2);
end
new_cnes = setdiff(1:nt, match(match > 0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Each reference cNE in both conditions (1)

T0 = R0.cNEs;   T1 = R1.cNEs;
member_rate0 = cellfun(@(m) mean(C0.rate(m)), C0.members(:));
member_rate1 = cellfun(@(m) mean(C1.rate(m)), C0.members(:));
silent = cellfun(@(m) sum(C1.rate(m) == 0), C0.members(:));
if nk == 0
    [member_rate0, member_rate1, silent] = deal(zeros(0,1));
end
roles_kept = 0;     roles_lost = 0;     roles_gained = 0;
for k=1:nk
    a = split_roles(T0.roles(k));     b = split_roles(T1.roles(k));
    roles_kept = roles_kept + sum(ismember(a, b));
    roles_lost = roles_lost + sum(~ismember(a, b));
    roles_gained = roles_gained + sum(~ismember(b, a));
end
roles_ok = described(R0.summary) && described(R1.summary);
if ~roles_ok
    [roles_kept, roles_lost, roles_gained] = deal(NaN);
end

excess_ratio = ratio(max(T1.excessRate, 0), T0.excessRate);
excess_ratio(T0.excessRate <= 0) = NaN;
cNEs = table((1:nk)', T0.nMembers, T0.memberCluster, silent, T0.rate, T1.rate, ratio(T1.rate, T0.rate), T0.chanceShare, T1.chanceShare, ...
    T0.atChance, T0.excessRate, T1.excessRate, excess_ratio, member_rate0, member_rate1, ...
    ratio(member_rate1, member_rate0), T0.memberParticipation, T1.memberParticipation, T0.largeCover, T1.largeCover, T0.position, ...
    T1.position, T0.roles, T1.roles, match, sim_match, best, 'VariableNames', {'cNE', 'nMembers', 'memberCluster', 'silentMembers', ...
    'rateReference', 'rateTest', 'rateRatio', 'chanceShareReference', 'chanceShareTest', 'atChance', 'excessReference', 'excessTest', ...
    'excessRatio', 'memberRateReference', 'memberRateTest', 'memberRateRatio', 'participationReference', ...
    'participationTest', 'largeCoverReference', 'largeCoverTest', 'positionReference', 'positionTest', 'rolesReference', 'rolesTest', ...
    'match', 'similarity', 'bestSimilarity'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output

s = struct();
s.nCNE = nk;
s.nCNETest = nt;
s.nFound = sum(match > 0);
s.nLost = nk - s.nFound;
s.nNew = numel(new_cnes);
s.foundShare = 100*s.nFound/nk;
s.medianSimilarity = median(best);
s.activeShare = 100*mean(T1.rate > 0);
s.rateRatio = median(cNEs.rateRatio, 'omitnan');
s.nAtChance = sum(cNEs.atChance);
s.excessRatio = median(cNEs.excessRatio(~cNEs.atChance), 'omitnan');
s.memberRateRatio = median(cNEs.memberRateRatio, 'omitnan');
s.participationChange = median(T1.memberParticipation - T0.memberParticipation, 'omitnan');
s.largeCoverChange = median(T1.largeCover - T0.largeCover, 'omitnan');
s.rolesKept = roles_kept;
s.rolesLost = roles_lost;
s.rolesGained = roles_gained;
s.referenceAgreement = R1.summary.referenceAgreement;
s.referenceAgreementP = R1.summary.referenceAgreementP;
s.referenceMatch = R1.summary.referenceMatch;
s.nLargeReference = R0.summary.nLarge;
s.nLargeTest = R1.summary.nLarge;
s.notDescribedReference = R0.summary.notDescribed;
s.notDescribedTest = R1.summary.notDescribed;
s.noBurstingReference = R0.summary.noBursting;
s.noBurstingTest = R1.summary.noBursting;
if nk == 0
    [s.foundShare, s.medianSimilarity, s.activeShare, s.rateRatio, s.excessRatio, s.memberRateRatio, s.participationChange, ...
        s.largeCoverChange] = deal(NaN);
end

comparison = struct('summary', s, 'cNEs', cNEs, 'similarity', sim, 'newCNEs', new_cnes, ...
    'reference', struct('name', reference.name, 'tRange', reference.tRange, 'cNE', C0, 'report', R0), ...
    'tracked', struct('cNE', C1, 'report', R1), 'detected', struct('cNE', C2, 'report', R2), ...
    'test', struct('name', test.name, 'tRange', test.tRange), 'parameters', parameters);
if nargout > 1
    fig = draw_comparison(comparison);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions

function c = complete(c, label)
% optional fields of a condition
if ~isfield(c, 'spikeTimestamp') || ~isfield(c, 'tRange')
    error("cNECompare: the %s condition needs the fields spikeTimestamp and tRange", label);
end
defaults = struct('NetworkBurstDetected', [], 'name', label, 'cNE', []);
f = fieldnames(defaults);
for i=1:numel(f)
    if ~isfield(c, f{i})
        c.(f{i}) = defaults.(f{i});
    end
end
c.name = string(c.name);
end

function trains = as_trains(spikeTimestamp)
if iscell(spikeTimestamp)
    trains = spikeTimestamp(:);
else
    trains = {spikeTimestamp.timeStamp}';
end
end

function r = ratio(a, b)
r = a./b;
r(b == 0 | isnan(b)) = NaN;
end

function ok = described(S)
ok = ~S.noBursting && ~S.notDescribed;
end

function r = split_roles(txt)
r = strip(split(string(txt), ","));
r = r(strlength(r) > 0);
end

function fig = draw_comparison(cmp)
% activation rate (left) and member firing rate (middle) of each reference cNE in both conditions; weight similarity with the cNEs
% detected in the test condition (right)
T = cmp.cNEs;   S = cmp.summary;    nk = height(T);
names = [cmp.reference.name, cmp.test.name];
fig = figure('Color', 'w', 'Position', [50 50 1500 460]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
head = sprintf(['%s vs %s: %d reference cNEs (%d at chance), %d found again, %d lost, %d new; median ratio test/reference: activations %.2f, ' ...
    'excess over chance %.2f, member rates %.2f'], names(1), names(2), nk, S.nAtChance, S.nFound, S.nLost, S.nNew, S.rateRatio, S.excessRatio, ...
    S.memberRateRatio);
if ~isnan(S.referenceAgreement)
    head = sprintf('%s\nlarge events of %s in the order of %s: agreement %.2f (p %.4f)', head, names(2), names(1), S.referenceAgreement, ...
        S.referenceAgreementP);
end
title(tl, head, 'FontSize', 10.5, 'Interpreter', 'none');
if nk == 0
    return
end
grey = [0.55 0.55 0.55];    blue = [0.15 0.40 0.75];
pairs = {[T.rateReference T.rateTest], [T.memberRateReference T.memberRateTest]};
ylabels = {'activations per minute', 'mean firing rate of the members (Hz)'};
for p=1:2
    ax = nexttile(tl);  hold(ax, 'on');
    V = pairs{p};
    for k=1:nk
        c = grey;
        if T.match(k) > 0
            c = blue;
        end
        if T.atChance(k)                                                    % reference activations mostly at chance
            plot(ax, [1 2], V(k,:), '--o', 'Color', c, 'MarkerFaceColor', 'w', 'MarkerSize', 4);
        else
            plot(ax, [1 2], V(k,:), '-o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 4);
        end
        text(ax, 2.05, V(k,2), sprintf('%d', T.cNE(k)), 'FontSize', 7.5, 'Color', c);
    end
    ax.XLim = [0.7 2.3];    ax.XTick = [1 2];   ax.XTickLabel = cellstr(names);     ax.TickLabelInterpreter = 'none';
    ax.FontSize = 9;    ax.TickDir = 'out';     ylabel(ax, ylabels{p});
    ax.YLim(1) = 0;
    if p == 1
        title(ax, {'each reference cNE (blue: found again by the separate detection;', 'dashed: reference activations mostly at chance)'}, ...
            'FontSize', 9, 'FontWeight', 'normal');
    end
end
ax = nexttile(tl);
if S.nCNETest == 0
    axis(ax, 'off');
    text(ax, 0.5, 0.5, sprintf('no cNEs detected in %s', names(2)), 'HorizontalAlignment', 'center', 'Interpreter', 'none');
    return
end
imagesc(ax, max(cmp.similarity, 0), [0 1]);
colormap(ax, interp1([0 0.5 1], [1 1 1; 0.62 0.77 0.96; 0.05 0.21 0.42], linspace(0, 1, 256)));    hold(ax, 'on');
for k=1:nk
    if T.match(k) > 0
        rectangle(ax, 'Position', [T.match(k) - 0.5, k - 0.5, 1, 1], 'EdgeColor', [0.85 0.3 0.1], 'LineWidth', 1.5);
    end
end
ax.XTick = 1:S.nCNETest;    ax.YTick = 1:nk;    ax.FontSize = 9;    ax.TickLength = [0 0];
xlabel(ax, sprintf('cNEs detected in %s', names(2)), 'Interpreter', 'none');
ylabel(ax, sprintf('cNEs of %s', names(1)), 'Interpreter', 'none');
title(ax, sprintf('weight similarity (matched: >= %.2f, outlined)', cmp.parameters.MATCH_SIMILARITY), 'FontSize', 9, 'FontWeight', 'normal');
colorbar(ax);
end
