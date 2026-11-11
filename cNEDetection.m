function cNE = cNEDetection(spikeTimestamp, tRange, inputs)

% function to detect coordinated neuronal ensembles (cNEs): groups of units that repeatedly fire together within the same time bin
% (October 2026) Reimplementation of the method of See et al. 2018 (eLife 7:e35587), which is based on Lopes-dos-Santos et al. 2013
% (J Neurosci Methods 220:149-166); the authors did not share their code, so this follows the description in their Methods:
%   1. spike counts in bins of BIN_SIZE msec, z-scored for each unit
%   2. number of cNEs = number of eigenvalues of the correlation matrix above the Marcenko-Pastur upper bound, (1+sqrt(N/B))^2
%      (N units, B bins)
%   3. fastICA on the projection of the z-scored counts onto the significant principal components -> IC weight of each unit in each cNE
%   4. cNE activity = projection of the IC onto the z-scored counts, R(b) = z(b)'*P*z(b), with P = w*w' and zero diagonal
%   5. null distributions from N_SHUFFLE matrices in which the binned spike train of each unit is circularly shifted independently:
%      - membership: PCA/ICA of each shuffled matrix (same number of components as the real data); a unit is a member of a cNE when its
%        weight is above mean + MEMBER_SD*sd of all the shuffled weights (and a negative member when it is below mean - MEMBER_SD*sd)
%      - activity: the ICs of the real data projected onto the shuffled matrices; a cNE is active in the bins above the
%        ACTIVITY_PERCENTILE percentile of that null activity (one threshold per cNE)
% Details not given in the paper, taken from Lopes-dos-Santos et al. 2013 (See et al. used their MATLAB code): each IC is scaled to unit
% length and its sign is set so that its largest absolute weight is positive (also the ICs of the shuffled matrices, before pooling).
% The units without spikes in the window are left out (weight = 0). The binned spike counts are kept as a sparse matrix and the z-scoring
% is done algebraically, so that long recordings with many units fit in memory; the results are the same as with the dense matrices
%
% (November 2026) Robust definition, also for large recordings (default parameters). Tested on 17 recording-conditions (9 to 1059 units)
% by comparing the cNEs of the odd and of the even 30-s blocks of the window: 86% of the cNEs were found in both halves and 68% were
% nearly identical (weight similarity >= 0.9), against 63% and 23% with the steps above; on the largest recordings, fastICA no longer
% fails to converge and every run takes seconds instead of hours:
%   a. only the units that fire at least MIN_RATE Hz in the window are analyzed: when a unit that rarely fires does fire, its z-scored
%      count is very large, so a few coincident spikes (e.g., in network bursts) made spurious components, more of them with less data
%   b. number of cNEs (step 2) = number of eigenvalues above the EIG_PERCENTILE percentile of the LARGEST eigenvalue of N_SHIFT_EIG
%      circularly shifted matrices (N_CNE_RULE = "shift"); the Marcenko-Pastur bound (N_CNE_RULE = "MP") let through 3-7 times more
%      components, and the number changed between halves of the same recording. Known limitation: a group of units that share slow
%      rate changes (e.g., seconds-long high-rate states) raises this bound for every component, which can hide a weaker ensemble
%      (synthetic example in "test_cNEDetection")
%   c. consensus: steps 1-3 and the membership null (step 5) are run on N_RESAMPLE resamples of the window, each made of a fraction
%      RESAMPLE_FRACTION of its RESAMPLE_BLOCK-sec blocks (drawn without replacement, joined in time order). The components with at least
%      MIN_MEMBERS members of all the resamples are clustered by the similarity of their weights (cosine; average linkage, cut at
%      CLUSTER_SIMILARITY); a cNE is a cluster found in at least MIN_RESAMPLES resamples, represented by its centrotype (its component
%      most similar to the others), with the weights and members of that component; cNE.stability = fraction of the resamples where it
%      was found. N_RESAMPLE = 0: a single run on the whole window (steps 1-5), and the cNEs with fewer than MIN_MEMBERS members are dropped
%   d. activation, on the whole window with the weights of the cNEs: a cNE is active in a bin when (1) its activity (step 4) is above the
%      ACTIVITY_PERCENTILE percentile of its null activity (See et al.) AND (2) more of its members fire in the bin (at least one spike)
%      than expected by chance, with at least MIN_FIRING members firing (MEMBER_COUNT_RULE = false: rule (1) only). Rule (1) alone was
%      crossed in bins where only 0 or 1 member fired: the activity squares the projection, so many non-members with small negative
%      weights firing together also raise it
%   e. chance for the rule (2) (COUNT_NULL = "jitter"): interval jitter. The bins are grouped in windows of JITTER_WINDOW sec, and each
%      spike of a member is moved to a random bin of its window, so the rate changes slower than the window are kept; a bin passes when
%      the probability of at least the observed number of members firing is below 1 - COUNT_PERCENTILE/100 (computed exactly, no
%      surrogates). The first version of the rule compared the count with circularly shifted trains (COUNT_NULL = "shift"), which also
%      removes the slow rate changes: cNEs whose members fire at high rates together for seconds (e.g., after network bursts) were
%      "active" in most bins of those periods. The window sets which rate changes count as a shared state (slower than the window) and
%      which as coordination (faster); network bursts count as coordination with any window of 2-5 sec. The window favors the activations
%      inside network bursts, which the cNEs are meant to describe: with 5-s windows (default), 35 of the 37 cNEs of the recordings with
%      large network events stayed active in as many of them as with the first version, and the activations of the cNEs above more than
%      1 s from any large event fell by 62-94%. 2-s windows removed more of those (87-98%), but the small cNEs (3-7 members) lost up to
%      all of their large events; windows of 0.5 sec or less kept only 1-3 bins per event
%   f. (November 2026) chance level of the activations: both rules are applied, with the thresholds of the real data, to N_CHANCE
%      surrogates in which every spike of every unit analyzed is moved to a random bin of its jitter window (the null of rule (2));
%      cNE.chanceEvents(k) = mean number of active bins of the cNE k on the surrogates, i.e., the activations expected without any
%      coordination faster than the jitter windows. The rules test each bin at 0.1%, so a cNE whose members rarely fire together can have
%      a large part of its activations at chance level: on 22 recordings (43 cNEs), a median of 1-5% of the baseline activations, up to 38%,
%      and 80% for one weak cNE (0.7 activations per minute); see "cNEAnalysis" (chance share, excess rate, flag)
% The steps of See et al. alone: struct('MIN_RATE', 0, 'N_CNE_RULE', "MP", 'N_RESAMPLE', 0, 'MIN_MEMBERS', 0, 'MEMBER_COUNT_RULE', false,
% 'BIN_SIZE', 10)
%
% (November 2026) Tracking (TEMPLATE = the output of "cNEDetection" for another window of the same units, e.g., the baseline before a
% drug): no new cNEs are detected; the cNEs of the template keep their weights and members, and their activation (step d) is computed on
% this window, with the activation settings of the template (bin size, rules (1) and (2), their percentiles and nulls) applied to the
% counts of this window (z-scored with this window's mean and sd, and nulls built from this window's spike trains). The units analyzed are
% the ones analyzed in the template that fire in this window (also below MIN_RATE: a drop of the rate is part of a change); a member that
% is silent here cannot fire, which is the same as leaving it out of rule (2). Used by "cNECompare"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   spikeTimestamp:     structure array with the field .timeStamp (spike times of each unit, in msec), e.g., the output of
%                       "ShapeSpikeTrainData" or "Reshape_spikeTimestamp"; or a cell array with the spike times of each unit, in msec
%   tRange:             [start end] of the analysis window, in sec (an incomplete last bin is dropped; with the consensus, an incomplete
%                       last block is left out of the resamples)
%   inputs:             (optional) structure that overwrites any of the PARAMETERS below, e.g., struct('N_RESAMPLE', 0, 'SEED', 1);
%                       struct('TEMPLATE', cNE_baseline): tracking of the cNEs of cNE_baseline in this window (see above)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (tracking: nCNE, weights, members, negativeMembers, memberThreshold, stability, similarity and nullWeight are those of the template;
%  eigenvalues are those of this window, lambdaMax is NaN, icaConverged is NaN and resamples is empty; the rest is computed on this window)

% cNE.nCNE                      number of cNEs detected
% cNE.nUnits                    number of units (all of them, also the ones left out)
% cNE.rate(i)                   firing rate of the unit i in the window, in Hz
% cNE.active(i)                 =true if the unit i is analyzed (spikes in the window and a rate of at least MIN_RATE Hz)
% cNE.binTime(b)                midpoint of the bin b, in msec
% cNE.eigenvalues               eigenvalues of the correlation matrix of the analyzed units (whole window), in descending order
% cNE.lambdaMax                 bound for the eigenvalues of the single run (the eigenvalues above it are the cNEs); NaN with the
%                               consensus (the bound of each resample is in cNE.resamples)
% cNE.lambdaMP                  Marcenko-Pastur bound for the whole window (for reference)
% cNE.weights(i,k)              IC weight of the unit i in the cNE k (unit length per cNE; =0 for the units not analyzed)
% cNE.memberThreshold(k,:)      [lower upper] thresholds for the IC weights of the cNE k (null mean -/+ MEMBER_SD*sd of its run)
% cNE.members{k}                units with a weight above the upper threshold (index into spikeTimestamp), sorted by weight
% cNE.negativeMembers{k}        units with a weight below the lower threshold (they tend to be silent during the cNE events)
% cNE.stability(k)              fraction of the resamples in which the cNE k was found (1 for a single run)
% cNE.similarity(k)             mean weight similarity (cosine) between the components of its cluster (NaN for a single run, or a cluster
%                               of one component)
% cNE.activity(k,b)             activity of the cNE k in the bin b (whole window)
% cNE.activityThreshold(k)      threshold of the activation rule (1) for the cNE k (ACTIVITY_PERCENTILE of the null activity)
% cNE.membersNeeded(k)          ("shift") members of the cNE k that must fire in a bin for the activation rule (2) (NaN with "jitter",
%                               where it changes from bin to bin, or if MEMBER_COUNT_RULE = false)
% cNE.activityEvents{k}         bins that pass the rule (1)
% cNE.memberCountEvents{k}      bins that pass the rule (2) (empty if MEMBER_COUNT_RULE = false)
% cNE.events{k}                 bins in which the cNE k is active (both rules; rule (1) only if MEMBER_COUNT_RULE = false)
% cNE.eventTimes{k}             midpoint of those bins, in msec
% cNE.chanceEvents(k)           (November 2026) active bins of the cNE k expected by chance: mean over N_CHANCE interval-jitter surrogates
%                               (step f; NaN if N_CHANCE = 0)
% cNE.nullWeight                mean and standard deviation of the IC weights of the shuffled matrices of the run of each cNE (.mean(k), .sd(k))
% cNE.icaConverged              fraction of the fastICA runs (real data + shuffles, of every resample) that converged
% cNE.resamples(r)              (consensus) resample r: .blocks (blocks of RESAMPLE_BLOCK sec of the window), .nCNE, .lambdaMax,
%                               .weights(i,c) and .members{c} of all its components (also the ones with fewer than MIN_MEMBERS members),
%                               .icaConverged
% cNE.parameters                parameters used (and the analysis window)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BIN_SIZE = 15;                   % bin size, in msec (See et al.: 10 ms)
MIN_RATE = 0.5;                  % minimum firing rate in the window, in Hz, for a unit to be analyzed (0: every unit with spikes, as See et al.)
N_CNE_RULE = "shift";            % number of cNEs: "shift" = bound from circularly shifted data; "MP" = Marcenko-Pastur bound (See et al.)
N_SHIFT_EIG = 100;               % ("shift") number of circularly shifted matrices for the bound
EIG_PERCENTILE = 99;             % ("shift") the bound is this percentile of the largest eigenvalue of the shifted matrices
N_RESAMPLE = 10;                 % consensus: number of resamples of the window (0: a single run on the whole window, as See et al.)
RESAMPLE_BLOCK = 30;             % consensus: length of the blocks, in sec
RESAMPLE_FRACTION = 0.8;         % consensus: fraction of the blocks in each resample
MIN_RESAMPLES = 7;               % consensus: a cNE must be found in at least this number of resamples
CLUSTER_SIMILARITY = 0.5;        % consensus: the clusters of components are cut at this weight similarity (cosine)
MIN_MEMBERS = 3;                 % the cNEs with fewer members are dropped (with the consensus: the components, before the clustering)
N_SHUFFLE = 100;                 % number of circularly shifted matrices for the null distributions (See et al.: 100)
MEMBER_SD = 1.5;                 % membership threshold: null mean +/- MEMBER_SD standard deviations of the IC weights (See et al.: 1.5)
ACTIVITY_PERCENTILE = 99.9;      % activation rule (1): percentile of the null cNE activity (See et al.: 99.9)
MEMBER_COUNT_RULE = true;        % activation rule (2), on top of the rule (1): enough members fire in the bin (false: rule (1) only, as See et al.)
COUNT_NULL = "jitter";           % rule (2), chance level of the number of members firing in a bin: "jitter" = interval jitter of each member's
                                 % spikes within windows of JITTER_WINDOW sec; "shift" = circular shifts of each member's train (first version)
JITTER_WINDOW = 5;               % rule (2), "jitter": length of the jitter windows, in sec (rounded to whole bins, from the start of the window)
COUNT_PERCENTILE = 99.9;         % rule (2): "jitter": a bin passes when P(count >= observed) < 1 - COUNT_PERCENTILE/100; "shift": percentile
                                 % of the counts on the surrogates
N_SHIFT_COUNT = 200;             % rule (2), "shift": number of surrogates (the binned train of each member circularly shifted independently)
MIN_FIRING = 2;                  % rule (2): at least this number of members must fire in the bin
N_CHANCE = 20;                   % chance level of the activations: number of interval-jitter surrogates (0: not computed)
ICA_CONTRAST = "logcosh";        % fastICA nonlinearity: "logcosh" (negentropy approximation), "pow3" (kurtosis) or "gauss"
ICA_MAX_ITER = 1000;             % maximum number of fastICA iterations
ICA_TOL = 1e-6;                  % fastICA convergence: 1 - |cos(angle)| between consecutive estimates of every component
SEED = 0;                        % seed for the random number generators (shifts, resamples and fastICA initialization); same seed + same inputs = same result
TEMPLATE = [];                   % tracking: output of "cNEDetection" for another window of the same units, whose cNEs are followed in this
                                 % window instead of detecting new ones ([] = detection)
apply_parameters('cNEDetection');                                           % values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")
if nargin > 2 && ~isempty(inputs)
    override_inputs(inputs);
end
tracking = ~isempty(TEMPLATE);
if tracking                                                                 % the activation settings of the template
    tp = TEMPLATE.parameters;
    BIN_SIZE = tp.BIN_SIZE;     N_SHUFFLE = tp.N_SHUFFLE;   ACTIVITY_PERCENTILE = tp.ACTIVITY_PERCENTILE;
    MEMBER_COUNT_RULE = tp.MEMBER_COUNT_RULE;   COUNT_NULL = tp.COUNT_NULL;     JITTER_WINDOW = tp.JITTER_WINDOW;
    COUNT_PERCENTILE = tp.COUNT_PERCENTILE;     N_SHIFT_COUNT = tp.N_SHIFT_COUNT;   MIN_FIRING = tp.MIN_FIRING;
end
N_CNE_RULE = string(N_CNE_RULE);
if ~ismember(N_CNE_RULE, ["shift" "MP"])
    error("cNEDetection: unknown N_CNE_RULE '%s' (use ""shift"" or ""MP"")", N_CNE_RULE);
end
COUNT_NULL = string(COUNT_NULL);
if ~ismember(COUNT_NULL, ["jitter" "shift"])
    error("cNEDetection: unknown COUNT_NULL '%s' (use ""jitter"" or ""shift"")", COUNT_NULL);
end

K_MSEC_SEC = 1/1000;                                                        % factor to rescale from msec to sec
jitter_bins = max(round(JITTER_WINDOW/(BIN_SIZE*K_MSEC_SEC)), 1);           % length of the jitter windows, in bins

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Binned spike counts (sparse matrix, bins x units) and units analyzed

if iscell(spikeTimestamp)
    trains = spikeTimestamp(:);
else
    trains = {spikeTimestamp.timeStamp}';
end
n_units = numel(trains);
edges = tRange(1)/K_MSEC_SEC : BIN_SIZE : tRange(2)/K_MSEC_SEC;             % bin edges, in msec
n_bin = numel(edges) - 1;                                                   % total number of bins
if n_bin < 2
    error("cNEDetection: the analysis window must contain at least 2 bins of %g msec", BIN_SIZE);
end
bin_time = edges(1:end-1)' + BIN_SIZE/2;                                    % midpoint of each bin, in msec

spike_bin = cell(n_units,1);    spike_unit = cell(n_units,1);
for i=1:n_units                                                             % bin of every spike of each unit
    b = discretize(double(trains{i}(:)), edges);
    spike_bin{i} = b(~isnan(b));                                            % (spikes outside the window are ignored)
    spike_unit{i} = i*ones(numel(spike_bin{i}),1);
end
counts = sparse(vertcat(spike_bin{:}), vertcat(spike_unit{:}), 1, n_bin, n_units);   % spike count of each unit in each bin (repeated entries are added)

rate = full(sum(counts,1))'/(n_bin*BIN_SIZE*K_MSEC_SEC);                    % firing rate of each unit in the window (complete bins), in Hz
[mu, sd] = count_stats(counts);                                             % mean and standard deviation of the count of each unit
active = sd > 0 & rate >= MIN_RATE;                                         % units analyzed: spikes in the window, at least MIN_RATE Hz
if tracking                                                                 % the units of the template that fire in this window
    if numel(TEMPLATE.active) ~= n_units
        error("cNEDetection: the template has %d units, spikeTimestamp has %d (tracking needs the same units, in the same order)", ...
              numel(TEMPLATE.active), n_units);
    end
    active = TEMPLATE.active(:) & sd > 0;
end
S = counts(:,active);   mu = mu(active);    sd = sd(active);
n_active = sum(active);
active_idx = find(active);
if ~tracking && n_active >= n_bin
    error("cNEDetection: the number of bins (%d) must be larger than the number of units analyzed (%d); use a longer window", n_bin, n_active);
end

stream = RandStream("twister", "Seed", SEED);                               % private random number generator (reproducible)
prm = struct('N_CNE_RULE', N_CNE_RULE, 'N_SHIFT_EIG', N_SHIFT_EIG, 'EIG_PERCENTILE', EIG_PERCENTILE, 'N_SHUFFLE', N_SHUFFLE, ...
             'MEMBER_SD', MEMBER_SD, 'ACTIVITY_PERCENTILE', ACTIVITY_PERCENTILE, 'ICA_CONTRAST', ICA_CONTRAST, 'ICA_MAX_ITER', ICA_MAX_ITER, ...
             'ICA_TOL', ICA_TOL);

cNE.nUnits = n_units;
cNE.rate = rate;
cNE.active = active;
cNE.binTime = bin_time;
cNE.lambdaMP = (1 + sqrt(n_active/n_bin))^2;
cNE.parameters = struct('BIN_SIZE', BIN_SIZE, 'MIN_RATE', MIN_RATE, 'N_CNE_RULE', N_CNE_RULE, 'N_SHIFT_EIG', N_SHIFT_EIG, ...
                        'EIG_PERCENTILE', EIG_PERCENTILE, 'N_RESAMPLE', N_RESAMPLE, 'RESAMPLE_BLOCK', RESAMPLE_BLOCK, ...
                        'RESAMPLE_FRACTION', RESAMPLE_FRACTION, 'MIN_RESAMPLES', MIN_RESAMPLES, 'CLUSTER_SIMILARITY', CLUSTER_SIMILARITY, ...
                        'MIN_MEMBERS', MIN_MEMBERS, 'N_SHUFFLE', N_SHUFFLE, 'MEMBER_SD', MEMBER_SD, 'ACTIVITY_PERCENTILE', ACTIVITY_PERCENTILE, ...
                        'MEMBER_COUNT_RULE', MEMBER_COUNT_RULE, 'COUNT_NULL', COUNT_NULL, 'JITTER_WINDOW', JITTER_WINDOW, ...
                        'COUNT_PERCENTILE', COUNT_PERCENTILE, 'N_SHIFT_COUNT', N_SHIFT_COUNT, ...
                        'MIN_FIRING', MIN_FIRING, 'N_CHANCE', N_CHANCE, 'ICA_CONTRAST', ICA_CONTRAST, 'ICA_MAX_ITER', ICA_MAX_ITER, ...
                        'ICA_TOL', ICA_TOL, 'SEED', SEED, 'tRange', tRange, 'TEMPLATE_RANGE', []);

if tracking
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Tracking: the cNEs of the template, on this window

    cNE.parameters.TEMPLATE_RANGE = TEMPLATE.parameters.tRange;
    cNE.eigenvalues = sort(eig(zcorr(S, mu, sd)), 'descend');
    cNE.lambdaMax = NaN;
    V = TEMPLATE.weights(active,:);
    to_col = zeros(n_units,1);  to_col(active_idx) = 1:n_active;           % column of S of each unit analyzed
    members = cellfun(@(m) nonzeros(to_col(m)), TEMPLATE.members, 'UniformOutput', false);             % (members silent here left out)
    negative_members = cellfun(@(m) nonzeros(to_col(m)), TEMPLATE.negativeMembers, 'UniformOutput', false);
    member_threshold = TEMPLATE.memberThreshold;
    null_mean = TEMPLATE.nullWeight.mean;   null_sd = TEMPLATE.nullWeight.sd;
    stability = TEMPLATE.stability;     similarity = TEMPLATE.similarity;
    R = cne_activity(S, mu, sd, V);
    activity_threshold = zeros(1,0);
    if size(V,2) > 0
        activity_threshold = activity_null(S, mu, sd, V, stream, N_SHUFFLE, ACTIVITY_PERCENTILE);
    end
    n_converged = NaN;  n_runs = 0;
    resamples = struct('blocks', {}, 'nCNE', {}, 'lambdaMax', {}, 'weights', {}, 'members', {}, 'icaConverged', {});
elseif n_active < 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% (November 2026) Fewer than 2 units analyzed (e.g., a silenced network): no correlations, no cNEs

    cNE.eigenvalues = ones(n_active,1);
    cNE.lambdaMax = NaN;
    V = zeros(n_active,0);  members = {};   negative_members = {};
    member_threshold = zeros(0,2);  null_mean = zeros(1,0);     null_sd = zeros(1,0);
    stability = zeros(1,0);     similarity = zeros(1,0);
    R = zeros(n_bin,0);     activity_threshold = zeros(1,0);
    n_converged = 0;    n_runs = 0;
    resamples = struct('blocks', {}, 'nCNE', {}, 'lambdaMax', {}, 'weights', {}, 'members', {}, 'icaConverged', {});
elseif N_RESAMPLE == 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Single run on the whole window (steps 1-5)

    run = single_run(S, mu, sd, stream, prm, true);
    keep = cellfun(@numel, run.members) >= MIN_MEMBERS;
    cNE.eigenvalues = run.eigenvalues;
    cNE.lambdaMax = run.lambdaMax;
    V = run.V(:,keep);
    members = run.members(keep);     negative_members = run.negativeMembers(keep);
    member_threshold = repmat(run.memberThreshold, sum(keep), 1);
    null_mean = repmat(run.nullWeight.mean, 1, sum(keep));      null_sd = repmat(run.nullWeight.sd, 1, sum(keep));
    stability = ones(1, sum(keep));     similarity = nan(1, sum(keep));
    R = run.R(:,keep);
    activity_threshold = run.activityThreshold(keep);
    n_converged = run.nConverged;    n_runs = run.nRuns;
    resamples = struct('blocks', {}, 'nCNE', {}, 'lambdaMax', {}, 'weights', {}, 'members', {}, 'icaConverged', {});
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Consensus of N_RESAMPLE resamples of the window

    cNE.eigenvalues = sort(eig(zcorr(S, mu, sd)), 'descend');
    cNE.lambdaMax = NaN;
    block_bins = round(RESAMPLE_BLOCK/(BIN_SIZE*K_MSEC_SEC));                % length of a block, in bins
    n_block = floor(n_bin/block_bins);                                      % complete blocks in the window
    n_pick = round(RESAMPLE_FRACTION*n_block);                              % blocks in each resample
    if n_pick < 2 || n_pick >= n_block
        error("cNEDetection: the window has %d blocks of %g sec, too few for resamples of %g%% of them; use a longer window or N_RESAMPLE = 0", ...
              n_block, RESAMPLE_BLOCK, 100*RESAMPLE_FRACTION);
    end
    W = zeros(n_active, 0);  M = {};  NM = {};  thr = zeros(0,2);  nw_mean = zeros(1,0);  nw_sd = zeros(1,0);  run_id = zeros(1,0);
    resamples = struct('blocks', cell(1, N_RESAMPLE), 'nCNE', [], 'lambdaMax', [], 'weights', [], 'members', [], 'icaConverged', []);
    n_converged = 0;    n_runs = 0;
    for r=1:N_RESAMPLE
        pick = sort(randperm(RandStream("twister", "Seed", 1000*r + SEED), n_block, n_pick));   % blocks of the resample, in time order
        rows = reshape((pick - 1)*block_bins + (1:block_bins)', [], 1);    % their bins, joined
        S_r = S(rows,:);
        [mu_r, sd_r] = count_stats(S_r);
        act_r = sd_r > 0;                                                   % (a unit analyzed may be silent in a resample: weight 0)
        idx_r = find(act_r);
        run = single_run(S_r(:,act_r), mu_r(act_r), sd_r(act_r), RandStream("twister", "Seed", 1000*SEED + r), prm, false);
        V_r = zeros(n_active, run.nCNE);    V_r(act_r,:) = run.V;
        mem_r = cellfun(@(m) idx_r(m), run.members, 'UniformOutput', false);             % index into the units analyzed
        neg_r = cellfun(@(m) idx_r(m), run.negativeMembers, 'UniformOutput', false);
        W_r = zeros(n_units, run.nCNE);     W_r(active,:) = V_r;
        resamples(r).blocks = pick;     resamples(r).nCNE = run.nCNE;   resamples(r).lambdaMax = run.lambdaMax;
        resamples(r).weights = W_r;     resamples(r).members = cellfun(@(m) active_idx(m), mem_r, 'UniformOutput', false);
        resamples(r).icaConverged = run.nConverged/max(run.nRuns, 1);
        n_converged = n_converged + run.nConverged;     n_runs = n_runs + run.nRuns;

        keep = cellfun(@numel, mem_r) >= MIN_MEMBERS;                       % components kept for the clustering
        W = [W, V_r(:,keep)];   M = [M, mem_r(keep)];   NM = [NM, neg_r(keep)];                                       %#ok<AGROW>
        thr = [thr; repmat(run.memberThreshold, sum(keep), 1)];                                                       %#ok<AGROW>
        nw_mean = [nw_mean, repmat(run.nullWeight.mean, 1, sum(keep))];     nw_sd = [nw_sd, repmat(run.nullWeight.sd, 1, sum(keep))];   %#ok<AGROW>
        run_id = [run_id, r*ones(1, sum(keep))];                                                                      %#ok<AGROW>
    end

    [groups, centro, similarity] = consensus_clusters(W, run_id, MIN_RESAMPLES, CLUSTER_SIMILARITY);
    stability = cellfun(@(g) numel(unique(run_id(g))), groups)/N_RESAMPLE;
    [~, order] = sortrows([-stability(:), -cellfun(@(c) numel(M{c}), num2cell(centro(:)))]);   % most stable first, then the largest
    centro = centro(order);     stability = stability(order);   similarity = similarity(order);
    V = W(:,centro);
    members = M(centro);    negative_members = NM(centro);
    member_threshold = thr(centro,:);
    null_mean = nw_mean(centro);    null_sd = nw_sd(centro);

    R = cne_activity(S, mu, sd, V);                                         % activity on the whole window, and its null
    if isempty(centro)
        activity_threshold = zeros(1,0);
    else
        activity_threshold = activity_null(S, mu, sd, V, stream, N_SHUFFLE, ACTIVITY_PERCENTILE);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Activation of each cNE on the whole window

n_cne = size(V,2);
activity_events = cell(1,n_cne);
for k=1:n_cne
    activity_events{k} = find(R(:,k) > activity_threshold(k));             % rule (1)
end
if MEMBER_COUNT_RULE
    [members_needed, count_events] = member_count_rule(S, members, stream, COUNT_NULL, N_SHIFT_COUNT, jitter_bins, COUNT_PERCENTILE, ...
                                                       MIN_FIRING);                                                  % rule (2)
    events = cellfun(@(a, c) intersect(a, c), activity_events, count_events, 'UniformOutput', false);
else
    members_needed = nan(1,n_cne);
    count_events = repmat({zeros(0,1)}, 1, n_cne);
    events = activity_events;
end
chance_events = nan(1,n_cne);                                               % (November 2026) chance level (own generator: the results
if N_CHANCE > 0                                                             % above do not change)
    chance_events = chance_level(S, mu, sd, V, members, activity_threshold, MEMBER_COUNT_RULE, COUNT_NULL, N_SHIFT_COUNT, jitter_bins, ...
                                 COUNT_PERCENTILE, MIN_FIRING, N_CHANCE, RandStream("twister", "Seed", 2^20 + SEED));
end

weights = zeros(n_units, n_cne);
weights(active,:) = V;
cNE.nCNE = n_cne;
cNE.weights = weights;
cNE.memberThreshold = reshape(member_threshold, n_cne, 2);
cNE.members = cellfun(@(m) active_idx(m), members, 'UniformOutput', false);
cNE.negativeMembers = cellfun(@(m) active_idx(m), negative_members, 'UniformOutput', false);
if tracking                                                                 % the cNEs as in the template (also the members silent here)
    cNE.weights = TEMPLATE.weights;
    cNE.members = TEMPLATE.members;     cNE.negativeMembers = TEMPLATE.negativeMembers;
end
cNE.stability = stability;
cNE.similarity = similarity;
cNE.activity = R';
cNE.activityThreshold = reshape(activity_threshold, n_cne, 1);
cNE.membersNeeded = members_needed;
cNE.activityEvents = activity_events;
cNE.memberCountEvents = count_events;
cNE.events = events;
cNE.eventTimes = cellfun(@(e) bin_time(e), events, 'UniformOutput', false);
cNE.chanceEvents = chance_events;
cNE.nullWeight = struct('mean', null_mean, 'sd', null_sd);
cNE.icaConverged = n_converged/max(n_runs, 1);
cNE.resamples = resamples;
cNE = orderfields(cNE, {'nCNE', 'nUnits', 'rate', 'active', 'binTime', 'eigenvalues', 'lambdaMax', 'lambdaMP', 'weights', 'memberThreshold', ...
                        'members', 'negativeMembers', 'stability', 'similarity', 'activity', 'activityThreshold', 'membersNeeded', ...
                        'activityEvents', 'memberCountEvents', 'events', 'eventTimes', 'chanceEvents', 'nullWeight', 'icaConverged', ...
                        'resamples', 'parameters'});
if n_converged < n_runs
    warning("cNEDetection: fastICA did not converge in %d of %d runs (ICA_MAX_ITER = %d)", n_runs - n_converged, n_runs, ICA_MAX_ITER);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function run = single_run(S, mu, sd, stream, prm, with_activity)

% steps 2-5 on the binned counts S (bins x units, every unit with spikes; mean mu and standard deviation sd of each unit)
% run.eigenvalues, run.lambdaMax, run.lambdaMP, run.nCNE; run.V (IC weights, units x components); run.members{k}, run.negativeMembers{k}
% (index into the columns of S, sorted by weight); run.memberThreshold, run.nullWeight; with_activity = true: run.R (activity, bins x
% components) and run.activityThreshold; run.nConverged of run.nRuns fastICA runs

[n_bin, n_active] = size(S);
[E, lambda] = sorted_eig(zcorr(S, mu, sd));                                 % principal components of the z-scored counts
run.eigenvalues = lambda;
run.lambdaMP = (1 + sqrt(n_active/n_bin))^2;                                % Marcenko-Pastur upper bound for the eigenvalues of independent units
[nz_bin, nz_unit, nz_count] = find(S);                                      % non-zero entries of the spike count matrix
if prm.N_CNE_RULE == "shift"                                                % bound from circularly shifted data: largest eigenvalue of each shift
    top_eig = zeros(prm.N_SHIFT_EIG, 1);
    for it=1:prm.N_SHIFT_EIG
        shift = randi(stream, n_bin, n_active, 1) - 1;
        top_eig(it) = max(eig(zcorr(shift_counts(nz_bin, nz_unit, nz_count, shift, n_bin, n_active), mu, sd)));
    end
    run.lambdaMax = prctile(top_eig, prm.EIG_PERCENTILE);
else
    run.lambdaMax = run.lambdaMP;
end
n_cne = sum(lambda > run.lambdaMax);
run.nCNE = n_cne;
if n_cne == 0                                                               % no components: nothing else to compute
    run.V = zeros(n_active,0);      run.members = {};   run.negativeMembers = {};
    run.memberThreshold = [NaN NaN];    run.nullWeight = struct('mean', NaN, 'sd', NaN);
    run.R = zeros(n_bin,0);         run.activityThreshold = zeros(1,0);
    run.nConverged = 0;     run.nRuns = 0;
    return
end

[V, converged] = ica_patterns(S, mu, sd, E(:,1:n_cne), lambda(1:n_cne), stream, prm.ICA_CONTRAST, prm.ICA_MAX_ITER, prm.ICA_TOL);   % units x components
if with_activity
    run.R = cne_activity(S, mu, sd, V);                                     % bins x components
end

% null distributions: each unit's binned spike train circularly shifted by a random number of bins
null_weights = zeros(n_active*n_cne, prm.N_SHUFFLE);                        % IC weights of every shuffled matrix
[n_null, rank_p, n_top] = null_rank(n_bin, prm.N_SHUFFLE, prm.ACTIVITY_PERCENTILE);
top_null = zeros(0, n_cne);                                                 % n_top largest null activity values of each component (descending)
n_converged = converged;
for it=1:prm.N_SHUFFLE
    shift = randi(stream, n_bin, n_active, 1) - 1;                          % shift of each unit, in bins (0 to n_bin-1)
    S_shuffle = shift_counts(nz_bin, nz_unit, nz_count, shift, n_bin, n_active);   % (same mean and sd of each unit)

    [E_shuffle, lambda_shuffle] = sorted_eig(zcorr(S_shuffle, mu, sd));    % membership: PCA/ICA with the same number of components
    [V_shuffle, converged] = ica_patterns(S_shuffle, mu, sd, E_shuffle(:,1:n_cne), lambda_shuffle(1:n_cne), stream, prm.ICA_CONTRAST, ...
                                          prm.ICA_MAX_ITER, prm.ICA_TOL);
    null_weights(:,it) = V_shuffle(:);
    n_converged = n_converged + converged;

    if with_activity                                                        % activity: ICs of the real data on the shuffled matrix
        top_null = maxk([top_null; maxk(cne_activity(S_shuffle, mu, sd, V), n_top, 1)], n_top, 1);
    end
end
null_mean = mean(null_weights(:));
null_sd = std(null_weights(:));
run.memberThreshold = null_mean + prm.MEMBER_SD*null_sd*[-1 1];
run.nullWeight = struct('mean', null_mean, 'sd', null_sd);
if with_activity
    run.activityThreshold = percentile_from_top(top_null, n_null, rank_p);
end

run.V = V;
run.members = cell(1,n_cne);    run.negativeMembers = cell(1,n_cne);
for k=1:n_cne
    [w, order] = sort(V(:,k), 'descend');
    run.members{k} = order(w > run.memberThreshold(2));
    run.negativeMembers{k} = order(w < run.memberThreshold(1));
end
run.nConverged = n_converged;
run.nRuns = prm.N_SHUFFLE + 1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [groups, centro, similarity] = consensus_clusters(W, run_id, min_runs, cut)

% clusters of the components (columns of W, unit length; run_id: resample of each one) by the similarity of their weights (cosine, negative
% values set to 0; average linkage, cut at the similarity "cut"); kept when found in at least min_runs resamples. groups{c}: components of
% the cluster c; centro(c): its centrotype (the component with the largest mean similarity to the others); similarity(c): mean similarity
% between its components

groups = {};    centro = zeros(1,0);    similarity = zeros(1,0);
n = size(W,2);
if n == 0
    return
end
if n == 1
    cl = 1;
else
    D = 1 - max(W'*W, 0);   D(1:n+1:end) = 0;   D = (D + D')/2;
    cl = cluster(linkage(squareform(D, 'tovector'), 'average'), 'cutoff', 1 - cut, 'criterion', 'distance');
end
for c=1:max(cl)
    idx = find(cl == c);
    if numel(unique(run_id(idx))) < min_runs
        continue
    end
    Sk = max(W(:,idx)'*W(:,idx), 0);
    [~, best] = max(mean(Sk, 2));
    groups{end+1} = idx;                                                    %#ok<AGROW>
    centro(end+1) = idx(best);                                              %#ok<AGROW>
    if numel(idx) > 1
        similarity(end+1) = mean(Sk(~eye(numel(idx))));                     %#ok<AGROW>
    else
        similarity(end+1) = NaN;                                            %#ok<AGROW>
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function threshold = activity_null(S, mu, sd, V, stream, n_shuffle, percentile)

% threshold of the activation rule (1): percentile of the activity of the ICs V on n_shuffle matrices in which each unit's binned spike
% train is circularly shifted independently

[n_bin, n_active] = size(S);
[nz_bin, nz_unit, nz_count] = find(S);
[n_null, rank_p, n_top] = null_rank(n_bin, n_shuffle, percentile);
top_null = zeros(0, size(V,2));
for it=1:n_shuffle
    shift = randi(stream, n_bin, n_active, 1) - 1;
    R_shuffle = cne_activity(shift_counts(nz_bin, nz_unit, nz_count, shift, n_bin, n_active), mu, sd, V);
    top_null = maxk([top_null; maxk(R_shuffle, n_top, 1)], n_top, 1);
end
threshold = percentile_from_top(top_null, n_null, rank_p);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function chance = chance_level(S, mu, sd, V, members, threshold, member_rule, count_null, n_shift, jitter_bins, percentile, min_firing, ...
                               n_chance, stream)

% (November 2026) mean number of bins in which each cNE passes the activation rules on n_chance surrogates in which every spike of every
% unit analyzed (columns of S) is moved to a random bin of its jitter window (windows of jitter_bins bins from the first bin, as in rule (2)):
% the activations expected by chance, given the rate changes slower than the windows. The thresholds of rule (1) are those of the real data;
% rule (2) is applied to each surrogate as to the real data (its chance level depends only on the spikes of each window, which are kept)

[n_bin, n_active] = size(S);
n_cne = size(V,2);
chance = zeros(1,n_cne);
if n_cne == 0 || n_active == 0
    return
end
[b, u, c] = find(S);
b = repelem(b, c);  u = repelem(u, c);                                      % one entry per spike
start = (ceil(b/jitter_bins) - 1)*jitter_bins;                              % last bin before the window of each spike
len = min(jitter_bins, n_bin - start);                                      % bins of that window (the last one may be shorter)
for it=1:n_chance
    S_jitter = sparse(start + ceil(rand(stream, numel(b), 1).*len), u, 1, n_bin, n_active);
    R = cne_activity(S_jitter, mu, sd, V);
    if member_rule
        [~, count_events] = member_count_rule(S_jitter, members, stream, count_null, n_shift, jitter_bins, percentile, min_firing);
    end
    for k=1:n_cne
        a = find(R(:,k) > threshold(k));
        if member_rule
            a = intersect(a, count_events{k});
        end
        chance(k) = chance(k) + numel(a);
    end
end
chance = chance/n_chance;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [needed, events] = member_count_rule(S, members, stream, count_null, n_shift, jitter_bins, percentile, min_firing)

% activation rule (2): for each cNE (members{k}: index into the columns of S), the number of members with at least one spike in each bin,
% against chance:
%   "jitter": interval jitter. The bins are grouped in windows of jitter_bins bins (from the first bin; the last window may be shorter),
%   and under the null each spike of a member lands in a random bin of its window, independently of the others. Member i then fires in
%   bin b with probability p_i(b) = 1 - (1 - 1/m)^s, where s = its spikes in the window of b and m = the bins of that window, and the
%   number of members firing is Poisson-binomial. events{k} = bins where P(count >= observed) < 1 - percentile/100, with at least
%   min_firing members firing; needed(k) = NaN (the count needed changes from bin to bin)
%   "shift": the same count on n_shift surrogates (the binned train of each member circularly shifted independently). needed(k) = the
%   smallest count above the percentile of the surrogate counts, and at least min_firing; events{k} = bins with at least needed(k) members

n_bin = size(S,1);
n_cne = numel(members);
needed = nan(1,n_cne);
events = cell(1,n_cne);
if count_null == "jitter"
    win = ceil((1:n_bin)'/jitter_bins);                                     % jitter window of each bin
    m_win = accumarray(win, 1);                                             % bins in each window
    to_win = sparse(win, 1:n_bin, 1, win(end), n_bin);                      % sums the bins of each window
    alpha = 1 - percentile/100;
end
for k=1:n_cne
    m = members{k};     n_mem = numel(m);
    if n_mem == 0
        events{k} = zeros(0,1);
        continue
    end
    [b, u] = find(S(:,m));                                                  % bins with spikes of each member
    n_firing = accumarray(b, 1, [n_bin 1]);                                 % members firing in each bin
    if count_null == "jitter"
        s_win = full(to_win*S(:,m));                                        % spikes of each member in each window
        p = 1 - (1 - 1./m_win(win)).^s_win(win,:);                          % P(member fires in the bin) under the null, bins x members
        events{k} = find(poisson_binomial_tail(p, n_firing) < alpha & n_firing >= min_firing);
        continue
    end
    h = zeros(n_mem + 1, 1);                                                % histogram of the surrogate counts (0 to n_mem)
    for s=1:n_shift
        shift = randi(stream, n_bin, n_mem, 1);
        n_surrogate = accumarray(mod(b - 1 + shift(u), n_bin) + 1, 1, [n_bin 1]);
        h = h + accumarray(n_surrogate + 1, 1, [n_mem + 1, 1]);
    end
    chance = find(cumsum(h)/sum(h) >= percentile/100, 1) - 1;               % percentile of the surrogate counts
    needed(k) = max(chance + 1, min_firing);
    events{k} = find(n_firing >= needed(k));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tail = poisson_binomial_tail(p, n_obs)

% P(N >= n_obs(b)) for every row b, where N = number of successes of independent trials with probabilities p(b,1), ..., p(b,n)
% (recursion on the distribution of N, one trial at a time)

[n_row, n] = size(p);
D = zeros(n_row, n + 1);    D(:,1) = 1;                                     % D(b, j+1) = P(N = j)
for i=1:n
    D(:,2:end) = D(:,2:end).*(1 - p(:,i)) + D(:,1:end-1).*p(:,i);
    D(:,1) = D(:,1).*(1 - p(:,i));
end
U = cumsum(D, 2, 'reverse');                                                % U(b, j+1) = P(N >= j)
tail = U(sub2ind([n_row, n + 1], (1:n_row)', min(n_obs, n) + 1));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mu, sd] = count_stats(S)

% mean and standard deviation (as "std") of each column of the sparse counts S (bins x units)

n_bin = size(S,1);
mu = full(sum(S,1))'/n_bin;
sd = sqrt(max(full(sum(S.^2,1))'/n_bin - mu.^2, 0)*n_bin/(n_bin-1));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S_shift = shift_counts(nz_bin, nz_unit, nz_count, shift, n_bin, n_units)

% counts with the train of each unit circularly shifted by shift(unit) bins (non-zero entries nz_bin, nz_unit, nz_count of the counts)

S_shift = sparse(mod(nz_bin - 1 + shift(nz_unit), n_bin) + 1, nz_unit, nz_count, n_bin, n_units);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n_null, rank_p, n_top] = null_rank(n_bin, n_shuffle, percentile)

% the percentile of n_null = n_bin*n_shuffle null values only needs their n_top largest values (rank_p: ascending rank, as in "prctile")

n_null = n_bin*n_shuffle;
rank_p = min(max(n_null*percentile/100 + 0.5, 1), n_null);
n_top = n_null - floor(rank_p) + 1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function threshold = percentile_from_top(top_null, n_null, rank_p)

% percentile from the largest null values (descending, one column per component), with linear interpolation, as "prctile"

lo = floor(rank_p);     hi = ceil(rank_p);
threshold = top_null(n_null-lo+1,:) + (rank_p-lo)*(top_null(n_null-hi+1,:) - top_null(n_null-lo+1,:));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = zcorr(S, mu, sd)

% correlation matrix of the z-scored spike counts (= corr of the columns of S), from the sparse counts S (bins x units)

n_bin = size(S,1);
C = (full(S'*S) - n_bin*(mu*mu')) ./ ((n_bin-1)*(sd*sd'));
C = (C + C')/2;                                                             % exactly symmetric

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E, lambda] = sorted_eig(C)

% eigenvectors (columns of E) and eigenvalues of the symmetric matrix C, in descending order of the eigenvalues

[E, D] = eig(C);
[lambda, order] = sort(diag(D), 'descend');
E = E(:,order);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = zproject(S, mu, sd, P)

% projection of the z-scored spike counts onto the columns of P: Z*P, with Z = (S - mu')./sd' (bins x columns of P)

Y = full(S*(P./sd)) - (mu./sd)'*P;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = cne_activity(S, mu, sd, V)

% activity of each cNE (columns of V) in each bin: R(b,k) = (z(b)'*v_k)^2 - sum_i v_ik^2*z_i(b)^2 (projection matrix with zero diagonal,
% Lopes-dos-Santos et al. 2013), with z(b) the z-scored counts of the bin b

V2 = V.^2./sd.^2;
diagonal = full(S.^2*V2) - 2*full(S*(V2.*mu)) + sum(V2.*mu.^2, 1);         % sum_i v_ik^2*z_i(b)^2, with z = (s - mu)/sd
R = zproject(S, mu, sd, V).^2 - diagonal;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V, converged] = ica_patterns(S, mu, sd, P, lambda, stream, contrast, max_iter, tol)

% IC weights (columns of V, one per component): fastICA on the z-scored counts projected onto the principal components P (eigenvalues
% lambda), as unmixing vectors in the space of the units; each one scaled to unit length, with its largest absolute weight positive

X = zproject(S, mu, sd, P) ./ sqrt(lambda');                                % whitened data (unit variance, uncorrelated), bins x components
[U, converged] = fast_ica(X, stream, contrast, max_iter, tol);              % rows of U: unmixing vectors in the whitened space
V = P * ((1./sqrt(lambda)) .* U');                                          % unmixing vectors in the space of the units: ic = Z*V

V = V ./ vecnorm(V);                                                        % unit length
[~, imax] = max(abs(V), [], 1);
V = V .* sign(V(sub2ind(size(V), imax, 1:size(V,2))));                     % largest absolute weight positive

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U, converged] = fast_ica(X, stream, contrast, max_iter, tol)

% symmetric fastICA (Hyvarinen 1999, IEEE Trans Neural Netw 10:626) on whitened data X (samples x components)
% U: orthogonal unmixing matrix (one component per row; the independent components are X*U')

[n_sample, n_comp] = size(X);
U = sym_decorrelation(randn(stream, n_comp, n_comp));                       % random orthogonal initialization
converged = false;
for it=1:max_iter
    Y = X*U';                                                               % current estimate of the independent components
    switch contrast
        case "logcosh"
            G = tanh(Y);                dG = 1 - G.^2;
        case "pow3"
            G = Y.^3;                   dG = 3*Y.^2;
        case "gauss"
            e = exp(-Y.^2/2);           G = Y.*e;       dG = (1 - Y.^2).*e;
        otherwise
            error("cNEDetection: unknown ICA_CONTRAST '%s' (use ""logcosh"", ""pow3"" or ""gauss"")", contrast);
    end
    U_new = sym_decorrelation((G'*X)/n_sample - mean(dG,1)'.*U);            % fixed-point update of all the components at once
    change = max(abs(abs(sum(U_new.*U, 2)) - 1));                           % 1 - |cos(angle)| between consecutive estimates
    U = U_new;
    if change < tol
        converged = true;
        break
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = sym_decorrelation(U)

% symmetric decorrelation: U = (U*U')^(-1/2)*U (the rows become orthonormal)

[E, D] = eig((U*U' + (U*U')')/2);
U = (E ./ sqrt(diag(D))') * E' * U;

end
