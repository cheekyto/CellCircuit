function cluster_neuronal_events = cluster_neuronal_events(NetworkBurstPlayer, NetworkBurstFeature, spikeTimestamp, range, option_plot)

%% cluster_neuronal_events.m
%
% Unsupervised clustering of neuronal events (families of network burst events) using mixed metrics:
%   - Scalar metrics (e.g., duration, amplitude, firing rate, ...): lower weight
%   - One population metric: the participating units of each event (variable length): higher weight
%
% (revised November 2026; user's choices: rules 1, 3, 4 and 5 below; the scalar metrics are kept, with their lower weight)
% The approach:
%   1. Participating units of each event (rule 1): a unit participates in an event when it fires significantly more spikes within the
%      event (start to last) than expected from its own firing rate outside the events of the condition (Poisson test, p <
%      PARTICIPATION_ALPHA, and at least PARTICIPATION_MIN_SPIKES spikes). Before, a unit participated when it fired at least one spike
%      within the window of the event (with the extra time around it of "eventSpikes"): in long events, most units fired by chance, so
%      events driven by different units looked alike, and all the large events had all the units
%   2. Normalize the scalar metrics (z-score) and compute their pairwise Euclidean distance, rescaled to [0,1] (as before)
%   3. Pairwise Jaccard distance between the participating units of the events (as before, with the units of step 1)
%   4. Composite distance: weighted average of 2 and 3 (w_scalar, w_population; as before)
%   5. Families (rule 3): average-linkage tree of the composite distance, cut where the events of a family are, on average, at least
%      MIN_SIMILARITY similar (composite distance <= 1 - MIN_SIMILARITY). Before, Ward linkage on a classical MDS projection, with 2 to
%      MAX_K families chosen with the silhouette score, a bootstrap of the events and a permutation test ("bootstrap_optimal_k",
%      "permutation_test_clustering"): it always looked for 2 or more families and created families of 1-2 events (overfitting)
%   6. Minimum size (rule 4): a family needs at least MIN_FAMILY_EVENTS events; the events of smaller groups belong to no family
%      (label 0). If no family is found, all the events have label 0 (in the per-family results, all the primary events are then a
%      single group: see "primaryFamilies")
%   7. Stability (rule 5): the participating units are resampled with replacement N_RESAMPLE times and the families are found again;
%      the stability of a family is the mean probability that two of its events are again in the same family
%   8. Visualize (option_plot): dendrogram with the cut, scalar metrics sorted by family, raster with the events coloured by family

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   cluster_neuronal_events.n_cluster                   number of families (revised November 2026: 0 when no family was found)
%   cluster_neuronal_events.label(i,1)                  family of the event i (revised November 2026: 0 = the event belongs to no family)
%                                (i,2)                  =1, event i is a primary event; =2, event i is a secondary event
%   cluster_neuronal_events.clusteringQuality           (revised November 2026) stability of the least stable family (0 to 1, see step 7);
%                                                       =0 when no family was found. Before: =0, the clustering was non-significant (one
%                                                       single family); =1, it seems to work; =2, the clustering is reliable
%   cluster_neuronal_events.stability(f)                (new, November 2026) stability of the family f (0 to 1)
%   cluster_neuronal_events.participation(i,j)          (new, November 2026) =true, the unit units(j) participates in the event i (step 1)
%   cluster_neuronal_events.units(j)                    (new, November 2026) units tested (participating units of the primary or secondary events)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_SCALAR = 7;
MIN_SECONDARY_EVENTS = 5;       % minimum number of secondary events to include them

colors_event = {'r', 'b', 'c', 'g', 'm'};

% --- Weighting ------------------------------------------------------------
% w_scalar + w_population = 1
% Increase w_population if the neuron identity is a stronger grouping signal.
w_scalar     = 0.35;            % (default: 0.35)
w_population = 0.65;            % (default: 0.65)

% --- Participating units of each event (new, November 2026) --------------
PARTICIPATION_ALPHA = 0.01;     % a unit participates in an event when the probability of firing as many spikes by chance is below this value (default: 0.01)
PARTICIPATION_MIN_SPIKES = 2;   % and it fires at least this number of spikes within the event (default: 2)

% --- Families (revised November 2026) -------------------------------------
LINKAGE_METHOD = 'average';     % linkage of the tree of the composite distance (before: 'ward', on an MDS projection)
MIN_SIMILARITY = 0.50;          % the events of a family are, on average, at least this similar (1 - composite distance) (default: 0.50)
MIN_FAMILY_EVENTS = 3;          % minimum number of events of a family; the events of smaller groups belong to no family (default: 3)

% --- Stability of the families (new, November 2026) -----------------------
N_RESAMPLE = 200;               % number of resamples of the participating units (default: 200)
RESAMPLE_SEED = 0;              % seed of the resamples (the same seed gives the same result)
apply_parameters('cluster_neuronal_events');                                % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

K_MSEC_SEC = 1/1000;                        % factor to rescale from msec to sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_events1 = NetworkBurstPlayer.nBurst;
n_events2 = NetworkBurstPlayer.Secondary.nBurst;
if isnan(n_events2) || n_events2 <= MIN_SECONDARY_EVENTS
    n_events2 = 0;
end
n_events = n_events1 + n_events2;

scalar_metrics = nan(n_events, N_SCALAR);

% --- Scalar metrics (n_events x n_scalar_metrics) -------------------------
for k=n_events1:-1:1
    scalar_metrics(k,1) = NetworkBurstFeature{1}.indivMetric(k).duration;
    scalar_metrics(k,2) = NetworkBurstFeature{1}.indivMetric(k).nSpikes;
    scalar_metrics(k,3) = NetworkBurstFeature{1}.indivMetric(k).nNeuronsAct;
    scalar_metrics(k,4) = NetworkBurstFeature{1}.indivMetric(k).isi;
    if k>1
        scalar_metrics(k,5) = NetworkBurstPlayer.IndivBurst(k).peak - NetworkBurstPlayer.IndivBurst(k-1).peak;
    else
        scalar_metrics(k,5) = max(median(scalar_metrics(:,5),'omitnan'), NetworkBurstPlayer.IndivBurst(k).peak - range(1));
    end
%    if k<n_events
%        scalar_metrics(k,6) = NetworkBurstPlayer.IndivBurst(k+1).peak - NetworkBurstPlayer.IndivBurst(k).peak;
%    end
    scalar_metrics(k,6) = NetworkBurstFeature{1}.indivMetric(k).tail;
    scalar_metrics(k,7) = NetworkBurstFeature{1}.indivMetric(k).amplitudePop;
%    scalar_metrics(k,8) = NetworkBurstFeature{1}.indivMetric(k).level;
end

if n_events2>MIN_SECONDARY_EVENTS
    for k=n_events2:-1:1
        idx = k+n_events1;
        scalar_metrics(idx,1) = NetworkBurstFeature{2}.indivMetric(k).duration;
        scalar_metrics(idx,2) = NetworkBurstFeature{2}.indivMetric(k).nSpikes;
        scalar_metrics(idx,3) = NetworkBurstFeature{2}.indivMetric(k).nNeuronsAct;
        scalar_metrics(idx,4) = NetworkBurstFeature{2}.indivMetric(k).isi;
        if k>1
            scalar_metrics(idx,5) = NetworkBurstPlayer.Secondary.IndivBurst(k).peak - NetworkBurstPlayer.Secondary.IndivBurst(k-1).peak;
        else
            scalar_metrics(idx,5) = max(scalar_metrics(:,5));
        end
        scalar_metrics(idx,6) = NetworkBurstFeature{2}.indivMetric(k).tail;
        scalar_metrics(idx,7) = NetworkBurstFeature{2}.indivMetric(k).amplitudePop;
%       scalar_metrics(idx,8) = 2;
    end
end
% (revised November 2026) "scalar_metrics2" (the scalar metrics plus the level of the primary events) was built but never used: removed

%scalar_metric_names = {'Duration','Nspikes','NActiveNeurons','isi','distance_pre','distance_post'};
scalar_metric_names = {'Duration','Nspikes','NActiveNeurons','isi', 'Distance_to_pre','tail', 'amplitude_Pop'};

% --- Population metric: participating units of each event -----------------
% (revised November 2026; rule 1) a unit participates in an event when it fires significantly more than its own rate outside the events
% (see the local function "event_participation"); before, when it fired at least one spike within the window of the event
total_neurons = size(NetworkBurstPlayer.timeStamp,2);
units = NetworkBurstPlayer.ChannelsIntegrated(:,1);                         % participating units of the primary or secondary events (as before)
starts = [NetworkBurstPlayer.IndivBurst.start]';                            % events to be clustered, in sec
lasts = [NetworkBurstPlayer.IndivBurst.last]';
if n_events2 > 0
    starts = [starts; [NetworkBurstPlayer.Secondary.IndivBurst.start]'];
    lasts = [lasts; [NetworkBurstPlayer.Secondary.IndivBurst.last]'];
end
starts_all = [NetworkBurstPlayer.IndivBurst.start]';                        % all the events of the condition (primary and secondary, also when
lasts_all = [NetworkBurstPlayer.IndivBurst.last]';                          % the secondary events are not clustered), to measure the rate of the
if ~isnan(NetworkBurstPlayer.Secondary.nBurst) && NetworkBurstPlayer.Secondary.nBurst > 0     % units outside them
    starts_all = [starts_all; [NetworkBurstPlayer.Secondary.IndivBurst.start]'];
    lasts_all = [lasts_all; [NetworkBurstPlayer.Secondary.IndivBurst.last]'];
end
pop_binary = event_participation(spikeTimestamp, units, starts, lasts, starts_all, lasts_all, range, PARTICIPATION_ALPHA, PARTICIPATION_MIN_SPIKES);
neuron_population = cell(n_events, 1);                                      % (as before) index of the participating units of each event
for k = 1:n_events
    neuron_population{k} = units(pop_binary(k,:));
end

%% =========================================================================
%  STEP 1 – Normalize scalar metrics
% =========================================================================
% Z-score each column; replace NaN with 0 (handles constant columns)
% scalar_norm = zscore(scalar_metrics);
scalar_norm = normalize(scalar_metrics);
scalar_norm(isnan(scalar_norm)) = 0;

% Pairwise Euclidean distance on normalized scalars (condensed form)
D_scalar = pdist(scalar_norm, 'euclidean');

% Rescale to [0,1] for comparability with Jaccard
D_scalar_norm = D_scalar ./ max(D_scalar + eps);

%% =========================================================================
%  STEP 2 – Families (revised November 2026; rules 3 and 4)
% =========================================================================
% Jaccard distance between the participating units (1 - |A∩B| / |A∪B|; two events without participating units: 1) and composite distance
% (inside "find_families"), average-linkage tree, cut at 1 - MIN_SIMILARITY, and groups with at least MIN_FAMILY_EVENTS events
[label1, Z] = find_families(pop_binary, D_scalar_norm, w_scalar, w_population, LINKAGE_METHOD, MIN_SIMILARITY, MIN_FAMILY_EVENTS);
n_families = max([label1; 0]);

%% =========================================================================
%  STEP 3 – Stability of the families (new, November 2026; rule 5)
% =========================================================================
% The participating units are resampled with replacement (the scalar metrics do not change) and the families are found again. The
% stability of a family is the mean, across all the pairs of its events, of the proportion of resamples in which both events are in the
% same family. The random numbers come from their own stream: the global random number generator is not used
stability = nan(n_families, 1);
if n_families > 0 && N_RESAMPLE > 0 && ~isempty(units)
    stream = RandStream('twister', 'Seed', RESAMPLE_SEED);
    n_units = numel(units);
    together = zeros(n_events);
    for b = 1:N_RESAMPLE
        pick = randi(stream, n_units, n_units, 1);                          % units of the resample (with replacement)
        label_b = find_families(pop_binary(:,pick), D_scalar_norm, w_scalar, w_population, LINKAGE_METHOD, MIN_SIMILARITY, MIN_FAMILY_EVENTS);
        together = together + (label_b == label_b' & label_b > 0);
    end
    together = together / N_RESAMPLE;
    for f = 1:n_families
        member = find(label1 == f);
        pairs = together(member, member);
        stability(f) = mean(pairs(~eye(numel(member))));
    end
end

if n_families > 0
    quality_clustering = min(stability);                                    % stability of the least stable family
else
    quality_clustering = 0;                                                 % no family
end

cluster_neuronal_events.n_cluster = n_families;
label = nan(n_events,2);
label(:,1) = label1;
label(1:n_events1, 2) = 1;
if n_events2>0
    label((n_events1+1):n_events, 2) = 2;
end
cluster_neuronal_events.label = label;
cluster_neuronal_events.clusteringQuality = quality_clustering;
cluster_neuronal_events.stability = stability;
cluster_neuronal_events.participation = pop_binary;
cluster_neuronal_events.units = units;

%% =========================================================================
%  STEP 4 – Visualisation
% =========================================================================

if option_plot

cmap = lines(max(n_families,1));   % color palette

% ----- Figure 1: Dendrogram, with the cut that defines the families ----------
if n_events > 1
figure('Name','Dendrogram','Color','w','Position',[100 100 700 420]);
dendrogram(Z, 0, 'ColorThreshold', 1 - MIN_SIMILARITY);
hold on; yline(1 - MIN_SIMILARITY, 'r--', sprintf('cut (similarity %.2f)', MIN_SIMILARITY));
title(sprintf('Families of events  (%s linkage; %d families, %d events in none)', LINKAGE_METHOD, n_families, sum(label1==0)));
xlabel('Event index'); ylabel('Composite distance');
set(gca,'FontSize',11); box off;
end

% ----- Figure 2: Scalar metric heatmap sorted by family -----------------
figure('Name','Metric heatmap','Color','w','Position',[100 100 820 480]);
[~, sort_idx] = sort(label1);
imagesc(scalar_norm(sort_idx, :));
colormap(redblue(256));
colorbar;
xlabel('Metric'); ylabel('Event (sorted by family; 0 = no family)');
title('Normalized scalar metrics  –  sorted by family');
xticks(1:numel(scalar_metric_names));
xticklabels(scalar_metric_names);
xtickangle(30);
set(gca,'FontSize',10);

% Overlay family boundaries
boundaries = find(diff(label1(sort_idx))) + 0.5;
for b = boundaries'
    yline(b, 'w-', 'LineWidth', 2);
end

% ----- Figure 3: Raster plot of participating neurons and Population Firing -----------------
figure('Name','Raster Plot and Population Firing','Color','w');

if ~isnan(NetworkBurstPlayer.nBurst) && NetworkBurstPlayer.nBurst>0 && numel(NetworkBurstPlayer.Channels(:,1))>0    % (revised October 2026: "&&"; with "||" it was always true)

    subplot(5,1,[1,2,3]);
    for i=1:numel(NetworkBurstPlayer.Channels(:,1))
        data = spikeTimestamp(NetworkBurstPlayer.Channels(i,1));
        n_unit = size(data.timeStamp,1);
        plot(data.timeStamp(:,1)*K_MSEC_SEC, ones(n_unit,1)+i,'|',"MarkerEdgeColor",[0 0 0]); hold on
    end
    ylabel('Neurons');  title('Selected units - Participating in network burst');  xlim([range(1) range(2)]); ylim([1 i+2])

    subplot(5,1,[4,5]);
    bar(NetworkBurstPlayer.time(:),NetworkBurstPlayer.spikesBin(:,1)); hold on;
    xlabel('Time (sec)');   ylabel('Number of spikes');  xlim([range(1) range(2)])

    for idx = 1:n_events
        if idx>n_events1
            peak_time = NetworkBurstPlayer.Secondary.IndivBurst(idx-n_events1).peak;
            peak_Amp = NetworkBurstPlayer.Secondary.IndivBurst(idx-n_events1).peakAmp;
        else
            peak_time = NetworkBurstPlayer.IndivBurst(idx).peak;
            peak_Amp = NetworkBurstPlayer.IndivBurst(idx).peakAmp;
        end
        if label1(idx) > 0                                                  % (revised November 2026) the colors repeat after 5 families; grey: no family
            color_k = colors_event{mod(label1(idx)-1, numel(colors_event))+1};
        else
            color_k = [0.6 0.6 0.6];
        end
        plot(peak_time, peak_Amp*1.1,'o', 'MarkerFaceColor',color_k); hold on
    end

end


%% =========================================================================
%  STEP 5 – Print family summary to console
% =========================================================================
fprintf('%-8s  %-10s  %s\n', 'Event', 'Family', 'Participating units');
fprintf('%s\n', repmat('-', 1, 55));
for i = 1:n_events
    fprintf('%-8d  %-10d  %d\n', i, label1(i), sum(pop_binary(i,:)));
end

fprintf('\n--- Family means (scalar metrics; family 0 = events in no family) ---\n');
fprintf('%-12s', '', scalar_metric_names{:});
fprintf('\n');
for k = unique(label1)'
    idx = label1 == k;
    mu  = mean(scalar_metrics(idx,:), 1, 'omitnan');
    fprintf('Family %-4d  ', k);
    fprintf('%-12.3f', mu);
    fprintf('\n');
end
fprintf('\nStability of the families: %s\n', mat2str(round(stability',2)));

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function participation = event_participation(spikeTimestamp, units, starts, lasts, starts_all, lasts_all, range, alpha, min_spikes)

% (new, November 2026; rule 1) participating units of each event: the unit units(j) participates in the event k when it fires at least
% "min_spikes" spikes from starts(k) to lasts(k) and the probability of firing as many spikes by chance is below "alpha" (Poisson
% distribution, with the rate of the unit outside all the events of the condition, "starts_all" to "lasts_all", within "range")
%   spikeTimestamp(ch).timeStamp(:,1):  spike times of the cluster ch, in msec
%   starts, lasts:          events to be clustered, in sec
%   participation(k,j):     =true, the unit units(j) participates in the event k

K_MSEC_SEC = 1/1000;
n_events = numel(starts);
participation = false(n_events, numel(units));
if n_events == 0 || isempty(units)
    return
end

% time covered by the events of the condition (overlapping events counted once), within the condition window
win = sortrows([max(starts_all(:), range(1)), min(lasts_all(:), range(2))]);
win = win(win(:,2) > win(:,1), :);
merged = zeros(0,2);
for k = 1:size(win,1)
    if ~isempty(merged) && win(k,1) <= merged(end,2)
        merged(end,2) = max(merged(end,2), win(k,2));
    else
        merged(end+1,:) = win(k,:); %#ok<AGROW>
    end
end
time_outside = (range(2) - range(1)) - sum(merged(:,2) - merged(:,1));     % time outside the events, in sec

duration = lasts(:) - starts(:);                                            % duration of each event, in sec
for j = 1:numel(units)
    spikes = spikeTimestamp(units(j)).timeStamp(:,1);                       % in msec
    spikes = spikes(spikes >= range(1)/K_MSEC_SEC & spikes < range(2)/K_MSEC_SEC);
    n_inside = sum(count_below(spikes, merged(:,2)/K_MSEC_SEC, 1) - count_below(spikes, merged(:,1)/K_MSEC_SEC, 0));   % spikes within the events
    rate_outside = (numel(spikes) - n_inside) / max(time_outside, eps);    % spikes per sec outside the events
    n_spikes = count_below(spikes, lasts(:)/K_MSEC_SEC, 1) - count_below(spikes, starts(:)/K_MSEC_SEC, 0);   % spikes within each event (start to last;
                                                                            % each event on its own: the events are not sorted, and can overlap)
    p = poisscdf(n_spikes - 1, rate_outside*duration, 'upper');             % probability of at least n_spikes spikes by chance
    participation(:,j) = n_spikes >= min_spikes & p < alpha;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [label, Z] = find_families(pop_binary, D_scalar_norm, w_scalar, w_population, linkage_method, min_similarity, min_events)

% (new, November 2026; rules 3 and 4) families of events: average-linkage tree of the composite distance, cut at 1 - min_similarity;
% the groups with at least min_events events are the families, numbered in the order of their first event; the other events get 0

n_events = size(pop_binary, 1);
label = zeros(n_events, 1);
Z = [];
if n_events < 2
    return
end
D_population = pdist(double(pop_binary), 'jaccard');                        % Jaccard distance: 1 - |A∩B| / |A∪B|
D_population(isnan(D_population)) = 1;                                      % two empty sets → max distance
D_composite = w_scalar * D_scalar_norm + w_population * D_population;
Z = linkage(D_composite, linkage_method);
group = cluster(Z, 'cutoff', 1 - min_similarity, 'criterion', 'distance');
[~, first] = unique(group, 'first');                                        % first event of each group
for g = group(sort(first))'                                                 % groups, in the order of their first event
    member = group == g;
    if sum(member) >= min_events
        label(member) = max(label) + 1;
    end
end
end
