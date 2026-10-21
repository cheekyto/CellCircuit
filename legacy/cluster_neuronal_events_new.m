%% cluster_neuronal_events.m
%
% Unsupervised clustering of neuronal events using mixed metrics:
%   - Scalar metrics (e.g., duration, amplitude, firing rate, ...)
%   - One population metric: array of participating neuron labels (variable length)
%
% The approach:
%   1. Normalize scalar metrics (z-score)
%   2. Compute pairwise Jaccard distance for the neuron population metric
%   3. Combine into a single composite distance matrix (weighted average)
%   4. Project composite distances to Euclidean space via classical MDS
%      (required for Ward linkage, which assumes Euclidean geometry)
%   5. Run agglomerative hierarchical clustering + evaluate k=1..MAX_K via
%      silhouette score to select the optimal number of clusters
%   6. Visualize: dendrogram, silhouette plot, and PCoA (metric MDS) scatter
%
% NOTE ON LINKAGE FALLBACK:
%   Even after MDS projection, Ward linkage occasionally produces a
%   non-monotonic tree (MATLAB warning) when the Euclidean embedding does
%   not fully capture the composite distance structure. The safe_linkage()
%   function below detects this at runtime and automatically falls back to
%   average linkage, which is valid for any distance matrix.
%
% -------------------------------------------------------------------------
% INPUTS  –– edit the "USER INPUT" section below
% -------------------------------------------------------------------------

clear; close all; clc;

%% =========================================================================
%  USER INPUT  ––  Replace this section with your actual data
% =========================================================================

n_events = 20;

rng(42);
scalar_metrics = [
    rand(n_events, 1) * 500,
    rand(n_events, 1) * 2,
    rand(n_events, 1) * 50,
    rand(n_events, 1),
    rand(n_events, 1) * 10,
    randn(n_events, 1),
    rand(n_events, 1) * 100,
    rand(n_events, 1) * 5
];

scalar_metric_names = {'Duration','Amplitude','FiringRate','Synchrony', ...
                       'SpatialSpread','LFP_Power','IEI','OscFreq'};

total_neurons = 30;
neuron_population = cell(n_events, 1);
for i = 1:n_events
    k = randi([3, 15]);
    neuron_population{i} = sort(randperm(total_neurons, k));
end

w_scalar     = 0.6;
w_population = 0.4;

MAX_K          = 3;
LINKAGE_METHOD = 'ward';

%% =========================================================================
%  STEP 1 – Normalize scalar metrics
% =========================================================================
scalar_norm = zscore(scalar_metrics);
scalar_norm(isnan(scalar_norm)) = 0;

D_scalar      = pdist(scalar_norm, 'euclidean');
D_scalar_norm = D_scalar ./ max(D_scalar + eps);

%% =========================================================================
%  STEP 2 – Jaccard distance for neuron population metric
% =========================================================================
pop_binary = false(n_events, total_neurons);
for i = 1:n_events
    pop_binary(i, neuron_population{i}) = true;
end

D_population = pdist(pop_binary, 'jaccard');
D_population(isnan(D_population)) = 1;

%% =========================================================================
%  STEP 3 – Composite distance matrix
% =========================================================================
D_composite = w_scalar * D_scalar_norm + w_population * D_population;
D_mat       = squareform(D_composite);

%% =========================================================================
%  STEP 4 – Project to Euclidean space via classical MDS
% =========================================================================
[Y_mds, eigvals] = cmdscale(D_mat);
Y_mds = Y_mds(:, eigvals > 1e-10);

n_neg = sum(eigvals < -1e-10);
if n_neg > 0
    fprintf(['Note: %d negative MDS eigenvalue(s) discarded ' ...
             '(expected for composite/non-Euclidean distances).\n\n'], n_neg);
end

var_explained = eigvals(eigvals > 1e-10) / sum(eigvals(eigvals > 0)) * 100;

%% =========================================================================
%  STEP 5 – Hierarchical clustering + silhouette-based model selection
% =========================================================================
[Z, method_used] = safe_linkage(Y_mds, LINKAGE_METHOD);

fprintf('\n%-10s  %-15s\n', 'K clusters', 'Mean silhouette');
fprintf('%s\n', repmat('-', 1, 28));

sil_scores    = zeros(1, MAX_K);
sil_scores(1) = NaN;

for k = 2:MAX_K
    labels_k      = cluster(Z, 'maxclust', k);
    sil_scores(k) = mean(silhouette(Y_mds, labels_k));
    fprintf('%-10d  %-15.4f\n', k, sil_scores(k));
end

[~, best_k] = max(sil_scores);
fprintf('\n>> Best K = %d  (mean silhouette = %.4f)\n\n', best_k, sil_scores(best_k));

event_labels = cluster(Z, 'maxclust', best_k);

%% =========================================================================
%  STEP 6 – Visualisation
% =========================================================================
cmap = lines(MAX_K);

% ----- Figure 1: Dendrogram ----------------------------------------------
figure('Name','Dendrogram','Color','w','Position',[100 100 700 420]);
dendrogram(Z, 0, 'ColorThreshold', 'default');
title(sprintf('Dendrogram  (%s linkage)', method_used));
xlabel('Event index'); ylabel('Distance');
set(gca,'FontSize',11); box off;

% ----- Figure 2: Silhouette scores per K ---------------------------------
figure('Name','Silhouette scores','Color','w','Position',[820 100 480 340]);
bar(2:MAX_K, sil_scores(2:end), 0.5, 'FaceColor',[0.35 0.60 0.85]);
hold on;
plot(best_k, sil_scores(best_k), 'r*', 'MarkerSize', 14, 'LineWidth', 2);
xlabel('Number of clusters K'); ylabel('Mean silhouette score');
title(sprintf('Model selection: silhouette score vs K  (%s linkage)', method_used));
xticks(2:MAX_K); ylim([0 1]); grid on; box off;
set(gca,'FontSize',11);

% ----- Figure 3: Per-event silhouette plot (best K) ----------------------
s_best = silhouette(Y_mds, event_labels);
figure('Name',sprintf('Silhouette  K=%d', best_k),'Color','w', ...
       'Position',[100 560 500 420]);
[~, srt] = sort(event_labels);
barh(s_best(srt), 'FaceColor', 'flat', 'CData', cmap(event_labels(srt),:));
xlabel('Silhouette value'); ylabel('Events (sorted by cluster)');
title(sprintf('Silhouette plot  (K = %d,  %s linkage)', best_k, method_used));
set(gca,'FontSize',11); box off;

% ----- Figure 4: PCoA scatter (first 2 MDS dimensions) ------------------
figure('Name','PCoA scatter','Color','w','Position',[620 560 580 460]);
hold on;
for k = 1:best_k
    idx = event_labels == k;
    scatter(Y_mds(idx,1), Y_mds(idx,2), 90, cmap(k,:), 'filled', ...
            'DisplayName', sprintf('Cluster %d (n=%d)', k, sum(idx)));
end
for i = 1:n_events
    text(Y_mds(i,1), Y_mds(i,2), sprintf('  %d',i), 'FontSize', 8, ...
         'Color', cmap(event_labels(i),:));
end
legend('Location','best'); grid on; box off;
xlabel(sprintf('PCoA 1  (%.1f%%)', var_explained(1)));
ylabel(sprintf('PCoA 2  (%.1f%%)', var_explained(2)));
title(sprintf('PCoA projection  –  K = %d clusters  (%s linkage)', ...
              best_k, method_used));
set(gca,'FontSize',11);

% ----- Figure 5: Scalar metric heatmap sorted by cluster -----------------
figure('Name','Metric heatmap','Color','w','Position',[100 100 820 480]);
[~, sort_idx] = sort(event_labels);
imagesc(scalar_norm(sort_idx, :));
colormap(redblue(256));
colorbar;
xlabel('Metric'); ylabel('Event (sorted by cluster)');
title(sprintf('Normalized scalar metrics  –  sorted by cluster  (%s linkage)', ...
              method_used));
xticks(1:numel(scalar_metric_names));
xticklabels(scalar_metric_names);
xtickangle(30);
set(gca,'FontSize',10);
boundaries = find(diff(event_labels(sort_idx))) + 0.5;
for b = boundaries'
    yline(b, 'w-', 'LineWidth', 2);
end

%% =========================================================================
%  STEP 7 – Print cluster summary to console
% =========================================================================
fprintf('%-8s  %-10s  %s\n', 'Event', 'Cluster', 'Neuron IDs');
fprintf('%s\n', repmat('-', 1, 55));
for i = 1:n_events
    fprintf('%-8d  %-10d  [%s]\n', i, event_labels(i), ...
            num2str(neuron_population{i}));
end

fprintf('\n--- Cluster means (scalar metrics) ---\n');
fprintf('%-12s', '', scalar_metric_names{:});
fprintf('\n');
for k = 1:best_k
    idx = event_labels == k;
    mu  = mean(scalar_metrics(idx,:), 1);
    fprintf('Cluster %-3d  ', k);
    fprintf('%-12.3f', mu);
    fprintf('\n');
end

%% =========================================================================
%  LOCAL FUNCTION – safe linkage with automatic fallback
% =========================================================================
function [Z, method_used] = safe_linkage(Y, preferred_method)
% Attempts linkage with preferred_method (e.g. 'ward').
% If MATLAB raises a non-monotonic-tree warning (which means the Euclidean
% MDS embedding did not fully resolve the composite distance geometry),
% falls back automatically to 'average' linkage, which is valid for any
% distance matrix and does not require a monotonic tree.
%
% Returns Z (linkage matrix) and method_used (string) so callers can
% report which method was actually applied.

    FALLBACK = 'average';

    % Temporarily promote the non-monotonic warning to an error so we can
    % catch it cleanly without it printing to the console
    warnID = 'stats:linkage:NonMonotonicTree';
    prev   = warning('error', warnID);

    try
        Z           = linkage(Y, preferred_method);
        method_used = preferred_method;
    catch
        warning(prev);   % restore warning state before the fallback call
        fprintf(['Warning: Ward linkage produced a non-monotonic tree ' ...
                 '(MDS embedding insufficient).\n' ...
                 'Falling back to ''%s'' linkage automatically.\n\n'], FALLBACK);
        Z           = linkage(Y, FALLBACK);
        method_used = FALLBACK;
    end

    warning(prev);   % always restore warning state on clean exit too
end

%% =========================================================================
%  LOCAL FUNCTION – diverging red-white-blue colormap
% =========================================================================
function c = redblue(n)
    half = floor(n/2);
    r1 = [linspace(0.7,1,half)', linspace(0.1,1,half)', linspace(0.1,1,half)'];
    r2 = [linspace(1,0.1,n-half)', linspace(1,0.1,n-half)', linspace(1,0.7,n-half)'];
    c  = [r1; r2];
end
