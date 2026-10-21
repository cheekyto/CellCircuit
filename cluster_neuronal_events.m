function cluster_neuronal_events = cluster_neuronal_events(NetworkBurstPlayer, NetworkBurstFeature, spikeTimestamp, range, option_plot)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   cluster_neuronal_events.n_cluster                   number of clusters   
%   cluster_neuronal_events.label(i,1)                  id of the cluster for unit i
%                                (i,2)                  =1, unit i is a primary prak; =2, unit i is a secondary peak
%   cluster_neuronal_events.clusteringQuality           =0, the clustering was non-significant (one single family)
%                                                       =1, the clustering seems to work, for the boostrapping and the permutations did not provide clear answer
%                                                       =2, the clustering is reliable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_SCALAR = 7;
MIN_SECONDARY_EVENTS = 5;       % minimum number of secondary events to include them

MIN_SIZE_CLUSTER = 0;           % minimum size of the clusters to run the boostrapping (default: 3)

colors_event = {'r', 'b', 'c', 'g', 'm'};

% --- Weighting ------------------------------------------------------------
% w_scalar + w_population = 1
% Increase w_population if the neuron identity is a stronger grouping signal.
w_scalar     = 0.35;            % (default: 0.35)
w_population = 0.65;            % (default: 0.65)

% --- Clustering settings --------------------------------------------------
MAX_K          = 4;       % maximum number of clusters to evaluate (default: 4)
LINKAGE_METHOD = 'ward';  % 'ward' | 'average' | 'complete' | 'single'

% --- Permutation settings --------------------------------------------------
ALPHA          = 0.05;    % significance threshold (default: 0.05)

% --- Integration settings --------------------------------------------------
MIN_DIFF = 0.10;            % maximum difference in "vote_frac" and "silhouette score" to consider to k-sizes non-different (default: 0.10)
SIL_THRESHOLD  = 0.40;      % AND have a mean silhouette >= 0.25 (0.25 is a modest but non-trivial structure threshold) (the higher the value, the more conservative)
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

scalar_metrics2 = nan(n_events, N_SCALAR+1);
scalar_metrics2(:,1:N_SCALAR) = scalar_metrics;
for k=1:n_events1
    scalar_metrics2(k,8) = NetworkBurstFeature{1}.indivMetric(k).level;
end
% (revised October 2026) the line "scalar_metrics2 = scalar_metrics2" (without ";") printed the whole matrix on every run

%scalar_metric_names = {'Duration','Nspikes','NActiveNeurons','isi','distance_pre','distance_post'};
scalar_metric_names = {'Duration','Nspikes','NActiveNeurons','isi', 'Distance_to_pre','tail', 'amplitude_Pop'};

% --- Population metric: cell array, one cell per event -------------------
% Each cell contains a numeric vector of neuron IDs participating in that event.
neuron_population = cell(n_events, 1);
total_neurons = size(NetworkBurstPlayer.timeStamp,2);
is_participant = false(1, total_neurons);                                  % =true, if neuron i is a participating unit (primary or secondary events)
is_participant(NetworkBurstPlayer.ChannelsIntegrated(:,1)) = true;
for k = 1:n_events1
    aux = [];
    for i=1:total_neurons
        if is_participant(i) && numel(NetworkBurstPlayer.timeStamp{k,i})>0         % participating neuron i fires during event k
            aux = vertcat(i,aux);
        end
    end
    neuron_population{k} =  aux;
    clear aux
end
if n_events2 > MIN_SECONDARY_EVENTS
    for k = 1:n_events2
        aux = []; idx = n_events1 + k;
        for i=1:total_neurons
            if is_participant(i) && numel(NetworkBurstPlayer.Secondary.timeStamp{k,i})>0   % participating neuron i fires during event k
                aux = vertcat(i,aux);
            end
        end
        neuron_population{idx} =  aux;
        clear aux
    end
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
%  STEP 2 – Jaccard distance for neuron population metric
% =========================================================================
% Represent each event as a binary membership vector (1 x total_neurons)
pop_binary = false(n_events, total_neurons);
for i = 1:n_events
    pop_binary(i, neuron_population{i}) = true;
end

% Jaccard distance: 1 - |A∩B| / |A∪B|
D_population = pdist(double(pop_binary), 'jaccard');
D_population(isnan(D_population)) = 1;  % two empty sets → max distance

%% =========================================================================
%  STEP 3 – Composite distance matrix
% =========================================================================
D_composite = w_scalar * D_scalar_norm + w_population * D_population;

% Square form kept for reference (e.g. inspection, export)
D_mat = squareform(D_composite);

%% =========================================================================
%  STEP 4 – Project to Euclidean space via classical MDS
% =========================================================================
% Ward linkage minimises within-cluster variance, which is only well-defined
% in Euclidean space. Passing a composite (non-Euclidean) distance matrix
% directly to Ward produces incorrect merges and a MATLAB warning.
% Classical MDS (cmdscale) finds the best Euclidean embedding of D_mat;
% dimensions with negative eigenvalues reflect non-Euclidean structure
% (expected here due to Jaccard) and are discarded.
%
% Y_mds is used for both clustering (Step 5) and the PCoA scatter (Fig 4),
% so this single projection serves both purposes.

[Y_mds, eigvals] = cmdscale(D_mat);
Y_mds = Y_mds(:, eigvals > 1e-10);   % keep only valid Euclidean dimensions

n_neg = sum(eigvals < -1e-10);
if n_neg > 0 && option_plot
    fprintf(['Note: %d negative MDS eigenvalue(s) discarded ' ...
             '(expected for composite/non-Euclidean distances).\n\n'], n_neg);
end

% Variance explained by each retained MDS dimension (for Figure 4 axis labels)
var_explained = eigvals(eigvals > 1e-10) / sum(eigvals(eigvals > 0)) * 100;

%% =========================================================================
%  STEP 5 – Hierarchical clustering + silhouette-based model selection
% =========================================================================

[Z, method_used] = safe_linkage(Y_mds, LINKAGE_METHOD);
%Z = linkage(Y_mds, LINKAGE_METHOD);   % Ward on Euclidean MDS coordinates

if option_plot
    fprintf('\n%-10s  %-15s\n', 'K clusters', 'Mean silhouette');
    fprintf('%s\n', repmat('-', 1, 28));
end

sil_scores = zeros(1, MAX_K);
sil_scores(1) = NaN;  % silhouette undefined for k=1

for k = 2:MAX_K
    labels_k    = cluster(Z, 'maxclust', k);
    s           = silhouette(Y_mds, labels_k);  % Euclidean silhouette on MDS coords
    sil_scores(k) = mean(s);
    if option_plot
        fprintf('%-10d  %-15.4f\n', k, sil_scores(k));
    end
end

[~, best_k] = max(sil_scores);
if option_plot
    fprintf('\n>> Best K = %d  (mean silhouette = %.4f)\n\n', best_k, sil_scores(best_k));
end

event_labels = cluster(Z, 'maxclust', best_k);

%% =========================================================================
% Bootstrapping and permutation analysis
% =========================================================================

size_cluster = nan(best_k,1);
idx = unique(event_labels);
for i=1:best_k
    size_cluster(i) = numel(find(event_labels==idx(i)));
end
[min_cluster, aux] = min(size_cluster);
idx_min_cluster = idx(aux);

if min_cluster > MIN_SIZE_CLUSTER
    idx_min_cluster = [];
    bootstrap_k = bootstrap_optimal_k(total_neurons, n_events, scalar_metrics, neuron_population, w_scalar, w_population, LINKAGE_METHOD, MAX_K, idx_min_cluster, best_k, option_plot);
else
    bootstrap_k = bootstrap_optimal_k(total_neurons, n_events, scalar_metrics, neuron_population, w_scalar, w_population, LINKAGE_METHOD, MAX_K, idx_min_cluster, best_k, option_plot);
end

permutation_test = permutation_test_clustering(total_neurons, n_events, scalar_metrics,neuron_population, w_scalar,w_population, LINKAGE_METHOD, MAX_K, 0);

% Final cluster assignment

quality_clustering = 0;
k = bootstrap_k.optimal;
vote = bootstrap_k.vote_ratio(k);
vote0 = bootstrap_k.best_k0.vote_ratio;
si = bootstrap_k.silhouette(k);
si0 = bootstrap_k.best_k0.silhouette;

if n_events2 > 0
   alpha_sig = ALPHA;       
   sil_sig =  SIL_THRESHOLD;
else                            % make it more conservative when there are no secondary peaks
    alpha_sig = ALPHA / 4;    
    sil_sig =  SIL_THRESHOLD * 1.2;
end

optimal_k = 1;          % default value                                      
quality_clustering = 0; % default value           
if best_k>1
    if k == best_k                                                          % the optimal k value from boostrap is the original one  
        if permutation_test.p_values(best_k) < ALPHA                        % the result of the permutation is statistically significant   
            optimal_k = best_k;
            quality_clustering = 2;
        else                                                                % the result of the permutation is statistically non-significant
            aux_k = find(permutation_test.p_values(:) < alpha_sig);         % the result of the permutation is statistically significant 
            if numel(aux_k) > 0                                             % at least one k-size with statistically significant result after permutation
                idx = aux_k(1);                                             % smaller k-size statistically significant
                vote0 = bootstrap_k.vote_ratio(k);
                si0 = bootstrap_k.silhouette(k);
                for i=idx:MAX_K                                             % across all k-sizes                                                                          
                    vote = bootstrap_k.vote_ratio(i);                           
                    si = bootstrap_k.silhouette(i);
                    if abs((vote-vote0)/((vote+vote0)/2)) < MIN_DIFF || abs((si-si0)/((si+si0)/2)) < MIN_DIFF   % the robutness of the clustering is similar
                        optimal_k = i;                                      % the smallest k_size with similar robutness than the statistically significant one (idx)
                        quality_clustering = 1;
                        break
                    end
                end
            else                                                            % there is no k_size with statistically significant result after permutation
                optimal_k = 1;                                              % =1
                quality_clustering = 0;            
            end
        end
    else                                                                    % the optimal k-size from boostrap is different from the original one 
        if k>best_k                                                         % the optimal k-size from boostrap is larger than the original one     
           aux_k = find(permutation_test.p_values(:) < alpha_sig);          % the result of the permutation is statistically significant 
            if numel(aux_k) > 0                                             % at least one k-size with statistically significant result after permutation
                idx = aux_k(1);                                             % smaller k-size statistically significant
                vote0 = bootstrap_k.vote_ratio(idx);
                si0 = bootstrap_k.silhouette(idx);
%                for i=best_k:idx                                            % across all k-sizes   
                for i=1:k 
                    vote = bootstrap_k.vote_ratio(i);                           
                    si = bootstrap_k.silhouette(i);
                    if abs((vote-vote0)/((vote+vote0)/2)) < MIN_DIFF || abs((si-si0)/((si+si0)/2)) < MIN_DIFF   % the robutness of the clustering is similar
                        optimal_k = i;                                      % the smallest k_size with similar robutness than the statistically significant one (idx)
                        quality_clustering = 1;
                        break
                    end
                end
            else                                                            % there is no k_size with statistically significant result after permutation
                optimal_k = 1;                                              % =1
                quality_clustering = 0;            
            end
        else                                                                % the optimal k-size from boostrap is smaller than the original one
            if k == 1                                                       % the boostrap did not find any meaningful cluster
                if bootstrap_k.best_k0.silhouette > sil_sig && permutation_test.p_values(best_k) < alpha_sig
                    optimal_k = best_k;
                    quality_clustering = 1;
                else
                    optimal_k = k;
                    quality_clustering = 0;
                end
            else
                optimal_k = k;
                if permutation_test.p_values(k) < alpha_sig                            % the result of the permutation is statistically significant   
                    quality_clustering = 1;
                else                                                                % the result of the permutation is statistically non-significant
                    quality_clustering = 0;
                end
            end
        end
    end
else                                                                        % the original k-size is 1
    if k>1                                                                  % originally single cluster, but boostrapping found a better clustering
       aux_k = find(permutation_test.p_values(:) < alpha_sig);                  % the result of the permutation is statistically significant 
        if numel(aux_k) > 0                                                 % at least one k-size with statistically significant result after permutation
            idx = aux_k(1);                                                 % smaller k-size statistically significant
            vote0 = bootstrap_k.vote_ratio(idx);
            si0 = bootstrap_k.silhouette(idx);
            for i=2:idx                                                     % across all k-sizes                                                                          
                vote = bootstrap_k.vote_ratio(i);                           
                si = bootstrap_k.silhouette(i);
                if abs((vote-vote0)/((vote+vote0)/2)) < MIN_DIFF || abs((si-si0)/((si+si0)/2)) < MIN_DIFF   % the robutness of the clustering is similar
                    optimal_k = i;                                          % the smallest k_size with similar robutness than the statistically significant one (idx)
                    quality_clustering = 1;
                    break
                end
            end
        else                                                                % there is no k_size with statistically significant result after permutation
            optimal_k = best_k;                                             % =1
            quality_clustering = 1;            
        end
    else                                                                    % single cluster (both original and after boostrapping) (no need of permutation)
        optimal_k = best_k;                                                 % =1
        quality_clustering = 2;
    end
end
    
cluster_neuronal_events.n_cluster = optimal_k;  
label = nan(n_events,2);
label(:,1) = cluster(Z, 'maxclust', optimal_k);
label(1:n_events1, 2) = 1;
if n_events2>0
    label((n_events1+1):n_events, 2) = 2;
end
cluster_neuronal_events.label = label;
cluster_neuronal_events.clusteringQuality = quality_clustering;

%% =========================================================================
%  STEP 6 – Visualisation
% =========================================================================

if option_plot

cmap = lines(MAX_K);   % color palette

if 1
% ----- Figure 1: Dendrogram ----------------------------------------------
figure('Name','Dendrogram','Color','w','Position',[100 100 700 420]);
dendrogram(Z, 0, 'ColorThreshold', 'default');
title(sprintf('Hierarchical Clustering Dendrogram  (%s linkage)', method_used));
xlabel('Event index'); ylabel('Distance');
set(gca,'FontSize',11); box off;
end

% ----- Figure 2: Silhouette scores per K ---------------------------------
figure('Name','Silhouette scores','Color','w');
bar(2:MAX_K, sil_scores(2:end), 0.5, 'FaceColor',[0.35 0.60 0.85]);
hold on;
plot(best_k, sil_scores(best_k), 'r*', 'MarkerSize', 14, 'LineWidth', 2);
xlabel('Number of clusters K'); ylabel('Mean silhouette score');
title('Model selection: silhouette score vs K');
xticks(2:MAX_K); ylim([0 1]); grid on; box off;
set(gca,'FontSize',11);

if 1
% ----- Figure 3: Per-event silhouette plot (best K) ----------------------
s_best = silhouette(Y_mds, event_labels);
figure('Name',sprintf('Silhouette  K=%d', best_k),'Color','w', ...
       'Position',[100 560 500 420]);
[~, srt] = sort(event_labels);
barh(s_best(srt), 'FaceColor', 'flat', 'CData', cmap(event_labels(srt),:));
xlabel('Silhouette value'); ylabel('Events (sorted by cluster)');
title(sprintf('Silhouette plot  (K = %d)', best_k));
set(gca,'FontSize',11); box off;
end

% ----- Figure 4: PCoA scatter (first 2 MDS dimensions) ------------------
% Y_mds was already computed in Step 4 — no redundant cmdscale call needed.
figure('Name','PCoA scatter','Color','w');
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
title(sprintf('PCoA projection  –  K = %d clusters', best_k));
set(gca,'FontSize',11);

if 1
% ----- Figure 5: Scalar metric heatmap sorted by cluster -----------------
figure('Name','Metric heatmap','Color','w','Position',[100 100 820 480]);
[~, sort_idx] = sort(event_labels);
imagesc(scalar_norm(sort_idx, :));
colormap(redblue(256));
colorbar;
xlabel('Metric'); ylabel('Event (sorted by cluster)');
title('Normalized scalar metrics  –  sorted by cluster');
xticks(1:numel(scalar_metric_names));
xticklabels(scalar_metric_names);
xtickangle(30);
set(gca,'FontSize',10);

% Overlay cluster boundaries
boundaries = find(diff(event_labels(sort_idx))) + 0.5;
for b = boundaries'
    yline(b, 'w-', 'LineWidth', 2);
end
end

% ----- Figure 6: Raster plot of participating neurons and Population Firing -----------------
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

    for k=1:cluster_neuronal_events.n_cluster
        aux = find(cluster_neuronal_events.label(:,1)==k);
        for i=1:numel(aux)
            idx = aux(i);
            if idx>n_events1
                peak_time = NetworkBurstPlayer.Secondary.IndivBurst(idx-n_events1).peak;
                peak_Amp = NetworkBurstPlayer.Secondary.IndivBurst(idx-n_events1).peakAmp;
            else
                peak_time = NetworkBurstPlayer.IndivBurst(idx).peak;
                peak_Amp = NetworkBurstPlayer.IndivBurst(idx).peakAmp;
            end
            plot(peak_time, peak_Amp*1.1,'o', 'MarkerFaceColor',colors_event{k}); hold on
        end
    end
           
end


%% =========================================================================
%  STEP 7 – Print cluster summary to console
% =========================================================================
fprintf('%-8s  %-10s  %s\n', 'Event', 'Cluster');
fprintf('%s\n', repmat('-', 1, 55));
for i = 1:n_events
    fprintf('%-8d  %-10d \n', i, event_labels(i));
end

fprintf('\n--- Cluster means (scalar metrics) ---\n');
fprintf('%-12s', '', scalar_metric_names{:});
fprintf('\n');
for k = 1:best_k
    idx = event_labels == k;
    mu  = mean(scalar_metrics(idx,:), 1, 'omitnan');
    fprintf('Cluster %-3d  ', k);
    fprintf('%-12.3f', mu);
    fprintf('\n');
end

end
