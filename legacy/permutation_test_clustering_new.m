%% permutation_test_clustering.m
%
% Tests whether the clustering structure found in your neuronal event data
% is statistically significant, or could have arisen by chance.
%
% NULL HYPOTHESIS: there is no genuine joint structure linking the scalar
%                  metrics (including their correlations with each other)
%                  to the neuron population identities — any apparent
%                  clustering is as good as pairing randomly-reshuffled
%                  scalar metrics with a randomly-assigned population.
%
% PERMUTATION LOGIC:
%   Each scalar-metric column is independently shuffled across events
%   (metric 1 gets one random permutation, metric 2 a different one, etc.),
%   fabricating surrogate scalar data with no preserved inter-metric
%   correlation structure. The population vectors are independently shuffled
%   too. The composite distance + MDS projection is rebuilt from scratch at
%   every permutation.
%
% NOTE ON LINKAGE FALLBACK:
%   Even after MDS projection, Ward linkage occasionally produces a
%   non-monotonic tree when the Euclidean embedding does not fully capture
%   the composite distance structure. safe_linkage() detects this at
%   runtime and automatically falls back to average linkage. Both the real
%   data and every permutation use the same safe_linkage() call, so the
%   null distribution and the real score are always computed under the same
%   method.
%
% -------------------------------------------------------------------------
% INPUTS  –– edit the "USER INPUT" section below
% -------------------------------------------------------------------------

clear; close all; clc;

%% =========================================================================
%  USER INPUT
% =========================================================================

n_events      = 20;
total_neurons = 30;

rng(42);
scalar_metrics = [rand(n_events,1)*500, rand(n_events,1)*2, ...
                  rand(n_events,1)*50,  rand(n_events,1),   ...
                  rand(n_events,1)*10,  randn(n_events,1),  ...
                  rand(n_events,1)*100, rand(n_events,1)*5];

neuron_population = cell(n_events, 1);
for i = 1:n_events
    k = randi([3, 15]);
    neuron_population{i} = sort(randperm(total_neurons, k));
end

w_scalar       = 0.6;
w_population   = 0.4;
MAX_K          = 3;
LINKAGE_METHOD = 'ward';
N_PERM         = 1000;
ALPHA          = 0.05;

%% =========================================================================
%  REAL DATA – composite distance + MDS + silhouette
% =========================================================================
[Y_mds, ~, n_neg] = build_mds(scalar_metrics, neuron_population, ...
                               n_events, total_neurons, w_scalar, w_population);
if n_neg > 0
    fprintf(['Note: %d negative MDS eigenvalue(s) discarded.\n' ...
             '(Expected for composite/non-Euclidean distances.)\n\n'], n_neg);
end

[Z_real, method_real] = safe_linkage(Y_mds, LINKAGE_METHOD);
fprintf('Real data linkage method: %s\n', method_real);

real_sil = nan(1, MAX_K);
for k = 2:MAX_K
    lbl         = cluster(Z_real, 'maxclust', k);
    real_sil(k) = mean(silhouette(Y_mds, lbl));
end
[best_real_sil, best_k] = max(real_sil);
fprintf('Real data  –  best K = %d,  silhouette = %.4f\n\n', ...
        best_k, best_real_sil);

%% =========================================================================
%  PERMUTATION TEST
% =========================================================================
perm_sil  = nan(N_PERM, MAX_K);
n_metrics = size(scalar_metrics, 2);

% Track how often the fallback was triggered across permutations
n_fallbacks = 0;

fprintf('Running %d permutations', N_PERM);
for p = 1:N_PERM
    if mod(p, 100) == 0, fprintf('.'); end

    % Independently permute each scalar-metric column
    scalar_perm = scalar_metrics;
    for m = 1:n_metrics
        scalar_perm(:, m) = scalar_metrics(randperm(n_events), m);
    end

    % Independently permute the population vectors
    pop_perm = neuron_population(randperm(n_events));

    % Rebuild composite distance + MDS on the fabricated pairing
    [Y_perm, ~, ~] = build_mds(scalar_perm, pop_perm, ...
                                n_events, total_neurons, w_scalar, w_population);
    if size(Y_perm, 2) == 0, continue; end

    [Z_perm, method_perm] = safe_linkage(Y_perm, LINKAGE_METHOD);
    if ~strcmp(method_perm, LINKAGE_METHOD), n_fallbacks = n_fallbacks + 1; end

    for k = 2:MAX_K
        lbl_perm       = cluster(Z_perm, 'maxclust', k);
        perm_sil(p, k) = mean(silhouette(Y_perm, lbl_perm));
    end
end
fprintf('\nDone.\n');
if n_fallbacks > 0
    fprintf('Note: fallback to average linkage occurred in %d / %d permutations.\n', ...
            n_fallbacks, N_PERM);
end
fprintf('\n');

%% =========================================================================
%  COMPUTE P-VALUES
% =========================================================================
fprintf('%-10s  %-16s  %-20s  %-10s  %s\n', ...
        'K', 'Real sil', 'Perm mean+-SD', 'p-value', 'Significant?');
fprintf('%s\n', repmat('-', 1, 72));

p_values = nan(1, MAX_K);
for k = 2:MAX_K
    p_val       = mean(perm_sil(:,k) >= real_sil(k), 'omitnan');
    p_values(k) = p_val;
    perm_mu     = mean(perm_sil(:,k), 'omitnan');
    perm_sd     = std(perm_sil(:,k), 'omitnan');
    sig_str = '';
    if p_val < ALPHA,  sig_str = '*';   end
    if p_val < 0.01,   sig_str = '**';  end
    if p_val < 0.001,  sig_str = '***'; end
    fprintf('%-10d  %-16.4f  %-20s  %-10.4f  %s\n', ...
            k, real_sil(k), sprintf('%.4f +- %.4f', perm_mu, perm_sd), ...
            p_val, sig_str);
end
fprintf('\n(* p<0.05  ** p<0.01  *** p<0.001)\n\n');

sig_ks = find(p_values < ALPHA);
if isempty(sig_ks)
    fprintf(['CONCLUSION: No K reached significance (all p >= %.2f).\n' ...
             'Clustering structure is not distinguishable from chance.\n'], ALPHA);
else
    fprintf(['CONCLUSION: K = %s reached significance.\n' ...
             'Cluster structure is unlikely to be due to chance.\n'], ...
             num2str(sig_ks, '%d '));
end

%% =========================================================================
%  VISUALISATION
% =========================================================================
for k = 2:MAX_K
    figure('Name', sprintf('Permutation test  K=%d', k), ...
           'Color','w','Position',[80+(k-2)*540, 100, 500, 380]);

    histogram(perm_sil(:,k), 30, ...
              'FaceColor',[0.65 0.65 0.65], 'EdgeColor','none', ...
              'Normalization','probability');
    hold on;

    xline(real_sil(k), 'r-', 'LineWidth', 2.5, ...
          'Label', sprintf('Real  (%.3f)', real_sil(k)), ...
          'LabelVerticalAlignment','bottom', 'FontSize', 10);

    thresh = prctile(perm_sil(:,k), (1-ALPHA)*100);
    xline(thresh, 'b--', 'LineWidth', 1.5, ...
          'Label', sprintf('alpha=%.2f threshold (%.3f)', ALPHA, thresh), ...
          'LabelVerticalAlignment','top', 'FontSize', 9);

    xl = xlim; yl = ylim;
    patch([thresh xl(2) xl(2) thresh], [0 0 yl(2) yl(2)], ...
          [1 0.8 0.8], 'FaceAlpha', 0.25, 'EdgeColor','none');

    xlabel('Silhouette score (permuted)');
    ylabel('Probability');
    title(sprintf('Null distribution  –  K = %d   (p = %.4f,  %s linkage)', ...
                  k, p_values(k), method_real));
    legend({'Null distribution','Real score', ...
            sprintf('alpha=%.2f threshold', ALPHA),'Rejection region'}, ...
           'Location','northwest');
    set(gca,'FontSize',11); box off;
end

%% =========================================================================
%  LOCAL FUNCTION – safe linkage with automatic fallback
% =========================================================================
function [Z, method_used] = safe_linkage(Y, preferred_method)
% Attempts linkage with preferred_method. If MATLAB raises a non-monotonic
% tree warning (Ward is inappropriate for this embedding), falls back to
% 'average' linkage automatically and prints a one-line notice.
% Returns Z and method_used so callers can report which was applied.

    FALLBACK = 'average';
    warnID   = 'stats:linkage:NonMonotonicTree';
    prev     = warning('error', warnID);

    try
        Z           = linkage(Y, preferred_method);
        method_used = preferred_method;
    catch
        warning(prev);
        fprintf(['  [safe_linkage] Non-monotonic tree detected: ' ...
                 'falling back from ''%s'' to ''%s''.\n'], ...
                 preferred_method, FALLBACK);
        Z           = linkage(Y, FALLBACK);
        method_used = FALLBACK;
    end

    warning(prev);
end

%% =========================================================================
%  LOCAL FUNCTION – build composite distance + MDS projection
% =========================================================================
function [Y_mds, var_explained, n_neg] = build_mds(scalar_m, neuro_pop, ...
                                          n_ev, total_n, w_sc, w_pop)
    sc_norm = zscore(scalar_m);
    sc_norm(isnan(sc_norm)) = 0;
    D_sc = pdist(sc_norm, 'euclidean');
    D_sc = D_sc ./ max(D_sc + eps);

    pop_bin = false(n_ev, total_n);
    for i = 1:n_ev
        if ~isempty(neuro_pop{i})
            pop_bin(i, neuro_pop{i}) = true;
        end
    end
    D_pop = pdist(pop_bin, 'jaccard');
    D_pop(isnan(D_pop)) = 1;

    D_comp = w_sc * D_sc + w_pop * D_pop;
    D_sq   = squareform(D_comp);

    [Y, eigvals] = cmdscale(D_sq);
    Y_mds  = Y(:, eigvals > 1e-10);
    n_neg  = sum(eigvals < -1e-10);

    if isempty(Y_mds)
        Y_mds = zeros(n_ev, 1);
        var_explained = [];
    else
        var_explained = eigvals(eigvals > 1e-10) / sum(eigvals(eigvals > 0)) * 100;
    end
end
