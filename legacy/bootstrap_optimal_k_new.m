%% bootstrap_optimal_k.m
%
% Tests the stability of the optimal number of clusters (K) across
% bootstrap resamples of your neuronal event data.
%
% CONSERVATIVE PHILOSOPHY:
%   Default to K = 1 (no clustering) unless there is clear, consistent
%   evidence for K >= 2. Two criteria must BOTH be met to accept K > 1:
%     (1) The candidate K wins in >= VOTE_THRESHOLD fraction of bootstraps
%     (2) The mean silhouette at that K is >= SIL_THRESHOLD
%
% NOTE ON LINKAGE FALLBACK:
%   Even after MDS projection, Ward linkage occasionally produces a
%   non-monotonic tree when the Euclidean embedding does not fully capture
%   the composite distance structure. safe_linkage() detects this at
%   runtime and automatically falls back to average linkage. The fallback
%   count across all bootstrap resamples is reported at the end.
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
N_BOOT         = 1000;

VOTE_THRESHOLD = 0.60;
SIL_THRESHOLD  = 0.25;

%% =========================================================================
%  ORIGINAL DATA – MDS + clustering (for reference)
% =========================================================================
[Y_orig, ~, n_neg] = build_mds(scalar_metrics, neuron_population, ...
                                n_events, total_neurons, w_scalar, w_population);
if n_neg > 0
    fprintf(['Note: %d negative MDS eigenvalue(s) discarded.\n' ...
             '(Expected for composite/non-Euclidean distances.)\n\n'], n_neg);
end

[Z_orig, method_orig] = safe_linkage(Y_orig, LINKAGE_METHOD);
fprintf('Original data linkage method: %s\n', method_orig);

orig_sil = nan(1, MAX_K);
for k = 2:MAX_K
    lbl         = cluster(Z_orig, 'maxclust', k);
    orig_sil(k) = mean(silhouette(Y_orig, lbl));
end
[~, orig_best_k] = max(orig_sil);
fprintf('Original data  –  best K = %d  (sil = %.4f)\n\n', ...
        orig_best_k, orig_sil(orig_best_k));

%% =========================================================================
%  BOOTSTRAP LOOP
% =========================================================================
boot_best_k  = zeros(N_BOOT, 1);
boot_sil_mat = nan(N_BOOT, MAX_K);
n_fallbacks  = 0;

fprintf('Running %d bootstrap resamples', N_BOOT);
for b = 1:N_BOOT
    if mod(b, 100) == 0, fprintf('.'); end

    idx  = randi(n_events, n_events, 1);
    sm_b = scalar_metrics(idx, :);
    np_b = neuron_population(idx);

    [Y_b, ~, ~] = build_mds(sm_b, np_b, n_events, total_neurons, ...
                              w_scalar, w_population);

    if size(Y_b, 2) == 0 || max(pdist(Y_b)) < eps
        boot_best_k(b) = 1;
        continue;
    end

    [Z_b, method_b] = safe_linkage(Y_b, LINKAGE_METHOD);
    if ~strcmp(method_b, LINKAGE_METHOD), n_fallbacks = n_fallbacks + 1; end

    sil_b = nan(1, MAX_K);
    for k = 2:MAX_K
        lbl_b = cluster(Z_b, 'maxclust', k);
        if numel(unique(lbl_b)) < 2
            sil_b(k) = 0;
        else
            sil_b(k) = mean(silhouette(Y_b, lbl_b));
        end
    end
    boot_sil_mat(b, 2:MAX_K) = sil_b(2:MAX_K);

    [best_sil_b, best_k_b] = max(sil_b);
    if isnan(best_sil_b) || best_sil_b <= 0
        boot_best_k(b) = 1;
    else
        boot_best_k(b) = best_k_b;
    end
end
fprintf('\nDone.\n');
if n_fallbacks > 0
    fprintf('Note: fallback to average linkage occurred in %d / %d resamples.\n', ...
            n_fallbacks, N_BOOT);
end
fprintf('\n');

%% =========================================================================
%  VOTE DISTRIBUTION AND CONSERVATIVE DECISION
% =========================================================================
vote_counts = histcounts(boot_best_k, 1:MAX_K+1);
vote_frac   = vote_counts / N_BOOT;

fprintf('%-10s  %-14s  %-18s\n', 'K', 'Vote fraction', 'Mean sil (boots)');
fprintf('%s\n', repmat('-',1,46));
for k = 1:MAX_K
    mu_sil = mean(boot_sil_mat(:,k), 'omitnan');
    if k == 1
        fprintf('%-10d  %-14.3f  %-18s\n', k, vote_frac(k), 'N/A (baseline)');
    else
        fprintf('%-10d  %-14.3f  %-18.4f\n', k, vote_frac(k), mu_sil);
    end
end
fprintf('\n');

chosen_k = 1;
for k = 2:MAX_K
    mu_sil = mean(boot_sil_mat(:,k), 'omitnan');
    if vote_frac(k) >= VOTE_THRESHOLD && mu_sil >= SIL_THRESHOLD
        chosen_k = k;
        break
    end
end

fprintf('=== CONSERVATIVE DECISION ===\n');
if chosen_k == 1
    fprintf(['Chosen K = 1  (no clustering).\n' ...
             'No candidate K met both criteria:\n' ...
             '  vote >= %.0f%%  AND  mean silhouette >= %.2f\n'], ...
             VOTE_THRESHOLD*100, SIL_THRESHOLD);
else
    fprintf(['Chosen K = %d.\n' ...
             'This K won %.1f%% of bootstrap votes and had mean sil = %.4f,\n' ...
             'both exceeding the thresholds (%.0f%% votes, sil >= %.2f).\n'], ...
             chosen_k, vote_frac(chosen_k)*100, ...
             mean(boot_sil_mat(:,chosen_k),'omitnan'), ...
             VOTE_THRESHOLD*100, SIL_THRESHOLD);
end

%% =========================================================================
%  VISUALISATION
% =========================================================================

% ----- Figure 1: K vote distribution -------------------------------------
figure('Name','Bootstrap K votes','Color','w','Position',[80 100 520 380]);
b_colors = repmat([0.75 0.75 0.75], MAX_K, 1);
b_colors(chosen_k, :) = [0.20 0.60 0.35];
bh = bar(1:MAX_K, vote_frac*100, 0.55, 'FaceColor','flat');
bh.CData = b_colors;
hold on;
yline(VOTE_THRESHOLD*100, 'r--', 'LineWidth', 1.8, ...
      'Label', sprintf('Vote threshold (%.0f%%)', VOTE_THRESHOLD*100), ...
      'LabelHorizontalAlignment','right', 'FontSize', 9);
xlabel('Number of clusters K');
ylabel('Bootstrap vote (%)');
title(sprintf('Bootstrap K stability  (N = %d,  %s linkage)', N_BOOT, method_orig));
xticks(1:MAX_K);
xticklabels(arrayfun(@(k) sprintf('K = %d', k), 1:MAX_K, 'UniformOutput', false));
ylim([0, max(vote_frac)*100*1.15 + 5]);
text(chosen_k, vote_frac(chosen_k)*100 + 1.5, ...
     sprintf('Chosen\nK=%d', chosen_k), ...
     'HorizontalAlignment','center','FontSize',10,'Color',[0.1 0.45 0.25]);
set(gca,'FontSize',11); box off; grid on;

% ----- Figure 2: Silhouette distribution per K ---------------------------
figure('Name','Bootstrap silhouette distribution','Color','w', ...
       'Position',[620 100 520 380]);
sil_data = boot_sil_mat(:, 2:MAX_K);
boxplot(sil_data, 'Labels', arrayfun(@(k) sprintf('K=%d',k), 2:MAX_K, ...
        'UniformOutput', false), 'Whisker', 1.5, 'Symbol', 'o');
hold on;
yline(SIL_THRESHOLD, 'r--', 'LineWidth', 1.8, ...
      'Label', sprintf('Sil threshold (%.2f)', SIL_THRESHOLD), ...
      'LabelHorizontalAlignment','right','FontSize',9);
mu_sils = mean(sil_data, 'omitnan');
plot(1:size(sil_data,2), mu_sils, 'rd', 'MarkerFaceColor','r','MarkerSize',8);
ylabel('Silhouette score');
xlabel('Number of clusters K');
title('Silhouette score across bootstrap resamples');
set(gca,'FontSize',11); box off; grid on;

% ----- Figure 3: Cumulative vote curve -----------------------------------
figure('Name','Cumulative vote stability','Color','w','Position',[80 520 520 340]);
running_vote = zeros(N_BOOT, MAX_K);
for b = 1:N_BOOT
    for k = 1:MAX_K
        running_vote(b,k) = sum(boot_best_k(1:b) == k) / b;
    end
end
cmap = lines(MAX_K);
hold on;
for k = 1:MAX_K
    plot(1:N_BOOT, running_vote(:,k)*100, 'Color', cmap(k,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('K = %d', k));
end
yline(VOTE_THRESHOLD*100, 'k--', 'LineWidth', 1.2, ...
      'Label', sprintf('%.0f%% threshold', VOTE_THRESHOLD*100), ...
      'LabelHorizontalAlignment','right','FontSize',9);
xlabel('Number of bootstrap resamples');
ylabel('Running vote fraction (%)');
title('Convergence of bootstrap votes');
legend('Location','east'); set(gca,'FontSize',11); box off; grid on;

%% =========================================================================
%  LOCAL FUNCTION – safe linkage with automatic fallback
% =========================================================================
function [Z, method_used] = safe_linkage(Y, preferred_method)
% Attempts linkage with preferred_method. If MATLAB raises a non-monotonic
% tree warning, falls back to 'average' linkage automatically.
% Returns Z and method_used so callers can report which was applied.

    FALLBACK = 'average';
    warnID   = 'stats:linkage:NonMonotonicTree';
    prev     = warning('error', warnID);

    try
        Z           = linkage(Y, preferred_method);
        method_used = preferred_method;
    catch
        warning(prev);
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
