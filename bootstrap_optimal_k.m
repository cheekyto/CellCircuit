function bootstrap_optimal_k = bootstrap_optimal_k(total_neurons, n_events, scalar_metrics,neuron_population, w_scalar,w_population, LINKAGE_METHOD, MAX_K, idx_min_cluster, best_k, option_plot)

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
% METHOD:
%   1. Resample events WITH replacement (bootstrap sample, same n as original)
%   2. Build composite distance matrix on the resample
%   3. Project to Euclidean space via classical MDS (cmdscale) so that
%      Ward linkage is mathematically valid (Ward requires Euclidean space)
%   4. Refit clustering on MDS coordinates + record winning K
%   5. Summarise the K-vote distribution and apply the conservative decision rule
%
% NOTE ON MDS PER RESAMPLE:
%   Each bootstrap resample produces a different distance matrix, so each
%   gets its own independent MDS projection. This is correct — we are
%   testing whether the cluster structure is stable across resamples of the
%   data, not across projections of the same distances.

% NOTE ON LINKAGE FALLBACK:
%   Even after MDS projection, Ward linkage occasionally produces a
%   non-monotonic tree when the Euclidean embedding does not fully capture
%   the composite distance structure. safe_linkage() detects this at
%   runtime and automatically falls back to average linkage. The fallback
%   count across all bootstrap resamples is reported at the end.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Analysis settings ----------------------------------------------------
% (revised October 2026) fixed seed for the random numbers, so that the result does not change from one run to the next
% (before, it depended on the state of the global random number generator, i.e., on any random number drawn before in the MATLAB session);
% the global generator is restored when the function ends, so the rest of the code is not affected
RNG_SEED = 0;                                                               % seed for the random number generator (change it to get a different, but still reproducible, result)
N_BOOT         = 1000;   % number of bootstrap resamples (>= 1000 recommended)

% Conservative decision thresholds
VOTE_THRESHOLD = 0.60;   % candidate K must win in >= 60% of bootstraps (the higher the value, the more conservative)
SIL_THRESHOLD  = 0.30;   % AND have a mean silhouette >= 0.25 (0.25 is a modest but non-trivial structure threshold) (the higher the value, the more conservative)
apply_parameters('bootstrap_optimal_k');                                    % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")
rng_state = rng;                                                            % state of the global random number generator before this function
restore_rng = onCleanup(@() rng(rng_state));                                %#ok<NASGU> restore it at the end (also if there is an error)
rng(RNG_SEED, 'twister');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% =========================================================================
%  BUILD MDS COORDINATES ON ORIGINAL DATA
% =========================================================================
Y_orig = build_mds(scalar_metrics,neuron_population,n_events,total_neurons, w_scalar, w_population);

% Silhouette of the original data at each K (for reference)
%Z_orig   = linkage(Y_orig, LINKAGE_METHOD);
[Z_orig, method_orig] = safe_linkage(Y_orig, LINKAGE_METHOD);
%fprintf('Original data linkage method: %s\n', method_orig);

orig_sil = nan(1, MAX_K);
for k = 2:MAX_K
    lbl        = cluster(Z_orig, 'maxclust', k);
    orig_sil(k) = mean(silhouette(Y_orig, lbl));
end
[~, orig_best_k] = max(orig_sil);
if option_plot
    fprintf('Original data  –  best K = %d  (sil = %.4f)\n\n', ...
            orig_best_k, orig_sil(orig_best_k));
end

%% =========================================================================
%  BOOTSTRAP LOOP
% =========================================================================
boot_best_k  = zeros(N_BOOT, 1);   % winning K for each resample
boot_sil_mat = nan(N_BOOT, MAX_K); % silhouette at each K for each resample
n_fallbacks  = 0;

if option_plot
    fprintf('Running %d bootstrap resamples', N_BOOT);
end
    
for b = 1:N_BOOT
    if mod(b, 100) == 0, fprintf('.'); end

    % Resample with replacement
    idx  = randi(n_events, n_events, 1);
    if numel(idx_min_cluster)>0
        idx(1:numel(idx_min_cluster))=idx_min_cluster(:);
    end

    sm_b = scalar_metrics(idx, :);
    np_b = neuron_population(idx);

    % Build MDS coordinates for this resample (own projection)
    Y_b = build_mds(sm_b, np_b, n_events, total_neurons, w_scalar, w_population);

    % Skip degenerate resamples (all points identical → zero distances)
    if size(Y_b, 2) == 0 || max(pdist(Y_b)) < eps
        boot_best_k(b) = 1;
        continue;
    end

%    Z_b   = linkage(Y_b, LINKAGE_METHOD);
    [Z_b, method_b] = safe_linkage(Y_b, LINKAGE_METHOD);
    if ~strcmp(method_b, LINKAGE_METHOD), n_fallbacks = n_fallbacks + 1; end

    sil_b = nan(1, MAX_K);

    for k = 2:MAX_K
        lbl_b = cluster(Z_b, 'maxclust', k);
        % Guard: silhouette requires at least 2 distinct labels
        if numel(unique(lbl_b)) < 2
            sil_b(k) = 0;
        else
            sil_b(k) = mean(silhouette(Y_b, lbl_b));
        end
    end
    boot_sil_mat(b, 2:MAX_K) = sil_b(2:MAX_K);

    % Select best K among k>=2 only if it beats k=1 (sil=0 baseline)
    [best_sil_b, best_k_b] = max(sil_b);
    if isnan(best_sil_b) || best_sil_b <= 0
        boot_best_k(b) = 1;
    else
        boot_best_k(b) = best_k_b;
    end
end
if option_plot
    fprintf('\nDone.\n\n');
end
if n_fallbacks > 0
    fprintf('Note: fallback to average linkage occurred in %d / %d resamples.\n', ...
            n_fallbacks, N_BOOT);
    fprintf('\n');
end

%% =========================================================================
%  VOTE DISTRIBUTION AND CONSERVATIVE DECISION
% =========================================================================
vote_counts = histcounts(boot_best_k, 1:MAX_K+1);
vote_frac   = vote_counts / N_BOOT;

if option_plot
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
end

% Conservative decision rule: smallest K that passes BOTH criteria
chosen_k = 1;
for k = 2:MAX_K
    mu_sil = mean(boot_sil_mat(:,k), 'omitnan');
    if vote_frac(k) >= VOTE_THRESHOLD && mu_sil >= SIL_THRESHOLD
        chosen_k = k;
        break
    end
end

if option_plot
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
end

bootstrap_optimal_k.optimal = chosen_k;
bootstrap_optimal_k.vote_ratio = vote_frac;
bootstrap_optimal_k.silhouette  = mean(boot_sil_mat, 1, 'omitnan');

bootstrap_optimal_k.best_k0.value = best_k;
bootstrap_optimal_k.best_k0.vote_ratio = vote_frac(best_k);
bootstrap_optimal_k.best_k0.silhouette  = mean(boot_sil_mat(:,best_k), 'omitnan');

%% =========================================================================
%  VISUALISATION
% =========================================================================

if option_plot
% ----- Figure 1: K vote distribution (bar chart) -------------------------
figure('Name','Bootstrap K votes','Color','w');
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
title(sprintf('Bootstrap K stability  (N = %d resamples)', N_BOOT));
xticks(1:MAX_K);
xticklabels(arrayfun(@(k) sprintf('K = %d', k), 1:MAX_K, 'UniformOutput', false));
ylim([0, max(vote_frac)*100*1.15 + 5]);
text(chosen_k, vote_frac(chosen_k)*100 + 1.5, ...
     sprintf('Chosen\nK=%d', chosen_k), ...
     'HorizontalAlignment','center','FontSize',10,'Color',[0.1 0.45 0.25]);
set(gca,'FontSize',11); box off; grid on;

if 0
% ----- Figure 2: Silhouette distribution per K (box plots) ---------------
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
end

% ----- Figure 3: Cumulative vote curve -----------------------------------
figure('Name','Cumulative vote stability','Color','w');
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

end

end
