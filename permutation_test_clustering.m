function permutation_test_clustering = permutation_test_clustering(total_neurons, n_events, scalar_metrics,neuron_population, w_scalar,w_population, LINKAGE_METHOD, MAX_K, option_plot)

%% permutation_test_clustering.m
%
% Tests whether the clustering structure found in your neuronal event data
% is statistically significant, or could have arisen by chance.
%
% NULL HYPOTHESIS: there is no genuine joint structure linking the scalar
%                  metrics (including their correlations with each other)
%                  to the neuron population identities — any apparent
%                  clustering is as good as pairing randomly-reshuffled
%                  scalar metrics with a randomly-assigned population
%
% IMPORTANT NOTE ON METHOD (please read):
%   Clustering algorithms (Ward linkage, silhouette) depend ONLY on the
%   pairwise distances between events, NOT on the order/labelling of the
%   events. Permuting the row order of a distance matrix or coordinate set
%   is therefore a NO-OP: it yields the exact same set of pairwise
%   distances, just relabelled, so the resulting silhouette score is
%   identical on every permutation (this was a bug in an earlier version
%   of this script, visible as perm_sd = 0 and perm_mu = real_sil).
%
%   This version uses a STRONG null: each of the 8 scalar-metric COLUMNS is
%   independently permuted across events (so metric 1 gets one random
%   shuffle, metric 2 gets a different random shuffle, etc.), fabricating a
%   surrogate scalar-metric dataset that destroys both (a) the correlation
%   structure among the scalar metrics themselves, and (b) any association
%   between the scalar metrics and the neuron population identity. This
%   surrogate dataset is then randomly paired with the (unshuffled) set of
%   population vectors. The composite distance + MDS projection is rebuilt
%   from scratch on this fabricated data at each permutation.
%
% METHOD:
%   1. Build composite distance + MDS projection on the REAL data; compute
%      real silhouette score for K = 2..MAX_K
%   2. Repeat N_PERM times:
%        - Independently permute each scalar-metric column across events
%        - Independently permute the order of population vectors
%        - Rebuild composite distance + MDS projection from this fabricated
%          pairing
%        - Refit clustering and record silhouette score
%   3. p-value = fraction of permuted scores >= real score
%   4. Visualise the null distribution vs. the real score

% NOTE ON LINKAGE FALLBACK:
%   Even after MDS projection, Ward linkage occasionally produces a
%   non-monotonic tree when the Euclidean embedding does not fully capture
%   the composite distance structure. safe_linkage() detects this at
%   runtime and automatically falls back to average linkage. Both the real
%   data and every permutation use the same safe_linkage() call, so the
%   null distribution and the real score are always computed under the same
%   method.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Analysis settings ----------------------------------------------------
% (revised October 2026) fixed seed for the random numbers, so that the result does not change from one run to the next
% (before, it depended on the state of the global random number generator, i.e., on any random number drawn before in the MATLAB session);
% the global generator is restored when the function ends, so the rest of the code is not affected
RNG_SEED = 0;                                                               % seed for the random number generator (change it to get a different, but still reproducible, result)
N_PERM         = 1000;  % number of permutations (>= 1000 recommended)
ALPHA          = 0.05;  % significance threshold
apply_parameters('permutation_test_clustering');                            % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")
rng_state = rng;                                                            % state of the global random number generator before this function
restore_rng = onCleanup(@() rng(rng_state));                                %#ok<NASGU> restore it at the end (also if there is an error)
rng(RNG_SEED, 'twister');

%% =========================================================================
%  REAL DATA – build composite distance + MDS projection + silhouette
% =========================================================================
[Y_mds, ~, n_neg] = build_mds_v2(scalar_metrics, neuron_population, ...
                               n_events, total_neurons, w_scalar, w_population);

if n_neg > 0 && option_plot
    fprintf(['Note: %d negative MDS eigenvalue(s) discarded.\n' ...
             'This is normal for composite/non-Euclidean distances.\n\n'], n_neg);
end

real_sil = nan(1, MAX_K);   % silhouette undefined for k=1
%Z_real   = linkage(Y_mds, LINKAGE_METHOD);
[Z_real, method_real] = safe_linkage(Y_mds, LINKAGE_METHOD);
for k = 2:MAX_K
    lbl         = cluster(Z_real, 'maxclust', k);
    real_sil(k) = mean(silhouette(Y_mds, lbl));
end

[best_real_sil, best_k] = max(real_sil);
if option_plot
fprintf('Real data   –  best K = %d,  silhouette = %.4f\n\n', ...
        best_k, best_real_sil);
end

%% =========================================================================
%  PERMUTATION TEST
% =========================================================================
% At each permutation:
%   - Each scalar-metric column is independently shuffled across events,
%     fabricating a surrogate scalar-metric dataset with no real
%     event-to-event correspondence (and no preserved inter-metric
%     correlation structure).
%   - The order of population vectors is independently shuffled.
%   - The two surrogate datasets are paired up, and the full composite
%     distance + MDS pipeline is rebuilt from scratch.

perm_sil = nan(N_PERM, MAX_K);
n_metrics = size(scalar_metrics, 2);

if option_plot
fprintf('Running %d permutations', N_PERM);
end

% Track how often the fallback was triggered across permutations
n_fallbacks = 0;

for p = 1:N_PERM
    if mod(p, 100) == 0, fprintf('.'); end

    % Independently permute each scalar-metric column
    scalar_perm = scalar_metrics;
    for m = 1:n_metrics
        scalar_perm(:, m) = scalar_metrics(randperm(n_events), m);
    end

    % Independently permute the population vectors
    pop_perm = neuron_population(randperm(n_events));

    % Rebuild composite distance + MDS projection on the fabricated pairing
    [Y_perm, ~, ~] = build_mds_v2(scalar_perm, pop_perm, ...
                                n_events, total_neurons, w_scalar, w_population);

    if size(Y_perm, 2) == 0
        continue;   % degenerate case: skip (leaves NaN for this permutation)
    end

    [Z_perm, method_perm] = safe_linkage(Y_perm, LINKAGE_METHOD);
%    Z_perm = linkage(Y_perm, LINKAGE_METHOD);
    if ~strcmp(method_perm, LINKAGE_METHOD), n_fallbacks = n_fallbacks + 1; end

    for k = 2:MAX_K
        lbl_perm       = cluster(Z_perm, 'maxclust', k);
        perm_sil(p, k) = mean(silhouette(Y_perm, lbl_perm));
    end
end
if option_plot
fprintf('\nDone.\n\n');

%% =========================================================================
%  COMPUTE P-VALUES
% =========================================================================
fprintf('%-10s  %-16s  %-16s  %-10s  %s\n', ...
        'K', 'Real sil', 'Perm mean+-SD', 'p-value', 'Significant?');
fprintf('%s\n', repmat('-', 1, 68));
end

p_values = nan(1, MAX_K);
perm_mu = nan(1, MAX_K);
for k = 2:MAX_K
    p_val        = mean(perm_sil(:,k) >= real_sil(k), 'omitnan');
    p_values(k)  = p_val;
    perm_mu(k)   = mean(perm_sil(:,k), 'omitnan');
    perm_sd      = std(perm_sil(:,k), 'omitnan');
    sig_str      = '';
    if p_val < ALPHA,  sig_str = '*';   end
    if p_val < 0.01,   sig_str = '**';  end
    if p_val < 0.001,  sig_str = '***'; end
    if option_plot
    fprintf('%-10d  %-16.4f  %-16s  %-10.4f  %s\n', ...
            k, real_sil(k), ...
            sprintf('%.4f +- %.4f', perm_mu(k), perm_sd), ...
            p_val, sig_str);
    end
end
if option_plot
fprintf('\n(* p<0.05  ** p<0.01  *** p<0.001)\n\n');
end

% Overall interpretation
sig_ks = find(p_values < ALPHA);
if option_plot
if isempty(sig_ks)
    fprintf(['CONCLUSION: No K reached significance (all p >= %.2f).\n' ...
             'The clustering structure is not distinguishable from chance.\n'], ALPHA);
else
    fprintf(['CONCLUSION: K = %s reached significance.\n' ...
             'Cluster structure is unlikely to be due to chance.\n'], ...
             num2str(sig_ks, '%d '));
end
end
permutation_test_clustering.p_values = p_values;
permutation_test_clustering.silhouette = perm_mu;

%% =========================================================================
%  VISUALISATION
% =========================================================================
if option_plot

for k = 2:MAX_K
    figure('Name', sprintf('Permutation test  K=%d', k),'Color','w');

    % Null distribution
    histogram(perm_sil(:,k), 30, ...
              'FaceColor', [0.65 0.65 0.65], 'EdgeColor','none', ...
              'Normalization','probability');
    hold on;

    % Real score
    xline(real_sil(k), 'r-', 'LineWidth', 2.5, ...
          'Label', sprintf('Real  (%.3f)', real_sil(k)), ...
          'LabelVerticalAlignment','bottom', 'FontSize', 10);

    % Alpha threshold (95th percentile of null)
    thresh = prctile(perm_sil(:,k), (1-ALPHA)*100);
    xline(thresh, 'b--', 'LineWidth', 1.5, ...
          'Label', sprintf('alpha=%.2f threshold (%.3f)', ALPHA, thresh), ...
          'LabelVerticalAlignment','top', 'FontSize', 9);

    % Shade the rejection region
    xl = xlim;
    yl = ylim;
    patch([thresh xl(2) xl(2) thresh], [0 0 yl(2) yl(2)], ...
          [1 0.8 0.8], 'FaceAlpha', 0.25, 'EdgeColor','none');

    xlabel('Silhouette score (permuted)');
    ylabel('Probability');
    title(sprintf('Null distribution  –  K = %d   (p = %.4f)', k, p_values(k)));
    legend({'Null distribution','Real score', ...
            sprintf('alpha=%.2f threshold', ALPHA), 'Rejection region'}, ...
           'Location','northwest');
    set(gca,'FontSize',11); box off;
end
end

end