%% =========================================================================
%  LOCAL FUNCTION – build composite distance matrix then project via MDS
% =========================================================================
function [Y_mds, var_explained, n_neg] = build_mds_v2(scalar_m, neuro_pop, ...
                                          n_ev, total_n, w_sc, w_pop)
    % --- Scalar part (normalised Euclidean) ---
    sc_norm = normalize(scalar_m);                     % (revised October 2026) as in "cluster_neuronal_events" (before: zscore, identical
                                                       %  without NaN; with a NaN, zscore turned the whole metric into NaN, i.e., 0 below)
    sc_norm(isnan(sc_norm)) = 0;
    D_sc = pdist(sc_norm, 'euclidean');
    D_sc = D_sc ./ max(D_sc + eps);

    % --- Population part (Jaccard) ---
    pop_bin = false(n_ev, total_n);
    for i = 1:n_ev
        if ~isempty(neuro_pop{i})
            pop_bin(i, neuro_pop{i}) = true;
        end
    end
    D_pop = pdist(double(pop_bin), 'jaccard');
    D_pop(isnan(D_pop)) = 1;

    % --- Composite distance ---
    D_comp = w_sc * D_sc + w_pop * D_pop;
    D_sq   = squareform(D_comp);

    % --- Project to Euclidean space via classical MDS ---
    [Y, eigvals] = cmdscale(D_sq);
    Y_mds = Y(:, eigvals > 1e-10);
    n_neg = sum(eigvals < -1e-10);

    if isempty(Y_mds)
        Y_mds = zeros(n_ev, 1);
        var_explained = [];
    else
        var_explained = eigvals(eigvals > 1e-10) / sum(eigvals(eigvals > 0)) * 100;
    end
end
