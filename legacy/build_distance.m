%% =========================================================================
%  LOCAL FUNCTION – build composite distance matrix
% =========================================================================
function [D_composite, D_square] = build_distance(scalar_metrics, neuron_population, n_events, total_neurons, w_scalar,w_population)
    
%% =========================================================================
%  STEP 1 – Normalize scalar metrics
% =========================================================================
% Z-score each column; replace NaN with 0 (handles constant columns)
scalar_norm = zscore(scalar_metrics);
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
    if ~isempty(neuron_population{i})
        pop_binary(i, neuron_population{i}) = true;
    end
end

% Jaccard distance: 1 - |A∩B| / |A∪B|
D_population = pdist(double(pop_binary), 'jaccard');
D_population(isnan(D_population)) = 1;  % two empty sets → max distance

%% =========================================================================
%  STEP 3 – Composite distance matrix
% =========================================================================
D_composite = w_scalar * D_scalar_norm + w_population * D_population;

% Convert to full square matrix for later use
D_square = squareform(D_composite);

end
