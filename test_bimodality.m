function results = test_bimodality(data)

% TEST_BIMODALITY  Assess whether a 1-D dataset is bimodal.
%
%   results = test_bimodality(data)
%

% (revised October 2026) fixed seed for the random numbers, so that the result does not change from one run to the next
% (before, it depended on the state of the global random number generator, i.e., on any random number drawn before in the MATLAB session);
% the global generator is restored when the function ends, so the rest of the code is not affected
RNG_SEED = 0;                                                               % seed for the random number generator (change it to get a different, but still reproducible, result)
THRESHOLD_BC = 0.555;          % (default: 0.555)
apply_parameters('test_bimodality');                                        % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")
rng_state = rng;                                                            % state of the global random number generator before this function
restore_rng = onCleanup(@() rng(rng_state));                                %#ok<NASGU> restore it at the end (also if there is an error)
rng(RNG_SEED, 'twister');

data = data(:);
n = numel(data);

%% --- 1. Bimodality coefficient ---
% BC > 0.555 (uniform distribution's value) suggests bimodality
s = skewness(data);
k = kurtosis(data); % note: MATLAB's kurtosis is NOT excess kurtosis
k_excess = k - 3;
BC = (s^2 + 1) / (k_excess + 3*(n-1)^2 / ((n-2)*(n-3)));

results.bimodality_coefficient = BC;
results.BC_suggests_bimodal = BC > THRESHOLD_BC;

%% --- 2. Gaussian Mixture Model comparison (1 vs 2 components) ---

%gm1 = fitgmdist(data, 1);
%gm2 = fitgmdist(data, 2, 'Replicates', 5, 'Options', statset('MaxIter', 500),'RegularizationValue', 1e-5);

gmOptions = statset('MaxIter', 1000, 'TolFun', 1e-6);

gm1 = fitgmdist(data, 1, 'Options', gmOptions, 'RegularizationValue', 1e-4);

% Suppress the specific convergence warning locally, but capture it
warnState = warning('off', 'stats:gmdistribution:FailedToConverge');
gm2 = fitgmdist(data, 2, 'Replicates', 10, 'Options', gmOptions, ...
    'RegularizationValue', 1e-4);
warning(warnState); % restore previous warning state

results.gm1_converged = gm1.Converged;
results.gm2_converged = gm2.Converged;

if ~gm2.Converged
    warning('test_bimodality:GMMNotConverged', ...
        '2-component GMM did not converge for this dataset (n=%d).', n);
end

results.AIC_1comp = gm1.AIC;
results.AIC_2comp = gm2.AIC;
results.BIC_1comp = gm1.BIC;
results.BIC_2comp = gm2.BIC;
results.GMM_suggests_bimodal = (gm2.AIC < gm1.AIC) && (gm2.BIC < gm1.BIC);
results.gm2 = gm2; % keep the fitted 2-component model (means, weights, etc.)

results.final = results.BC_suggests_bimodal && results.GMM_suggests_bimodal;

end