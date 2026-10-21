function Mutual_Information = Mutual_Information(neuron1, neuron2) 

% Function to calculate the mutual information between the activity of two neurons
% Estimates MI(X;Y) in bits from paired discrete observations.
% The joint distribution P(X,Y) is estimated from the empirical histogram of (counts1, counts2) pairs.  Marginals P(X) and P(Y) are obtained by summing the joint.  Zero-count bins are skipped (0·log2(0) ≡ 0).
%   MI is computed as:
%     MI(X;Y) = sum_x sum_y P(x,y) * log2( P(x,y) / (P(x)*P(y)) )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   neurons1:   array with number of spikes for each bin, neuron 1  
%   neurons2:   array with number of spikes for each bin, neuron 2 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Mutual_Information

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ensure column vectors
neuron1 = neuron1(:);
neuron2 = neuron2(:);

if any(neuron1 < 0) || any(neuron2 < 0)
    error('Spike counts must be non-negative integers!!!');
end

if numel(neuron1) ~= numel(neuron2)
    error('Spike-time cell arrays must contain the same number of bins/length');
end

% Joint frequency table (revised October 2026: built in a single pass with "accumarray", instead of comparing every pair 
% of values across all bins; the result only differs by floating-point round-off)
N     = numel(neuron1);             % total number of observations (bins)
if all(neuron1 == round(neuron1)) && all(neuron2 == round(neuron2))    % integer spike counts (the usual case): the count itself is the index (count+1)
    ix = neuron1 + 1;   iy = neuron2 + 1;                   % (the values that do not occur give empty rows/columns, which do not contribute to MI)
else                                                        % other values: index of each value within the sorted list of unique values
    [~, ~, ix] = unique(neuron1);
    [~, ~, iy] = unique(neuron2);
end
jointCount = accumarray([ix iy], 1);

% Estimate probabilities
Pxy = jointCount / N;             % joint
Px  = sum(Pxy, 2);                % marginal X
Py  = sum(Pxy, 1);                % marginal Y

% Calculate Mutual Information in bits (only the non-zero joint probabilities; the marginals are then also > 0)
PxPy = Px * Py;                   % product of the marginals, for every pair of values
nz = Pxy > 0;
MI = sum(Pxy(nz) .* log2(Pxy(nz) ./ PxPy(nz)));
if N == 0                         % no bins: same result as before (zero)
    MI = 0;
end

Mutual_Information = MI;

end

