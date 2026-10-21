function Entropy = Entropy(spikeTrain)

%   Function to calculate the Entropy of a spike train

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   spikeTrain:   array with number of spikes per bin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Entropy:    entropy in bits of the count-valued random variable per bin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Notes:
%      - Each bin contains an integer spike count: 0,1,2,...
%      - The entropy is computed from the empirical distribution of counts.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ensure column vector
spikeTrain = spikeTrain(:);

% possible count values (alphabet) observed, and empirical probability distribution p(k)
% (revised October 2026: counted in a single pass with "accumarray", instead of one comparison per value; 
%  the result only differs by floating-point round-off)
N = numel(spikeTrain);                          % number of bins
if all(spikeTrain == round(spikeTrain))         % integer spike counts (the usual case): the count itself is the index (count+1)
    idx_val = spikeTrain + 1;
else                                            % other values: index of each value within the sorted list of unique values
    [~, ~, idx_val] = unique(spikeTrain);
end
pk = accumarray(idx_val, 1) / N;                % empirical probability of each value
pk = pk(pk > 0);                                % only the observed values

% Shannon entropy in bits:
% H(X) = - sum_k p(k) log2 p(k)
Hbits = -sum(pk .* log2(pk));
if N == 0                                       % no bins: same result as before (zero)
    Hbits = 0;
end
Entropy = Hbits;

end