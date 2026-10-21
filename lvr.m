function local_variation_ratio = lvr(data)

% function to calculate the local variation (Shinomoto et al., 2003) of a sequence of inter-spike intervals:
%   LV = 3/(n-1) * sum_i [ (I(i)-I(i+1)) / (I(i)+I(i+1)) ]^2
% (revised October 2026: vectorized, without the loop across the intervals; the result only differs by floating-point round-off)

%   data:                   inter-spike intervals (vector)
%   local_variation_ratio:  local variation (NaN when there are less than MIN_N intervals)

MIN_N = 3;
apply_parameters('lvr');                                                    % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

nPoints = numel(data);
if nPoints < MIN_N
    local_variation_ratio = nan;
else
    data = data(:);                                                         % column vector
    aux = 3*(data(1:end-1)-data(2:end)).^2 ./ (data(1:end-1)+data(2:end)).^2;   % one term per pair of consecutive intervals (same expression as in the original loop)
    local_variation_ratio = sum(aux) / (nPoints-1);
end

end
