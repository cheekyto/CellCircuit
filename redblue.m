%% =========================================================================
%  LOCAL FUNCTION – diverging red-white-blue colormap
% =========================================================================
function c = redblue(n)
    % Red-white-blue diverging colormap (n levels)
    half = floor(n/2);
    r1 = [linspace(0.7,1,half)', linspace(0.1,1,half)', linspace(0.1,1,half)'];
    r2 = [linspace(1,0.1,n-half)', linspace(1,0.1,n-half)', linspace(1,0.7,n-half)'];
    c  = [r1; r2];
end
