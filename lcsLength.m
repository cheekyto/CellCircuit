% =========================================================
% Helper function: longest common subsequence (LCS) length
% =========================================================
function len = lcsLength(s1, s2)
    m = length(s1);
    n = length(s2);
    dp = zeros(m+1, n+1);
    for i = 1:m
        for j = 1:n
            if s1(i) == s2(j)
                dp(i+1, j+1) = dp(i, j) + 1;
            else
                dp(i+1, j+1) = max(dp(i, j+1), dp(i+1, j));
            end
        end
    end
    len = dp(m+1, n+1);
end
