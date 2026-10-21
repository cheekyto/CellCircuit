% =========================================================
% Helper function: similarity based on longest common
% subsequence length vs average string length
% =========================================================
function score = stringSimilarity(s1, s2)
    s1 = char(s1);
    s2 = char(s2);
    lcsLen = lcsLength(s1, s2);
    score  = (2 * lcsLen) / (length(s1) + length(s2));  % Sørensen-Dice coefficient
end
