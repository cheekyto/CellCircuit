function [base_cond, by_name] = baselineCondition(file_names)

% (new, November 2026; user's choice) index of the baseline among the recording conditions of a recording, identified by the file name of
% each condition instead of its position. The baseline is the condition whose file name (without the extension ".brw" / ".bxr") is the
% start of the file name of every other condition, e.g. "rec_D3.brw" is the baseline of "rec_D3_ttx_00.brw" and "rec_D3_01.brw".
% Before, "Analyze_spikeTrain" always took the first condition (the first row of recording_data.csv) as the baseline: when the TTX file was
% listed first, the TTX condition was saved as "Baseline" and the baseline as "unknown".
% When no condition fits this rule (e.g., the 3Brain names "<name>_Baseline" and "<name>_TTX"), the first condition is the baseline, as before.
% If two conditions have the same name without extension, the first of them is taken.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   file_names:         file name of each recording condition (cell or string array; e.g., rec_conditions(:,1) in "Analyze_spikeTrain")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   base_cond:          index of the baseline condition
%   by_name:            =true, if the baseline was identified by its file name (or there is a single condition)
%                       =false, if no condition fits the rule and the first condition was taken, as before

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = string(file_names(:));
n = numel(names);
base_cond = 1;
by_name = n == 1;
if n < 2
    return
end

stems = strings(n,1);
for k=1:n
    aux = split(names(k), ".b");                                            % name without the extension (as in "Analyze_spikeTrain", which assumes
    stems(k) = aux(1);                                                      % that the raw data file was a .brw or .bxr file)
end

for k=1:n                                                                   % the first condition (in the order of the list) that fits the rule
    others = names([1:k-1, k+1:n]);
    if strlength(stems(k)) > 0 && all(startsWith(others, stems(k)))
        base_cond = k;
        by_name = true;
        return
    end
end

end
