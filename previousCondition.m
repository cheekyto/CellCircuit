function prev = previousCondition(windows, base_cond, excluded)

% (new, November 2026; user's choice) the condition recorded just before each condition, for the comparison of their coordinated
% neuronal ensembles ("Analyze_cNE", "Compile_MEAproject"; e.g., with a baseline and three drugs, drug 3 is compared with the baseline and
% with drug 2): the condition whose analysis window starts last before the start of the window of the condition (the order of the
% recording, whatever the order of the list of conditions). 0 = no such comparison: the condition recorded before is the baseline (the
% comparison with the baseline already covers it), the condition is the baseline or the first one, it has no window, or it is excluded

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   windows:    analysis window [start end] of each condition, in sec: a cell array (e.g., {PerCondition.window}) or a matrix with one row
%               per condition; [] or NaN = no window (the condition is left out)
%   base_cond:  index of the baseline condition
%   excluded:   (optional) logical vector, true = the condition is left out (it has no previous condition and is no one's previous
%               condition; e.g., the TTX conditions, left out of the cNE analysis)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   prev:       row vector, prev(c) = index of the condition recorded just before the condition c (0 = none, see above)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(windows)
    windows = num2cell(windows, 2);
end
n = numel(windows);
start = nan(1, n);
if nargin < 3 || isempty(excluded)
    excluded = false(1, n);
end
for c=1:n
    w = windows{c};
    if numel(w) >= 1 && isfinite(w(1)) && ~excluded(c)
        start(c) = w(1);
    end
end
prev = zeros(1, n);
for c=1:n
    if c == base_cond || isnan(start(c))
        continue
    end
    before = find(start < start(c));
    if isempty(before)
        continue
    end
    [~, k] = max(start(before));
    if before(k) ~= base_cond
        prev(c) = before(k);
    end
end
end
