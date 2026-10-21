%% =========================================================================
%  LOCAL FUNCTION – safe linkage with automatic fallback
% =========================================================================
function [Z, method_used] = safe_linkage(Y, preferred_method)
% Attempts linkage with preferred_method (e.g. 'ward').
% If MATLAB raises a non-monotonic-tree warning (which means the Euclidean
% MDS embedding did not fully resolve the composite distance geometry),
% falls back automatically to 'average' linkage, which is valid for any
% distance matrix and does not require a monotonic tree.
%
% Returns Z (linkage matrix) and method_used (string) so callers can
% report which method was actually applied.

    FALLBACK = 'average';
    apply_parameters('safe_linkage');                                       % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

    % Temporarily promote the non-monotonic warning to an error so we can
    % catch it cleanly without it printing to the console
    warnID = 'stats:linkage:NonMonotonicTree';
    prev   = warning('error', warnID);

    try
        Z           = linkage(Y, preferred_method);
        method_used = preferred_method;
    catch
        warning(prev);   % restore warning state before the fallback call
        fprintf(['Warning: Ward linkage produced a non-monotonic tree ' ...
                 '(MDS embedding insufficient).\n' ...
                 'Falling back to ''%s'' linkage automatically.\n\n'], FALLBACK);
        Z           = linkage(Y, FALLBACK);
        method_used = FALLBACK;
    end

    warning(prev);   % always restore warning state on clean exit too
end
