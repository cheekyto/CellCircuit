function [out, out_id] = parameter_set(action, P)

% function to keep the parameter set that is active while "run_pipeline" runs (revised October 2026)
% While a set is active, every pipeline function replaces the values written in its PARAMETERS section by the values of the set
% (see "apply_parameters"). Without an active set (e.g., a stage run by hand), the values written in the code are used, as before

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   action:     'set', activate the set P (structure with one field per function, as built by "pipeline_parameters")
%               'get', return the active set ([] if none), and its id
%               'clear', deactivate the set
%               'active', return true if a set is active
%               'id', return a number that changes each time a set is activated (used by "apply_parameters")
%   P:          parameter set, only for 'set'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent current id
if isempty(id)
    id = 0;
end

out = [];
out_id = id;
switch action
    case 'set'
        current = P;
        id = id + 1;
    case 'get'
        out = current;
    case 'clear'
        current = [];
    case 'active'
        out = ~isempty(current);
    case 'id'
        out = id;
    otherwise
        error("Unknown action '%s' (use 'set', 'get', 'clear', 'active' or 'id')", action);
end

end
