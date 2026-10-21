function apply_parameters(group, subgroup)

% function to overwrite the PARAMETERS defined inside the code of a pipeline function with the values of the active parameter set
% (revised October 2026). It must be called by the pipeline function itself, at the end of its PARAMETERS section
% The set is activated by "run_pipeline" (see "parameter_set"); without an active set nothing changes, and the values written in the
% code are used (e.g., a stage run by hand, or the tests)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   group:      name of the calling function (field of the parameter set, e.g., "NetworkBurstDetection")
%   subgroup:   (optional) field of the group with the parameters that depend on an option (e.g., "method3")

% NOTE: every parameter of the set must already exist in the calling function (to catch typos in "pipeline_parameters"); it is checked
%       the first time a group is used after a set is activated. The fields that are structures are subgroups, and they are skipped

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent cache_id cache

[P, id] = parameter_set('get');
if isempty(P)                                                               % no active set: the values written in the code are used
    return
end
if ~isequal(cache_id, id)                                                   % a new set was activated (or cleared): empty the cache
    cache_id = id;
    cache = struct();
end
if nargin > 1
    key = [group '__' subgroup];
else
    key = group;
end

if ~isfield(cache, key)                                                     % first use of this group with the active set
    names = {};     vals = {};
    if isfield(P, group)
        values = P.(group);
        if nargin > 1
            if isfield(values, subgroup)
                values = values.(subgroup);
            else
                values = struct();
            end
        end
        names = fieldnames(values);
        is_value = ~cellfun(@(n) isstruct(values.(n)), names);             % the structures are subgroups
        names = names(is_value);
        vals = cellfun(@(n) values.(n), names, 'UniformOutput', false);
        for k=1:numel(names)                                                % the parameter must exist in the calling function (typos)
            if ~evalin('caller', sprintf('exist(''%s'',''var'')', names{k}))
                error("The parameter '%s' of the parameter set (%s) is not defined in the function '%s'", names{k}, strrep(key, '__', '.'), group);
            end
        end
    end
    cache.(key) = {names, vals};
end

entry = cache.(key);
names = entry{1};   vals = entry{2};
for k=1:numel(names)                                                        % replace the values in the calling function
    assignin('caller', names{k}, vals{k});
end

end
