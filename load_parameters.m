function P = load_parameters(source, changes)

% function to build the parameter set used by "run_pipeline": read a parameter file (or take a set saved before), apply a few changes,
% and copy the linked parameters (revised October 2026)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   source:     path of a parameter file (a script that builds the structure P, as "pipeline_parameters.m"),
%               or a structure (e.g., the set saved for a recording in "<file_name>_parameters.mat")
%   changes:    (optional) values changed for this run, as {'Group.NAME', value, 'Group.NAME', value, ...} (or one pair per row),
%               or a structure with the same fields as the set (e.g., struct('NetworkBurstDetection', struct('MIN_CHAN', 8)))
%               Only parameters that exist in the set can be changed (to catch typos); a linked parameter (see P.links) cannot be
%               changed directly: change the parameter it is linked to

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   P:          parameter set, with one field per function (P.<function>.<parameter>), the changes applied, and the linked
%               parameters copied

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    changes = {};
end

if isstruct(source)                                                         % a set saved before
    P = source;
else                                                                        % a parameter file
    file = char(source);
    if ~isfile(file)
        found = which(file);                                                % a file on the MATLAB path (e.g., "pipeline_parameters.m")
        if isempty(found)
            error("The parameter file '%s' was not found", file);
        end
        file = found;
    end
    P = read_parameter_file(file);
end
links = cell(0,2);
if isfield(P, 'links')
    links = P.links;
end

% apply the changes
list = changes_list(changes);
for k=1:size(list,1)                                                        % across all the changed parameters
    name = char(list{k,1});
    parts = strsplit(name, '.');
    if any(strcmp(links(:,1), name))
        source_name = links{strcmp(links(:,1), name), 2};
        error("The parameter '%s' is linked to '%s' (the same value): change '%s' instead", name, source_name, source_name);
    end
    if numel(parts) < 2 || ~has_field(P, parts)
        error("The parameter '%s' does not exist in the parameter set (check the name: 'Function.PARAMETER', or 'Function.method3.PARAMETER')", name);
    end
    P = setfield(P, parts{:}, list{k,2});
end

% copy the linked parameters
for k=1:size(links,1)
    source_parts = strsplit(links{k,2}, '.');
    target_parts = strsplit(links{k,1}, '.');
    if ~has_field(P, source_parts)
        error("The linked parameter '%s' does not exist in the parameter set", links{k,2});
    end
    P = setfield(P, target_parts{:}, getfield(P, source_parts{:}));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = read_parameter_file(file)

% function to run the parameter file (a script) and return the structure P that it builds

P = [];
run(file);                                                                  % the script defines P
if ~isstruct(P)
    error("The parameter file '%s' must define the structure P (see 'pipeline_parameters.m')", file);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function list = changes_list(changes)

% function to convert the changes to a list {name, value} (one row per changed parameter)

if isempty(changes)
    list = cell(0,2);
elseif isstruct(changes)                                                    % structure: one row per leaf of the structure
    list = flatten(changes, '');
elseif iscell(changes)
    if size(changes,2) == 2 && size(changes,1) > 1                          % one pair per row
        list = changes;
    else                                                                    % {'name', value, 'name', value, ...}
        if mod(numel(changes), 2) ~= 0
            error("The changes must be pairs {'Function.PARAMETER', value, ...}");
        end
        list = reshape(changes(:), 2, [])';
    end
else
    error("The changes must be a cell array {'Function.PARAMETER', value, ...} or a structure");
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function list = flatten(S, prefix)

% function to list the leaves of a structure, as {'a.b.c', value}

list = cell(0,2);
names = fieldnames(S);
for k=1:numel(names)
    name = names{k};
    if ~isempty(prefix)
        name = [prefix '.' names{k}];
    end
    value = S.(names{k});
    if isstruct(value)
        list = [list; flatten(value, name)];                                %#ok<AGROW>
    else
        list(end+1,:) = {name, value};                                      %#ok<AGROW>
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = has_field(S, parts)

% function to check that the nested field S.(parts{1}).(parts{2})... exists

ok = true;
for k=1:numel(parts)
    if ~isstruct(S) || ~isfield(S, parts{k})
        ok = false;
        return
    end
    S = S.(parts{k});
end

end
