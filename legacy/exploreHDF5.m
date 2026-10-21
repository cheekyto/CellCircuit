function exploreHDF5(filename)
% EXPLOREHDF5 Recursively explores all groups and attributes in an HDF5 file
%   exploreHDF5('yourfile.h5') prints the full structure with attribute values

    fprintf('\n========================================\n');
    fprintf('HDF5 File Structure: %s\n', filename);
    fprintf('========================================\n\n');

    info = h5info(filename);
    
    % --- Root-level attributes ---
    fprintf('[ROOT: /]\n');
    printAttributes(filename, '/', info.Attributes);
    
    % --- Recursively explore all groups ---
    for i = 1:length(info.Groups)
        exploreGroup(filename, info.Groups(i), 1);
    end

    fprintf('\n========================================\n');
    fprintf('Exploration complete.\n');
    fprintf('========================================\n');
end


function exploreGroup(filename, group, depth)
% Recursively prints group info, attributes, and datasets

    indent = repmat('  ', 1, depth);  % indentation for readability
    
    fprintf('\n%s[GROUP: %s]\n', indent, group.Name);
    
    % --- Attributes of this group ---
    printAttributes(filename, group.Name, group.Attributes);
    
    % --- Datasets in this group ---
    if ~isempty(group.Datasets)
        fprintf('%s  Datasets:\n', indent);
        for i = 1:length(group.Datasets)
            ds = group.Datasets(i);
            dsPath = [group.Name, '/', ds.Name];
            fprintf('%s    - %s  [Size: %s | Type: %s]\n', ...
                indent, ds.Name, ...
                mat2str(ds.Dataspace.Size), ...
                ds.Datatype.Class);
            
            % Also print dataset-level attributes if any
            if ~isempty(ds.Attributes)
                fprintf('%s      Dataset attributes:\n', indent);
                printAttributes(filename, dsPath, ds.Attributes);
            end
        end
    end
    
    % --- Recurse into subgroups ---
    for i = 1:length(group.Groups)
        exploreGroup(filename, group.Groups(i), depth + 1);
    end
end


function printAttributes(filename, location, attributes)
% Prints attribute names and their values for a given HDF5 location

    if isempty(attributes)
        return;
    end
    
    indent = '      ';  % fixed indent for attribute lines
    fprintf('%sAttributes:\n', indent);
    
    for i = 1:length(attributes)
        attrName = attributes(i).Name;
        
        try
            val = h5readatt(filename, location, attrName);
            
            % Format value based on type for clean printing
            if ischar(val) || isstring(val)
                fprintf('%s  %-30s = "%s"\n', indent, attrName, val);
            elseif isnumeric(val) && isscalar(val)
                fprintf('%s  %-30s = %g\n', indent, attrName, val);
            elseif isnumeric(val)
                fprintf('%s  %-30s = [%s]\n', indent, attrName, mat2str(val));
            elseif iscell(val)
                fprintf('%s  %-30s = {cell: %s}\n', indent, attrName, mat2str(size(val)));
            else
                fprintf('%s  %-30s = [unrecognized type: %s]\n', indent, attrName, class(val));
            end
            
        catch ME
            fprintf('%s  %-30s = [ERROR reading: %s]\n', indent, attrName, ME.message);
        end
    end
end