function override_inputs(inputs)

% function to overwrite the INPUTS / PARAMETERS defined inside the code of a pipeline function (e.g., "QualityControl_HyperSynch")
% with the values provided in the structure "inputs" (each field must have the same name than the variable to be replaced)
% It allows to run the pipeline from a script (e.g., the test harness in the folder "tests") without editing the code
% It must be called by the pipeline function itself, right after its sections INPUTS and PARAMETERS
% (revised October 2026)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs.(name):        new value for the variable "name" in the workspace of the calling function
%                       (the field "system" is skipped, because it is used earlier to set up the MATLAB path)

% NOTE: variables that are derived from an input inside the INPUTS section (e.g., "ellipse_cropped" in "QualityControl_General")
%       are NOT recalculated, so they need to be provided as well

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Replace the variables in the workspace of the calling function

names = fieldnames(inputs);                                                 % list of variables to be replaced
for k=1:numel(names)                                                        % across all variables to be replaced
    if strcmp(names{k}, 'system')                                           % "system" was already used by the calling function to set up the MATLAB path
        continue
    end
    if ~evalin('caller', sprintf('exist(''%s'',''var'')', names{k}))        % the variable must already exist in the calling function (to catch typos in the field names)
        stack = dbstack(1);                                                 % name of the calling function, for the error message
        error("The input '%s' is not defined in the function '%s'", names{k}, stack(1).name);
    end
    assignin('caller', names{k}, inputs.(names{k}));                        % replace the value of the variable in the calling function
end

end
