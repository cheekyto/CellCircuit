function cropped = cropped_datasets(folder, file_name)

% function to find the datasets written by "SelectElectrodes" among the recordings found in folders (November 2026), so that the recordings
% table does not list them as recordings of their own: "<recording>_<organoid_label>", in the folder of the recording "<recording>", whose
% file "<recording>_<organoid_label>_QCmanual.mat" holds spikeTimestamp_BeforeCropping (only "SelectElectrodes" writes it). The cropping
% is run, and run again (e.g., with other parameters), from the row of the recording itself (its file "<recording>_QCmanual.mat" is not
% changed by the cropping); the cropped dataset is written again at each run

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   folder:         folders of the datasets found (string array)
%   file_name:      their names, as in the column file_name of the recordings table (same size as folder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   cropped:        =true for the datasets written by "SelectElectrodes" (same size as file_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = string(folder);    file_name = string(file_name);
cropped = false(size(file_name));
for k = 1:numel(file_name)
    others = file_name(folder == folder(k) & file_name ~= file_name(k));
    if ~any(startsWith(file_name(k), others + "_"))                         % no recording "<recording>" in the same folder
        continue
    end
    file = fullfile(folder(k), file_name(k) + "_QCmanual.mat");
    if isfile(file)
        try
            cropped(k) = ismember('spikeTimestamp_BeforeCropping', who('-file', file));
        catch                                                               % (a file that cannot be read: kept as a recording)
        end
    else                                                                    % only the later files of the pipeline ("_curated", "_integrated")
        cropped(k) = ~isfile(fullfile(folder(k), file_name(k) + "_raw.mat")) && ~isfile(fullfile(folder(k), file_name(k) + "_cleaned.mat"));
    end
end

end
