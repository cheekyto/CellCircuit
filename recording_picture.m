function [picture, candidates] = recording_picture(folder, file_name)

% function to find the photo of the chip of a recording, for the column "picture" of the recordings table (November 2026): the image file
% in the folder of the recording (.jpg, .jpeg, .heic, .heif, .png, .tif, .tiff; not in its subfolders), without the figures written by
% the pipeline (names that start with the name of the recording, or with "_crop_" in them) and the hidden files. The photos (.jpg,
% .jpeg, .heic, .heif) are preferred to the .png and .tif images. When there are several, none is chosen (the user writes it)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   folder:         folder of the recording
%   file_name:      (optional) common name of the files of the recording (its figures are not photos of the chip)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   picture:        file name of the photo ("" when there is none, or several)
%   candidates:     all the image files found (string array)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2, file_name = ""; end
picture = "";   candidates = strings(0, 1);
if ~isfolder(folder)
    return
end
f = dir(folder);
f = f(~[f.isdir]);
names = string({f.name})';
if isempty(names)
    return
end
ext = lower(regexp(names, '\.[^.]+$', 'match', 'once'));
ext(ismissing(ext)) = "";
photo_ext = [".jpg", ".jpeg", ".heic", ".heif"];
keep = ismember(ext, [photo_ext, ".png", ".tif", ".tiff"]) & ~startsWith(names, ".") & ~contains(names, "_crop_");
if strlength(string(file_name)) > 0
    keep = keep & ~startsWith(names, string(file_name));
end
candidates = names(keep);
photos = candidates(ismember(ext(keep), photo_ext));
if isscalar(photos)
    picture = photos;
elseif isempty(photos) && isscalar(candidates)
    picture = candidates;
end

end
