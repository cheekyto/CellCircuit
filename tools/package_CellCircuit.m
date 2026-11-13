function file = package_CellCircuit(version, output_file)

% BUILDS THE INSTALL FILE OF CELLCIRCUIT FOR USERS WITHOUT CODING EXPERIENCE (new, November 2026)
% It packages the pipeline as a MATLAB toolbox (.mltbx). Users double-click the file (or drag it into MATLAB) to install it: the
% pipeline is added to their MATLAB path, and the window "CellCircuit" appears in the Apps tab (it can also be opened by typing
% CellCircuit). Installing a newer file replaces the older version.
%
%   package_CellCircuit("2.1.0")                         writes CellCircuit_2.1.0.mltbx in the folder of the repository
%   package_CellCircuit("2.1.0", "/path/CellCircuit.mltbx")
%
% (revised, November 2026) Only the committed code is packaged: the files come from a clean copy of the last commit (git archive), so
% uncommitted changes and files of the working folder never reach the package (before, the files were read from the working folder, so
% uncommitted edits of committed files, e.g., of the documentation, were packaged too); the uncommitted changes are listed.
% Included: the functions of the main folder (with the window, CellCircuit.m and common_parameters.m), the recordings template, README,
% LICENSE and the documentation. Not included: tests, tools, legacy functions, data, and npy-matlab (each user sets its folder in the
% window, tab Settings). The version and the commit are written in CellCircuit_version.txt, which "run_pipeline" saves in the history of
% each recording when the code is not a git repository. Tested with MATLAB R2025b.

if nargin < 1 || strlength(string(version)) == 0
    error("Give the version of the package, e.g., package_CellCircuit(""2.1.0"")");
end
version = string(version);
root = string(fileparts(fileparts(mfilename('fullpath'))));                 % the folder of the repository
if nargin < 2
    output_file = fullfile(root, "CellCircuit_" + version + ".mltbx");
end

% clean copy of the last commit
[status, commit] = system(sprintf('git -C "%s" rev-parse --short HEAD', root));
if status ~= 0
    error("git is needed: the package is built from the last commit of the repository");
end
commit = strtrim(string(commit));
copy = string(tempname);
zip_file = copy + ".zip";
cleanup = onCleanup(@() remove_copy(copy, zip_file));
[status, out] = system(sprintf('git -C "%s" archive --format=zip -o "%s" HEAD', root, zip_file));
if status ~= 0
    error("The copy of the last commit could not be made: %s", out);
end
unzip(zip_file, copy);

% what is left out: the uncommitted changes of the working folder
[~, dirty] = system(sprintf('git -C "%s" status --porcelain', root));
if strlength(strtrim(string(dirty))) > 0
    fprintf("Not packaged (uncommitted changes of the working folder; the package has the commit %s):\n%s", commit, dirty);
end

version_text = "CellCircuit " + version + " (commit " + commit + ")";
version_file = fullfile(copy, "CellCircuit_version.txt");
fid = fopen(version_file, 'w');   fprintf(fid, '%s\n', version_text);   fclose(fid);

main_m = string({dir(fullfile(copy, '*.m')).name})';
files = [fullfile(copy, main_m); version_file; fullfile(copy, ["recordings_template.csv", "README.md", "LICENSE"])'; ...
         fullfile(copy, "docs", "FUNCTION_REFERENCE.md"); fullfile(copy, "OrganoidsFromPhoto_chip.mat")];   % (pattern of the chip, November 2026)
files = files(isfile(files));
if ~any(main_m == "CellCircuit.m")
    error("The last commit has no CellCircuit.m");
end

opts = matlab.addons.toolbox.ToolboxOptions(char(copy), "CellCircuit-MEA-pipeline");
opts.ToolboxName = "CellCircuit";
opts.ToolboxVersion = version;
opts.Summary = "Multielectrode-array (MEA) recordings of brain organoids: firing, bursting and network bursting, to Excel tables.";
opts.Description = "Run the whole pipeline from a window (Apps tab: CellCircuit, or type CellCircuit): choose the recordings, " + ...
    "the stages and the parameters, follow the progress, and open the Excel results. Needs the Signal Processing Toolbox, the Statistics " + ...
    "and Machine Learning Toolbox, and npy-matlab (to read Kilosort data).";
opts.AuthorName = "Sung Min Yang";
opts.ToolboxFiles = cellstr(files);
opts.ToolboxMatlabPath = {char(copy)};
opts.AppGalleryFiles = {char(fullfile(copy, "CellCircuit.m"))};
opts.OutputFile = char(output_file);
matlab.addons.toolbox.packageToolbox(opts);
file = string(output_file);
fprintf("Package written: %s (%s, %d files)\n", file, version_text, numel(files));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function remove_copy(copy, zip_file)
if isfolder(copy)
    rmdir(copy, 's');
end
if isfile(zip_file)
    delete(zip_file);
end
end
