function file = package_CellCircuit(version, output_file)

% BUILDS THE INSTALL FILE OF CELLCIRCUIT FOR USERS WITHOUT CODING EXPERIENCE (new, November 2026)
% It packages the pipeline as a MATLAB toolbox (.mltbx). Users double-click the file (or drag it into MATLAB) to install it: the
% pipeline is added to their MATLAB path, and the window "CellCircuit" appears in the Apps tab (it can also be opened by typing
% CellCircuit). Installing a newer file replaces the older version.
%
%   package_CellCircuit("2.1.0")                         writes CellCircuit_2.1.0.mltbx in the folder of the repository
%   package_CellCircuit("2.1.0", "/path/CellCircuit.mltbx")
%
% Included: the functions of the main folder that are tracked by git (commit them first), the window (CellCircuit.m, common_parameters.m),
% the recordings template, README, LICENSE and the documentation. Not included: tests, legacy functions, data, and npy-matlab (each user
% sets its folder in the window, tab Settings). The version and the git commit are written in CellCircuit_version.txt, which "run_pipeline"
% saves in the history of each recording when the code is not a git repository. Tested with MATLAB R2025b.

if nargin < 1 || strlength(string(version)) == 0
    error("Give the version of the package, e.g., package_CellCircuit(""2.1.0"")");
end
version = string(version);
root = string(fileparts(fileparts(mfilename('fullpath'))));                 % the folder of the repository
if nargin < 2
    output_file = fullfile(root, "CellCircuit_" + version + ".mltbx");
end

% files: the .m files of the main folder tracked by git, plus the window
[status, out] = system(sprintf('git -C "%s" ls-files', root));
if status ~= 0
    error("git is needed to list the files of the package (only the committed files are packaged)");
end
tracked = splitlines(string(strtrim(out)));
main_m = string({dir(fullfile(root, '*.m')).name})';
included = main_m(ismember(main_m, tracked) | ismember(main_m, ["CellCircuit.m", "common_parameters.m"]));
left_out = setdiff(main_m, included);
if ~isempty(left_out)
    fprintf("Not packaged (not committed): %s\n", strjoin(left_out, ", "));
end
[~, commit] = system(sprintf('git -C "%s" rev-parse --short HEAD', root));
[~, dirty] = system(sprintf('git -C "%s" status --porcelain --untracked-files=no', root));
version_text = "CellCircuit " + version + " (commit " + strtrim(string(commit));
if ~isempty(strtrim(dirty))
    version_text = version_text + ", with uncommitted changes";
end
version_text = version_text + ")";
version_file = fullfile(root, "CellCircuit_version.txt");
fid = fopen(version_file, 'w');   fprintf(fid, '%s\n', version_text);   fclose(fid);
cleanup = onCleanup(@() delete(version_file));

files = [fullfile(root, included); version_file; fullfile(root, ["recordings_template.csv", "README.md", "LICENSE"])'; ...
         fullfile(root, "docs", "FUNCTION_REFERENCE.md")];
files = files(isfile(files));

opts = matlab.addons.toolbox.ToolboxOptions(char(root), "CellCircuit-MEA-pipeline");
opts.ToolboxName = "CellCircuit";
opts.ToolboxVersion = version;
opts.Summary = "Multielectrode-array (MEA) recordings of brain organoids: firing, bursting and network bursting, to Excel tables.";
opts.Description = "Run the whole pipeline from a window (Apps tab: CellCircuit, or type CellCircuit): choose the recordings, " + ...
    "the stages and the parameters, follow the progress, and open the Excel results. Needs the Signal Processing Toolbox, the Statistics " + ...
    "and Machine Learning Toolbox, and npy-matlab (to read Kilosort data).";
opts.AuthorName = "Sung Min Yang";
opts.ToolboxFiles = cellstr(files);
opts.ToolboxMatlabPath = {char(root)};
opts.AppGalleryFiles = {char(fullfile(root, "CellCircuit.m"))};
opts.OutputFile = char(output_file);
matlab.addons.toolbox.packageToolbox(opts);
file = string(output_file);
fprintf("Package written: %s (%s, %d files)\n", file, version_text, numel(files));

end
