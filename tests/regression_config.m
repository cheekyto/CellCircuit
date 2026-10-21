function cfg = regression_config()

% Configuration for the regression tests of the pipeline (used by "run_regression", "make_baseline" and "run_stage")
% Edit the values below to use a different reference recording or a different location for the temporary files
% (revised October 2026)

cfg.repo_dir = string(fileparts(fileparts(mfilename('fullpath'))));         % folder of the repository (parent of the folder "tests")
cfg.ref_folder = "/path/to/reference_recording";                           % reference recording, spike-sorted with Kilosort (its files are only read and copied, never modified)
cfg.base_name = "reference_recording_name";                                % common name of the data files of the reference recording
cfg.npy_folder = "/path/to/npy-matlab";                                     % npy-matlab library (readNPY / writeNPY), needed by "ShapeSpikeTrainData"
env_names = struct('ref_folder', 'CELLCIRCUIT_REF_FOLDER', 'base_name', 'CELLCIRCUIT_BASE_NAME', 'npy_folder', 'CELLCIRCUIT_NPY_FOLDER');
fields = fieldnames(env_names);                                             % the three values above can also be set with environment variables
for k=1:numel(fields)                                                       % (e.g., setenv('CELLCIRCUIT_REF_FOLDER', "...") before "run_regression")
    if ~isempty(getenv(env_names.(fields{k})))
        cfg.(fields{k}) = string(getenv(env_names.(fields{k})));
    end
end
cfg.work_root = string(fullfile(tempdir, 'CellCircuit_regression'));        % folder for the copies of the data, the versions of the code and the reports (~1 GB)
if ~isempty(getenv('CELLCIRCUIT_REGRESSION_DIR'))                           % the location can also be set with an environment variable
    cfg.work_root = string(getenv('CELLCIRCUIT_REGRESSION_DIR'));
end
cfg.matlab_exe = string(fullfile(matlabroot, 'bin', 'matlab'));             % MATLAB executable, to run each version of the code in a fresh MATLAB process
cfg.n_iteration_hs = 5;                                                     % number of shuffling iterations for "QualityControl_HyperSynch" (500 in the real analysis; reduced here to save time)
cfg.hs_n_clusters = 100;                                                    % number of clusters kept for "QualityControl_HyperSynch" and "QualityControl_Manual" (the first ones; =0, all)
                                                                            % with all the clusters (~800), the original "QualityControl_HyperSynch" would take weeks
cfg.analyze_conditions = 1;                                                 % recording conditions analyzed by "Analyze_spikeTrain" (=-1, all; the network burst detection takes ~30 min per condition)
cfg.tolerance = 1e-10;                                                      % relative tolerance to compare numeric values (to account for floating-point round-off)
cfg.max_report = 200;                                                       % maximum number of differences listed in the report (per output file)

end
