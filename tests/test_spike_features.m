function test_spike_features(code_dirs, out_file, qcmanual_file, t_window)

% Function-level test for "SpikeFeatures" (and "Entropy", "Mutual_Information", "sttc", "lvr", "bin_spike_time")
% Run it once with the old and once with the new version of the code (each in a fresh MATLAB process), then compare both
% output files with "compare_outputs" (see "test_core_functions" for an example)
% (revised October 2026)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dirs:            cell array with the folders of the code to be tested (the first one has priority)
% out_file:             .mat file to save the results (and "<out_file>_timing.mat" with the running time)
% qcmanual_file:        file with the variable "spikeTimestamp" (e.g., the "_curated_QCmanual.mat" file used by "Analyze_spikeTrain")
% t_window:             analysis window, in sec (default: [300 900])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4, t_window = [300 900]; end
for k=numel(code_dirs):-1:1                                                 % added in reverse order, so that the first folder has priority
    addpath(char(code_dirs{k}));
end
work_dir = tempname;    mkdir(work_dir);    cd(work_dir);                   % neutral current folder

load(qcmanual_file, 'spikeTimestamp');
spikeTrain = Reshape_spikeTimestamp(spikeTimestamp, t_window);              % as in "Analyze_spikeTrain"
t0 = tic;
SpikeFeature = SpikeFeatures(spikeTrain, t_window);
elapsed = toc(t0);

save(out_file, 'SpikeFeature');
[p, f] = fileparts(out_file);
save(fullfile(p, [f '_timing.mat']), 'elapsed');
fprintf('SpikeFeatures: %.2f s for %g clusters\n', elapsed, numel(spikeTrain));

end
