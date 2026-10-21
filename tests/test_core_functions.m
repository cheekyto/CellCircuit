function test_core_functions(code_dirs, out_file, spike_file)

% Function-level test for the core functions shared by all the stages (BurstDetection, BurstFeatures, lvr, Reshape_spikeTimestamp)
% Run it once with the old and once with the new version of the code (each in a fresh MATLAB process), then compare both
% output files with "compare_outputs"
% (revised October 2026)
%
% Example (from a terminal):
%   matlab -batch "addpath('<repo>/tests'); test_core_functions({'<old code folder>'}, 'core_old.mat', '<..._raw.mat>')"
%   matlab -batch "addpath('<repo>/tests'); test_core_functions({'<repo>'}, 'core_new.mat', '<..._raw.mat>')"
%   matlab -batch "addpath('<repo>/tests'); [same, report] = compare_outputs('core_old.mat', 'core_new.mat', 1e-10)"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dirs:            cell array with the folders of the code to be tested (the first one has priority)
% out_file:             .mat file to save the results (and "timing_file" with the running times, next to it)
% spike_file:           .mat file with a structure array "spikeTimestamp" (e.g., the output of "ShapeSpikeTrainData")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up the MATLAB path (only the code under test)

for k=numel(code_dirs):-1:1                                                 % added in reverse order, so that the first folder has priority
    addpath(char(code_dirs{k}));
end
work_dir = tempname;    mkdir(work_dir);    cd(work_dir);                   % neutral current folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spike trains to be tested: all clusters of a real recording + synthetic edge cases

load(spike_file, 'spikeTimestamp');
n_real = numel(spikeTimestamp);
trains = cell(n_real, 1);
for i=1:n_real
    trains{i} = spikeTimestamp(i).timeStamp;                                % real spike trains, in msec
end

rng(0);                                                                     % synthetic spike trains (reproducible)
synthetic = {
    zeros(0,1);                                                             % empty spike train
    [];                                                                     % empty spike train (0x0, as saved by "ShapeSpikeTrainData")
    1000;                                                                   % 1 spike
    [1000; 1100];                                                           % 2 spikes
    [1000; 1100; 1150];                                                     % 3 spikes (minimum burst)
    (0:100:10000)';                                                         % regular 10 Hz (no bursts: isi = 100 < 200, but no peak isi < 150)
    (0:50:10000)';                                                          % regular 20 Hz (one long burst)
    cumsum([0; 200; 150; 300; 149; 315; 300])+1000;                         % isi exactly at the thresholds
    cumsum([0; 100*ones(30,1); 1000; 20*ones(5,1); 600; 20*ones(5,1)]);     % long burst (maturation of the threshold) + bursts close to each other
    cumsum(exprnd(300, 2000, 1));                                           % Poisson spike train (~3 Hz)
    sort([cumsum(exprnd(2000, 200, 1)); reshape(cumsum(exprnd(2000, 200, 1)) + cumsum(exprnd(15, 200, 6), 2), [], 1)])};   % bursty spike train
trains = [trains; synthetic];
n_trains = numel(trains);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the functions

BD = cell(n_trains,1);  BF = cell(n_trains,1);  LV = nan(n_trains,1);  err = strings(n_trains,1);
t_bd = 0;   t_bf = 0;   t_lv = 0;
for i=1:n_trains
    try
        t0 = tic;   BD{i} = BurstDetection(trains{i}, 1);   t_bd = t_bd + toc(t0);
        t0 = tic;   BF{i} = BurstFeatures(BD{i});           t_bf = t_bf + toc(t0);
    catch ME
        err(i) = string(ME.message);                                        % the old version crashes with an empty spike train
    end
    t0 = tic;   LV(i) = lvr(diff(trains{i}));   t_lv = t_lv + toc(t0);
end

windows = [0 300; 300 900; 1e5 1e5+10];                                     % in sec (the last one has no spikes)
R = cell(size(windows,1),1);
t0 = tic;
for w=1:size(windows,1)
    R{w} = Reshape_spikeTimestamp(spikeTimestamp, windows(w,:));
end
t_reshape = toc(t0);

save(out_file, 'BD', 'BF', 'LV', 'err', 'R');
timing = struct('BurstDetection', t_bd, 'BurstFeatures', t_bf, 'lvr', t_lv, 'Reshape_spikeTimestamp', t_reshape);
[p, f] = fileparts(out_file);
save(fullfile(p, [f '_timing.mat']), 'timing');
fprintf('BurstDetection %.2f s | BurstFeatures %.2f s | lvr %.3f s | Reshape_spikeTimestamp %.3f s\n', t_bd, t_bf, t_lv, t_reshape);
fprintf('Crashes: %g (of %g spike trains)\n', sum(strlength(err) > 0), n_trains);

end
