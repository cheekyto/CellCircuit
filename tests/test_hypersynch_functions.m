function test_hypersynch_functions(code_dirs, out_file, cleaned_file, n_pairs, n_shuffle)

% Function-level test for "HyperSyncPair" and "ShuffleSpikeTrain" (methods 1, 2 and 3)
% Run it once with the old and once with the new version of the code (each in a fresh MATLAB process), then compare both
% output files with "compare_outputs" (see "test_core_functions" for an example)
% (revised October 2026)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dirs:            cell array with the folders of the code to be tested (the first one has priority)
% out_file:             .mat file to save the results (and "<out_file>_timing.mat" with the running times)
% cleaned_file:         output of "QualityControl_General" (variables "spike_plot" and "discarded_cluster")
% n_pairs:              number of random pairs of clusters to test "HyperSyncPair" (default: 2000)
% n_shuffle:            number of shuffling iterations to test "ShuffleSpikeTrain" (default: 5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4, n_pairs = 2000; end
if nargin < 5, n_shuffle = 5; end
for k=numel(code_dirs):-1:1                                                 % added in reverse order, so that the first folder has priority
    addpath(char(code_dirs{k}));
end
work_dir = tempname;    mkdir(work_dir);    cd(work_dir);                   % neutral current folder

load(cleaned_file, 'spike_plot', 'discarded_cluster');
st = spike_plot(discarded_cluster == 0);                                    % good clusters (as in "QualityControl_HyperSynch")
n_ch = numel(st);
t_step = st(1).t_step;
range = st(1).window/1e-3;                                                  % in msec
BD = cell(n_ch,1);  BF = cell(n_ch,1);  sig = cell(n_ch,1);
for i=1:n_ch
    BD{i} = BurstDetection(st(i).timeStamp, 1);
    BF{i} = BurstFeatures(BD{i});
    sig{i} = [BD{i}.spikes(:,1) BD{i}.spikes(:,3)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HyperSyncPair on random pairs, on pairs of a cluster with itself (every spike is hypersynchronized) and on edge cases

rng(1);
pairs = [randi(n_ch, n_pairs, 1), randi(n_ch, n_pairs, 1)];
pairs(1:10,2) = pairs(1:10,1);                                              % a cluster against itself
HS = cell(n_pairs + 4, 1);
t0 = tic;
for p=1:n_pairs
    HS{p} = HyperSyncPair(sig{pairs(p,1)}, sig{pairs(p,2)}, t_step);
end
t_pair = toc(t0);
s = sig{pairs(1,1)};
HS{n_pairs+1} = HyperSyncPair(s, zeros(0,2), t_step);                       % empty channel 2
HS{n_pairs+2} = HyperSyncPair(zeros(0,2), s, t_step);                       % empty channel 1
HS{n_pairs+3} = HyperSyncPair(s, flipud(s), t_step);                        % channel 2 not sorted
HS{n_pairs+4} = HyperSyncPair(s, [s(:,1)+0.5*t_step s(:,2); s(:,1)-0.5*t_step s(:,2)], t_step);  % several spikes of channel 2 close to each spike of channel 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ShuffleSpikeTrain (methods 1, 2 and 3) for a few clusters

nb = cellfun(@(x) x.nBurst, BF);
[~, order] = sort(nb, 'descend');
test_ch = order([1 2 round(n_ch/4) round(n_ch/2) n_ch]);                    % clusters with many, some and no bursts
SH = cell(numel(test_ch), 3);   t_shuffle = zeros(1,3);
for m=1:3
    t0 = tic;
    for c=1:numel(test_ch)
        i = test_ch(c);
        SH{c,m} = ShuffleSpikeTrain(sig{i}(:,1), range, t_step, m, n_shuffle, 0, BD{i}, BF{i});
    end
    t_shuffle(m) = toc(t0);
end

save(out_file, 'HS', 'SH', 'pairs', 'test_ch');
timing = struct('HyperSyncPair', t_pair, 'ShuffleSpikeTrain', t_shuffle);
[p, f] = fileparts(out_file);
save(fullfile(p, [f '_timing.mat']), 'timing');
fprintf('HyperSyncPair: %.2f s for %g pairs | ShuffleSpikeTrain (methods 1/2/3): %.2f / %.2f / %.2f s\n', t_pair, n_pairs, t_shuffle);

end
