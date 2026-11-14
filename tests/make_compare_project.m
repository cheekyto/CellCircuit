function info = make_compare_project(project)

% synthetic project for the tests of "CompareNeurons" and of tab 4 of "CellCircuit" (new, November 2026): three groups of 8 organoids
% (A, B, C), with a recordings table (folder, file_name, group, batch, include), the group folders with the "_Analyzed.mat" files
% (PerCondition with TTX listed before the baseline, spikeTimestamp, NeuronScore) and the compiled tables (Firing, NETbursting1).
% Known effects: the units of B are broader (higher waveform score 1), the rate of the network bursts follows the mean of that score of
% each organoid, and the broader units start the network bursts earlier. C_org8 has no saved scores and no templates; A_org1, A_org2 and
% B_org1 have no network bursts. info: names, baseline rate and firing rate of each organoid, the organoid without scores, the number of
% bursting organoids of A and B
mkdir(project);
R = NeuronScores(struct('relativeRate', 0, 'lv', 1, 'cv', 1, 'burstIndex', 0), [], [], [], NaN).reference;
D.a1 = R.WF_AXES(:,1) / norm(R.WF_AXES(:,1));
a2 = R.WF_AXES(:,2) - (D.a1' * R.WF_AXES(:,2)) * D.a1;
D.a2 = a2 / norm(a2);
D.b = R.FP_AXIS / norm(R.FP_AXIS);
groups = ["A", "B", "C"];
shift = [0 1 0];                                                            % B: broader units
not_bursting = ["A_org1", "A_org2", "B_org1"];
rows = strings(0, 5);
info.names = strings(0, 1);  info.rate = [];  info.firing = [];
info.n_bursting = 0;  info.n_role_organoids = 0;
for g=1:3
    mkdir(fullfile(project, groups(g)));
    Firing = table();  NETbursting1 = table();
    for k=1:8
        name = groups(g) + "_org" + k;
        rec = name + "_curated";
        n = 60 + randi(40);
        z1 = shift(g) + 2 * randn(n, 1);
        [NS, ~] = synth(R, D, z1);
        nb = 20 * ~ismember(name, not_bursting);
        P = conditions(name, n, nb, z1);
        spikeTimestamp = struct('cluster', num2cell(1:n));                  %#ok<NASGU>
        file = fullfile(project, groups(g), rec + "_Analyzed.mat");
        if name == "C_org8"                                                 % no saved scores, no templates
            P(2).spikeTrain = arrayfun(@(u) struct('timeStamp', sort(rand(1, 300) * 600e3)), 1:n);
            P(2).window = [0 600];
            PerCondition = P;                                               %#ok<NASGU>
            save(file, 'PerCondition', 'spikeTimestamp');
            info.no_scores_file = file;  info.no_scores_name = rec;
        else
            NS.cluster = (1:n)';  NS.condition = 'Baseline';
            NeuronScore = NS;  PerCondition = P;                            %#ok<NASGU>
            save(file, 'PerCondition', 'spikeTimestamp', 'NeuronScore');
        end
        rate = NaN;  part = NaN;  dur = NaN;
        if nb > 0
            rate = 10^(0.3 + 0.5 * mean(z1) + 0.03 * randn);                % follows the composition
            part = 0.5 + 0.1 * rand;  dur = 0.5 + rand;
        end
        firing = 10^(0.2 * randn);
        if groups(g) ~= "C" && nb > 0
            info.n_bursting = info.n_bursting + 1;
            info.n_role_organoids = info.n_role_organoids + 1;
        end
        info.names(end+1) = rec;  info.rate(end+1) = rate;  info.firing(end+1) = firing;
        Firing = [Firing; table({char(rec); char(rec)}, {char(name + "_ttx_00.brw"); char(name + ".brw")}, {'+TTX'; 'Baseline'}, ...
            [0.01; firing], 'VariableNames', {'Organoid', 'File', 'Condition', 'Firing_Rate'})];   %#ok<AGROW>
        NETbursting1 = [NETbursting1; table({char(rec); char(rec)}, {char(name + "_ttx_00.brw"); char(name + ".brw")}, ...
            {'+TTX'; 'Baseline'}, [0; nb], [NaN; rate], [NaN; part], [NaN; dur], 'VariableNames', {'Organoid', 'File', 'Condition', ...
            'Number_NETburst_Events', 'Rate_NETburst', 'Participation_NETburst', 'Duration_NETburst'})];   %#ok<AGROW>
        rows(end+1, :) = [fullfile(project, "raw", name), name, groups(g), "b" + mod(k, 3), "1"];   %#ok<AGROW>
    end
    save(fullfile(project, groups(g) + "_Compiled.mat"), 'Firing', 'NETbursting1');
end
writetable(array2table(rows, 'VariableNames', {'folder', 'file_name', 'group', 'batch', 'include'}), fullfile(project, 'recordings.csv'));
end

function P = conditions(name, n, nb, z1)
% TTX (listed first) and baseline; primary network bursts of the baseline: 70% of the units take part, the broader ones start earlier
P = struct('file', {char(name + "_ttx_00.brw"), char(name + ".brw")}, 'condition', {'+TTX', 'Baseline'}, ...
    'NetworkBurstDetected', {struct('nBurst', 0), []}, 'NetworkBurstStructure', {[], []}, 'spikeTrain', {[], []}, 'window', {[], []});
N.nBurst = nb;
if nb > 0
    part = find(rand(n, 1) < 0.7);
    N.Channels = part;
    [~, order] = sort(-z1(part) + 0.5 * randn(numel(part), 1));            % broader units (higher z1) first
    rank = zeros(numel(part), 1);  rank(order) = 1:numel(part);
    N.timeStamp = cell(nb, n);
    for u=part'
        N.timeStamp(rand(nb, 1) < 0.8, u) = {1};
    end
    P(2).NetworkBurstStructure = struct('structure', [rank, rank * 5]);
else
    N.nBurst = NaN;                                                         % (rejected by the final check: no bursts)
    N.timeStamp = cell(0, n);
end
P(2).NetworkBurstDetected = N;
end

function [NS, z1] = synth(R, D, z1)
% scores of a synthetic organoid: waveform measures along D.a1 (latent z1) and D.a2, firing pattern along D.b
n = numel(z1);
T = R.WF_MEAN + R.WF_SD .* ([z1, 1.5 * randn(n, 1)] * [D.a1 D.a2]' + 0.3 * randn(n, numel(R.WF_MEAN)) + 0.1 * randn(1, numel(R.WF_MEAN)));
for m=1:numel(R.WF_MEASURES)
    switch R.WF_TRANSFORM{m}
        case 'log',      W.(R.WF_MEASURES{m}) = exp(T(:,m));
        case 'log(-x)',  W.(R.WF_MEASURES{m}) = -exp(T(:,m));
        otherwise,       W.(R.WF_MEASURES{m}) = T(:,m);
    end
end
t = linspace(-1, 1, numel(R.NOISE_SHAPE));
W.waveform = -exp(-(t / 0.2).^2) + 0.3 * sin(3 * t) + 0.02 * randn(n, numel(t));
W.nonSomatic = false(n, 1);
W.parameters.samplingRate = R.REF_SAMPLING_RATE;
Ft = R.FP_MEAN + R.FP_SD .* ((2 * randn(n, 1)) * D.b' + 0.5 * randn(n, 4));
FP.relativeRate = nan(n, 1);  FP.lv = Ft(:,2);  FP.cv = exp(Ft(:,3));  FP.burstIndex = Ft(:,4);
FP.rateOutside = 10.^(Ft(:,1) + randn);
FP.enoughSpikes = rand(n, 1) > 0.1;
FP.nSpikesOutside = 50 + 100 * FP.enoughSpikes;
FP.parameters.MIN_REFERENCE = 5;
FP = FiringPatternFeatures(FP, [], [], []);
NS = NeuronScores(FP, [], [], W, NaN);
end

