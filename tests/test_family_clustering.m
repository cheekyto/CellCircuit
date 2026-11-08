function passed = test_family_clustering(code_dir)

% Function-level test of the families of network burst events ("cluster_neuronal_events", "primaryFamilies") on synthetic spike trains
% (November 2026) It does not compare two versions of the code: it checks that
%   1. the participating units of each event are the units that fire above their own rate (planted units found; units that only fire
%      at their usual rate in a long event are not participating units, as they were with the rule "at least one spike")
%   2. two planted ensembles, each repeated, are two families; the events of a third ensemble that appears only twice belong to no
%      family (label 0); the two families are fully stable when the units are resampled
%   3. same seed = same result
%   4. "primaryFamilies": the families, then one group (label 0) with the remaining events, in no family; with no family, all the primary
%      events are that single group
%
% Example (from a terminal):
%   matlab -batch "addpath('<repo>/tests'); test_family_clustering('<repo>')"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code_dir:             folder with the code to be tested

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% passed:               =true if all the checks passed (each check is printed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(char(code_dir));
passed = true;
N_UNITS = 40;
RANGE = [0 600];                                                            % window of the synthetic condition, in sec
GROUPS = {1:12, 20:30, 33:37};                                              % planted ensembles
N_REPEAT = [8 6 2];                                                         % number of events of each ensemble (the third one: too few for a family)
rs = RandStream('twister', 'Seed', 1);

% events: 0.5 s each, separated by 20 s; one long event (4 s) where no unit fires above its usual rate
group_of_event = [repelem(1:numel(GROUPS), N_REPEAT), 0];
group_of_event = group_of_event([randperm(rs, numel(group_of_event)-1), numel(group_of_event)]);
n_events = numel(group_of_event);
starts = 20*(1:n_events)' - 0.25;   lasts = starts + 0.5;   lasts(end) = starts(end) + 4;   % peaks every 20 s (see "synthetic_inputs")

% spike trains: Poisson background at 2 Hz; in each event, the units of its ensemble fire 8 extra spikes
spikes = cell(N_UNITS, 1);
for i = 1:N_UNITS
    t = RANGE(2)*rand(rs, poissrnd(2*RANGE(2)), 1);
    for k = find(group_of_event > 0)
        if ismember(i, GROUPS{group_of_event(k)})
            t = [t; starts(k) + 0.5*rand(rs, 8, 1)]; %#ok<AGROW>
        end
    end
    spikes{i} = sort(t)*1000;                                               % in msec
end

[NBP, NBF, spikeTimestamp] = synthetic_inputs(spikes, starts, lasts);
clustering = cluster_neuronal_events(NBP, NBF, spikeTimestamp, RANGE, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Participating units of each event

found = 0;  planted = 0;  false_units = 0;  tested = 0;
for k = find(group_of_event > 0)
    members = GROUPS{group_of_event(k)};
    found = found + sum(clustering.participation(k, members));
    planted = planted + numel(members);
    others = setdiff(1:N_UNITS, members);
    false_units = false_units + sum(clustering.participation(k, others));
    tested = tested + numel(others);
end
passed = check(passed, found == planted, sprintf('1a. planted units found in their events: %d of %d', found, planted));
passed = check(passed, false_units/tested < 0.02, sprintf('1b. other units counted as participating: %d of %d (%.1f%%)', false_units, tested, 100*false_units/tested));
n_long = sum(clustering.participation(end, :));
n_fired = sum(cellfun(@(t) any(t >= starts(end)*1000 & t <= lasts(end)*1000), spikes));
passed = check(passed, n_long <= 1, sprintf('1c. long event without ensemble: %d participating units (%d units fired at least one spike)', n_long, n_fired));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Families

label = clustering.label(:,1);
same_family = @(g) numel(unique(label(group_of_event == g))) == 1 && all(label(group_of_event == g) > 0);
passed = check(passed, clustering.n_cluster == 2 && same_family(1) && same_family(2) && ...
    label(find(group_of_event == 1, 1)) ~= label(find(group_of_event == 2, 1)), ...
    sprintf('2a. two families, one per repeated ensemble (found %d families)', clustering.n_cluster));
passed = check(passed, all(label(group_of_event == 3) == 0) && label(end) == 0, '2b. the events of the ensemble seen twice and the long event belong to no family (label 0)');
passed = check(passed, all(clustering.stability == 1) && clustering.clusteringQuality == 1, sprintf('2c. stability of the families: %s', mat2str(clustering.stability')));
passed = check(passed, isequal(label(1), 1) && all(diff(arrayfun(@(f) find(label == f, 1), 1:clustering.n_cluster)) > 0), '2d. families numbered in the order of their first event');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Reproducibility

clustering2 = cluster_neuronal_events(NBP, NBF, spikeTimestamp, RANGE, 0);
passed = check(passed, isequaln(clustering, clustering2), '3. same seed = same result');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. primaryFamilies

[f1, n1] = primaryFamilies([2 0 1 2 0]');
[f2, n2] = primaryFamilies([0 0 0]');
[f3, n3] = primaryFamilies([1 1 1]');
passed = check(passed, isequal(f1, [1; 2; 0]) && n1 == 2 && isequal(f2, 0) && n2 == 0 && isequal(f3, 1) && n3 == 1, ...
    '4. primaryFamilies: families, then the remaining events (label 0) as one group, not counted as a family; no family -> one group (label 0); old results (no label 0) unchanged');

if passed
    fprintf('\nAll the checks passed\n');
else
    fprintf('\nSOME CHECKS FAILED\n');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NBP, NBF, spikeTimestamp] = synthetic_inputs(spikes, starts, lasts)
% minimal inputs of "cluster_neuronal_events": primary events only; all the units are participating units; the 7 scalar metrics are
% the same for every event, so the families depend only on the units
n_events = numel(starts);   n_units = numel(spikes);
NBP.nBurst = n_events;
NBP.Secondary.nBurst = 0;
NBP.timeStamp = cell(n_events, n_units);
NBP.ChannelsIntegrated = (1:n_units)';
for k = 1:n_events
    NBP.IndivBurst(k).start = starts(k);
    NBP.IndivBurst(k).last = lasts(k);
    NBP.IndivBurst(k).peak = starts(k) + 0.25;                              % (the same interval between peaks for every event: the 7 metrics are equal)
end
metric = struct('duration', 0.5, 'nSpikes', 100, 'nNeuronsAct', 10, 'isi', 5, 'tail', 0, 'amplitudePop', 1);
NBF{1}.indivMetric = repmat(metric, n_events, 1);
spikeTimestamp = struct('timeStamp', cellfun(@(t) t(:), spikes, 'UniformOutput', false));
end

function passed = check(passed, ok, text)
if ok
    fprintf('  ok    %s\n', text);
else
    fprintf('  FAIL  %s\n', text);
    passed = false;
end
end
