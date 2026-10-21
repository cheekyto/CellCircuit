function mergeSpikeTrains = mergeSpikeTrains(units, spikeTimestamp)

% function to merge the spike trains of the selected clusters in a single array, sorted in time (revised October 2026)
% It is used by "refinePeakDetection" to bin all the clusters at once (instead of one call of "bin_spike_time" per cluster); it is built
% once in "PopulationActivityDetection" and used for all the peaks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   units(j):           =0, if channel j should be ignored;     =1 (or >0), if channel j should be incoporated
%   spikeTimestamp:     structure array with the field .timeStamp (time stamp of all the spikes, in msec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   mergeSpikeTrains.time(i):       time stamp of the spike i, in msec (all the selected clusters, ascending order)
%   mergeSpikeTrains.unit(i):       cluster of the spike i (1..n_clusters, in the order of the selected clusters, as the "spikeData(k)=[]" loops)
%   mergeSpikeTrains.nSpikes(j):    total number of spikes of the selected cluster j (to normalize, as in "bin_spike_time")
%   mergeSpikeTrains.n_clusters:    total number of selected clusters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

selected = find(units(:)' ~= 0);                                            % selected clusters, in the original order
n_clusters = numel(selected);
trains = cell(n_clusters,1);
for j=1:n_clusters
    trains{j} = spikeTimestamp(selected(j)).timeStamp(:);
end
nSpikes = cellfun(@numel, trains);
time = vertcat(trains{:}, zeros(0,1));
unit = repelem((1:n_clusters)', nSpikes);
[time, order] = sort(time);                                                 % (the order of the spikes with the same time does not matter: only counts)

mergeSpikeTrains.time = time;
mergeSpikeTrains.unit = unit(order);
mergeSpikeTrains.nSpikes = nSpikes;
mergeSpikeTrains.n_clusters = n_clusters;

end
