function HS_spikes = HyperSyncPair(signal1, signal2, t_step)

% function to identify hypersynchronized pairs of spikes from the activity recorded in two channels of an MEA system
% (revised October 2026: same results, much faster -> binary search instead of comparing every pair of spikes; it uses "count_below")
% It is expected that individual spike trains were filtered based on the absolute refractory period

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% signal1(i,1):         time stamp of the spike i recorded in channel 1, in msec (reference channel)
% signal1(i,2):         indicates if the spike i belongs to a burst (=0, No; =1, Yes)

% signal2(i,1):         time stamp of the spike i recorded in channel 1, in msec (study channel)
% signal2(i,2):         indicates if the spike i belongs to a burst (=0, No; =1, Yes)

% t_step:               temporal resolution - sampling rate, in ms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HS_spikes.burst.nspikes_norm          proportion of spikes displaying hypersynchronicity within bursts 
% HS_spikes.burst.nspikes_hs            total number of spikes displaying hypersynchronicity within bursts, ch1 + ch2 
% HS_spikes.burst.nspikes_total         total number of spikes within bursts, ch1 + ch2 
% HS_spikes.nonburst.nspikes_norm       proportion of spikes displaying hypersynchronicity outside of bursts 
% HS_spikes.nonburst.nspikes_hs         total number of spikes displaying hypersynchronicity outside of bursts, ch1 + ch2 
% HS_spikes.nonburst.nspikes_total      total number of spikes outside of bursts, ch1 + ch2 
% HS_spikes.total.nspikes_norm          proportion of spikes displaying hypersynchronicity (entire recording)
% HS_spikes.total.nspikes_hs            total number of spikes displaying hypersynchronicity, ch1 + ch2 (entire recording)
% HS_spikes.total.nspikes_total         total number of spikes, ch1 + ch2 (entire recording)
% HS_spikes.threshold                   threshold used to identify hypersynchronized spikes, in msec
% HS_spikes.data{k}                     (i,1): time stamp of the spike i, in channel k
%                                       (i,2): # of hypersynchornized pair for spike i, channel K (=0 if it is not a HS spike)
%                                       (i,3): =1 if spike i belongs to a burst in channel k, =0 if it is outside bursts in channel k
% HS_spikes.total.ch_X.nspikes_norm     proportion of spikes displaying hypersynchronicity in channel X (1 or 2)
% HS_spikes.total.ch_X.nspikes_hs       total number of spikes displaying hypersynchronicity in channel X (1 or 2)
% HS_spikes.total.ch_X.nspikes_total    total number of spikes in channel X (1 or 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JITTER = 3;                             % threshold to define the hypersynchronicity, relative to t_step (default: 3.0)
MAX_JITTER = 1.0;                       % maximum threshold to define the hypersynchronicity, in msec (default: 1.0)
apply_parameters('HyperSyncPair');                                          % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify all hypersynchronized pair of spikes

n_spikes1 = size(signal1,1);    n_spikes2 = size(signal2,1);                % total number of spikes in channel 2
data1 = zeros(n_spikes1,2);     data2 = zeros(n_spikes2,2);                 % initialize the variable data
data1(:,1) = signal1(:,1);      data2(:,1) = signal2(:,1);                  % data(i,1): time stamp of the spike i, in msec                                
threshold = min(MAX_JITTER,JITTER*t_step);                                  % set the threshold for hypersynchronicity (it can NOT be larger than MAX_JITTER), in msec
HS_spikes.threshold = threshold;                                            % load the threshold into the output variable

% (revised October 2026) Same result as the original loop across all spikes in channel 1, where the distance to ALL spikes in channel 2 was
% calculated each time (n_spikes1 x n_spikes2 operations, plus a "clear" that took most of the time). Now, for each spike in channel 1,
% the few candidate spikes in channel 2 are found by binary search, and the original criterion abs(t2-t1)<threshold is checked only for them
t1 = signal1(:,1);                                                          % time stamp of the spikes in channel 1, in msec
[t2, order2] = sort(signal2(:,1));                                          % time stamp of the spikes in channel 2 (sorted, in case they were not), in msec
margin = max(1e-6*threshold, 16*eps(max(abs([t1; t2; 1]))));               % small margin, so that no candidate is missed because of round-off (the exact criterion is checked below)
first2 = count_below(t2, t1 - threshold - margin, 1) + 1;                   % first candidate in channel 2 for each spike i in channel 1 (index in the sorted channel 2)
last2 = count_below(t2, t1 + threshold + margin, 0);                        % last candidate in channel 2 for each spike i in channel 1
n_cand = max(last2 - first2 + 1, 0);                                        % number of candidates for each spike i in channel 1 (usually 0, 1 or 2)
idx1 = repelem((1:n_spikes1)', n_cand);                                     % all pairs (spike in channel 1, candidate spike in channel 2)
idx2 = repelem(first2, n_cand) + ((1:sum(n_cand))' - repelem(cumsum(n_cand)-n_cand, n_cand) - 1);
is_hs = abs(t2(idx2) - t1(idx1)) < threshold;                               % original criterion: the two spikes are too close (hypersynchronized)
idx1 = idx1(is_hs);     idx2 = order2(idx2(is_hs));                         % hypersynchronized pairs (index of the spike in channel 2, back in its original order)

has_pair = false(n_spikes1,1);      has_pair(idx1) = true;                  % spikes in channel 1 with at least one spike in channel 2 that is too close
pair_id = cumsum(has_pair);                                                 % id of the hypersynchronized pair, which increases with i (as in the original loop)
n_HSpairs = sum(has_pair);                                                  % total number of hypersynchronized pairs (it does not matter if it is one or more spikes in channel 2)
data1(has_pair,2) = pair_id(has_pair);                                      % label the spike i from channel 1, with the id of the pair
last1 = accumarray(idx2(:), idx1(:), [n_spikes2 1], @max);                  % for each spike in channel 2, the LAST spike in channel 1 that is too close (as in the original loop, where the id was overwritten)
data2(last1>0,2) = pair_id(last1(last1>0));                                 % label the spikes from channel 2 with the same id (shared by a single spike in channel 1 and all spikes in channel 2 that are too close to it)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify the bursts in both channels

data1(:,3) = signal1(:,2);                                                  % data1(i,3) indicates if the spike i in channel 1 belongs to a burst in channel 1
data2(:,3) = signal2(:,2);                                                  
HS_spikes.data{1} =data1;                                                   % load the variable "data" into the output variable
HS_spikes.data{2} =data2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Count the number of hypersynchronized spikes 

n_total = n_spikes1 + n_spikes2;                                            % total number of spikes across channel 1 and 2
aux1_hs = find(data1(:,2)>0);                                               % index of hypersynchronized spikes in channel 1
aux2_hs = find(data2(:,2)>0);                                               % index of hypersynchronized spikes in channel 2
n_total_burst = numel(find(data1(:,3)>0)) + numel(find(data2(:,3)>0));      % total number of spikes within bursts (either in channel 1 or 2)
n_total_nonburst = numel(find(data1(:,3)==0)) + numel(find(data2(:,3)==0)); % total number of spikes outside bursts (either in channel 1 or 2)

if numel(aux1_hs)+numel(aux2_hs) > 0 && n_HSpairs > 0                       % there is at least one pair of hypersynchronized spikes
    aux1_hs_burst = find(data1(aux1_hs,3)>0);                               % index of hypersynchronized spikes within burst in channel 1
    aux2_hs_burst = find(data2(aux2_hs,3)>0);                               % index of hypersynchronized spikes within burst in channel 2
    n_burst_synch = numel(aux1_hs_burst)+numel(aux2_hs_burst);              % total number of hypersynchronized spikes within bursts
    aux1_hs_nonburst = find(data1(aux1_hs,3)==0);                           % index of hypersynchronized spikes outside bursts in channel 1
    aux2_hs_nonburst = find(data2(aux2_hs,3)==0);                           % index of hypersynchronized spikes outside bursts in channel 2
    n_nonburst_synch = numel(aux1_hs_nonburst)+numel(aux2_hs_nonburst);     % total number of hypersynchronized spikes outside bursts
else
    n_burst_synch = 0;                                                      % no hypersynchronized spikes within bursts
    n_nonburst_synch = 0;                                                   % no hypersynchronized spikes outside bursts
end

HS_spikes.burst.nspikes_norm = n_burst_synch/n_total_burst;                 % it is NaN when there is no bursting in any of both channels     
HS_spikes.burst.nspikes_hs = n_burst_synch;
HS_spikes.burst.nspikes_total = n_total_burst;
HS_spikes.nonburst.nspikes_norm = n_nonburst_synch/n_total_nonburst;
HS_spikes.nonburst.nspikes_hs = n_nonburst_synch;
HS_spikes.nonburst.nspikes_total = n_total_nonburst;
HS_spikes.total.nspikes_norm = (numel(aux1_hs)+numel(aux2_hs))/n_total;
HS_spikes.total.nspikes_hs = numel(aux1_hs)+numel(aux2_hs);
HS_spikes.total.nspikes_total = n_total;

HS_spikes.total.ch_1.nspikes_norm = numel(aux1_hs)/n_spikes1;
HS_spikes.total.ch_1.nspikes_hs = numel(aux1_hs); 
HS_spikes.total.ch_1.nspikes_total = n_spikes1;
HS_spikes.total.ch_2.nspikes_norm = numel(aux2_hs)/n_spikes2;
HS_spikes.total.ch_2.nspikes_hs = numel(aux2_hs); 
HS_spikes.total.ch_2.nspikes_total = n_spikes2;

end

