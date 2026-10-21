function QualityControl_General(inputs)

% function to run Quality Control assessment in a MEA recording, based on general features of the activity pattern
% It removes clusters based on: 1) too many spikes that violate the refractory period criteria, 2) excessive (tonic) firing level, 3) too low level of firing, 4) responsiveness to TTX treatment, 5) cropping of the electrode map based on a picture
% (revised October 2026: faster removal of echo spikes; fixed the check that rescues the outliers with too much activity ("outlier2");
%  no crash when no cluster is discarded or when all clusters are silent; "discarded_signal" is now a column cell array)
% (optional) inputs:    structure to overwrite any variable of the sections INPUTS and PARAMETERS (e.g., inputs.folder, inputs.option_plot); see "override_inputs"
%                       when it is not provided, the values defined in the code are used

if ~parameter_set('active'), clc; end; close all;      % (revised October 2026) no "clc" while "run_pipeline" runs (it would erase its progress messages)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system = 3;                 % =1, Ubuntu; =2, Yang's PC; =3, Yang's MacBook; =0, do not modify the MATLAB path (e.g., test harness)
if nargin < 1, inputs = struct(); end                                       % no inputs provided -> use the values defined in the code
if isfield(inputs, 'system'), system = inputs.system; end                   % the system can also be provided as an input
switch system
    case 1      % Anlysis PCP - Ubuntu
        addpath(genpath('/home/cortex/Documents/Yang'));
    case 2      % Yang's PC - GitHub
        addpath(genpath('C:\Users\Yang\OneDrive - Harvard University\Documents\GitHub\MEA-analysis\'));
    case 3      % Yang's MacBook
        addpath(genpath('/Users/sungminyang/Documents/GitHub/CellCircuit'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = '/path/to/recording_folder';

condition = -1;                                     % =-1, the entire recording (including all conditions)
                                                    % otherwise, the list of conditions (it must be in order)
option_plot = 1;                                    % = 1, plot the raster plot and some metrics
                                                    % = 0, do not plot (just save the output)
option_pause = 1;                                   % = 1, wait for a key press before moving on to the next file
                                                    % = 0, do not wait (e.g., batch processing)
condition_3Brain = 1 ;                               % = 0, the data does not come directly from 3Brain (e.g., spike sorting done by Kilosort)
                                                    % = 1, the data comes directly from the 3Brain platform (no intermediate spike sorting) 
                                                    
option_picture = 0;                                 % there is a picture to estimate the location of the organoid
if  option_picture                                  % crop the electrode map using an ellipse with the following features
    x0 = 34;                                        % coordinate 'x' of the center, in pixels (left to right)
    y0 = 42;                                        % coordinate 'y' of the center, in pixels (bottom to up)  
    radiusx = 23;                                   % radius in the axis 'x', in pixels
    radiusy = 23;                                   % radius in the axis 'y', in pixels
    ellipse_cropped = [x0 y0 radiusx radiusy];      % ellipse to be cropped
    CROPPING_SAFETY = 1.10;                         % proportion of the actual ellipse to be cropped (default: 1.10)
else
    ellipse_cropped = nan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REFRACTORY_PERIOD = 1.0;                        % absolute refractory period for neuronal firing, in msec (default: 1.0)
REF_PERIOD_THRESHOLD = 5;                       % threshold for the percentage of spikes that violate the refractory period criteria; above this threshold, the cluster will be discarded (default: 5%)
REF_PERIOD_MINSPIKES = 3;                       % threshold for the absolute number of spikes that violate the refractory period criteria (default: 3); to account for the low activity units

ACTIVITY_THRESHOLD = 0.03;                      % minimum level of spiking frequency for an individual electrode to be considered active, in Hz (0.02 Hz = 1/50 sec; 0.03 Hz = 1/33 sec)
CV_THRESHOLD = 3.5;                             % threshold for coefficient of variation to separate phasic vs tonic activity (default: 3.5)
LVR_THRESHOLD = 1.00;                           % threshold for local variation ratio to separate phasic vs tonic activity (default: 1.00)
PERCENTILE = 98;                                % percentile to identify the outliers with too high firing rate, among ONLY the clusters with tonic activity (default: 98)
PERCENTILE_TOP = 99.5;                          % percentile to identify the outliers with too high firing rate, among ALL clusters (default: 99.5)
ISI_RATIO = 8;                                  % threshold for the ratio of isi during bursting vs non-bursting; the lower, the more permissive (default: 8)
N_SPIKES_BURST = 0.90;                          % threshold for the proportion of spikes in bursts; if the cluster is mostly bursting, probably is a neuron (default: 0.90)

TTX_EFFECT = 0.15;                              % minimum effect of TTX to be considered; proportion relative to the maximum frequency (0.10 - 0.20)

BURST_METHOD = 1;                               % method for the burst detection (1)
MIN_LENGTH = 90;                                % minimum duration of a recording condition, in sec (default: 90)
N_PLOT = 20;                                    % number of clusters for raster plots (20)
SHIFT=0.1;                                      % distance between rows in the raster plots (0.1)
LABEL_POS = 0.6;                                % position for a message in the raster plot (0.6)
MATRIX_SIZE_X = 64;                             % horizontal length of the electrode array (64*64=4096)
MATRIX_SIZE_Y = 64;                             % vertical length of the electrode array (64*64=4096)
K_MSEC_SEC = 1/1000;                            % factor to rescale from msec to sec

apply_parameters('QualityControl_General');                                 % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")
override_inputs(inputs);                        % overwrite the INPUTS / PARAMETERS above with the ones provided as input (if any)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save the following variables in the file sprintf('%s\\%s_cleaned.mat', folder,file_data):
%   - structure array spikeTimestamp (same structure than the output generated by "ShapeSpikeTrainData"), ONLY with good clusters   
%   - removedSpikes(k,1):       percentage of spikes removed for violating the refractory period in cluster k; since it refers to the originally set of clusters, it needs to be used together with the variable "spikeTimestamp" generated by "ShapeSpikeTrainData"
%   - echo_cluster(k,1):        it indicates if cluster k was discarded for having too many echo spikes (=1, discarded; =0, good)
%   - silent_cluster(k,1):      it indicates if cluster k was discarded for having low level of activity (=1, discarded; =0, good)
%   - noisy_cluster(k,1):       it indicates if cluster k was discarded for having high level of tonic activity (=1, discarded; =0, good)
%   - ttx_cluster(k,1):         it indicates if cluster k showed effect upon TTX treatment
%   - picture_cluster(k,1):     it indicates if cluster k where positively cropped based on a picture of the setup
%   - discarded_cluster(k,1):   it indicates if cluster k was finally discarded (=1, discarded; =0, good)
%   - n_discarded:              total number of discarded clusters (independent of the reason)
%   - n_cluster0:               total number of originally identified clusters
%   - ellipse_cropped:          parameters that define the ellipse to be cropped
%   - discarded_signal:         time stamp of units from discarded clusters    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if condition_3Brain
    file_data = dir(fullfile(folder, '*_3Brain_raw.mat'));                  % files in the folder that went through "ShapeSpikeTrainData"
else
    file_data = dir(fullfile(folder, '*_raw.mat'));                         % files in the folder that went through "ShapeSpikeTrainData"
end
n_files = size(file_data,1);                                                % # files to assess
if n_files<1
    error('There are no data files in this folder!!!')
end

for i_file = 1:n_files                                                      % across all data files in the folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data

fprintf(sprintf("\n Processing file %s ... \n\n", file_data(i_file).name));

file_input = fullfile(folder, sprintf('%s', file_data(i_file).name));    % path to the data file
load(file_input, 'Indiv_Cond');       
if condition_3Brain
    Indiv_Cond_3Brain = Indiv_Cond;
    clear Indiv_Cond
    for j=1:size(Indiv_Cond_3Brain,1)
        Indiv_Cond{j} = Indiv_Cond_3Brain{j,1};                             % the second and third column contain the waveform information
    end
end 
if condition==-1                                                            % take the whole recording (across all diferent conditions)
    load(file_input, 'spikeTimestamp','rec_conditions');       % load the original data
    n_cluster = size(spikeTimestamp,2);                                     % total number of clusters
    n_conditions = size(rec_conditions,1);                                  % total number of conditions to be extracted

else                                                                        % take only the conditions listed in "condition" 
    load(file_input, 'rec_conditions');                        % load the original data
    original_rec_conditions = rec_conditions;
    original_Indiv_Cond = Indiv_Cond;
    clear rec_conditions Indiv_Cond
    n_conditions = numel(condition);                                        % total number of conditions to be extracted
    rec_conditions = cell(n_conditions,1);
    
    spikeTimestamp = original_Indiv_Cond{condition(1)};                     % initialize with the first condition
    rec_conditions(1) = original_rec_conditions(condition(1));              % label of the first condition
    Indiv_Cond{1} = original_Indiv_Cond{condition(1)};                      % load the data for the first condition
    n_cluster = size(spikeTimestamp,2);                                     % total number of clusters/units
    if n_conditions>1                                                                   % when there is more than one recording condition
        for cond=2:n_conditions                                                         % across all recording conditions
            rec_conditions(cond) = original_rec_conditions(condition(cond));            % label of the condition "cond"
            Indiv_Cond{cond} = original_Indiv_Cond{condition(cond)};                    % data for the condition "cond"    
            for i=1:n_cluster
                if condition_3Brain
                    shift = spikeTimestamp(i).recLength/K_MSEC_SEC;
                else
                    shift = 0;
                end
                spikeTimestamp(i).window(2) = Indiv_Cond{cond}(i).window(2);                                        % extend the end of the window    
                spikeTimestamp(i).recLength = spikeTimestamp(i).recLength + Indiv_Cond{cond}(i).recLength;          % add the duration of the condition cond
                spikeTimestamp(i).nUnits = spikeTimestamp(i).nUnits + Indiv_Cond{cond}(i).nUnits;                   % add the number of units of the condition cond
                spikeTimestamp(i).timeStamp = vertcat(spikeTimestamp(i).timeStamp, Indiv_Cond{cond}(i).timeStamp+shift);  % append the spike time recorded during condition cond

            end
        end
    end
end
spike_plot = spikeTimestamp;                                                % keep the original data to plot

fprintf('The conditions to be extracted are: \n')
for i=1:n_conditions
    fprintf('%s \n', string(rec_conditions(i)));
    if Indiv_Cond{i}(1).recLength < MIN_LENGTH                              % if the recording condition is not long enough in duration
        error('This condition is too short in duration!! \n')
    end
end
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Eliminate the spikes that violate the refractory period and identify the clusters with too many of them

echo_cluster = zeros(n_cluster,1);
removedSpikes = zeros(n_cluster,1);

for i=1:n_cluster                                                           % across all identified clusters
    isi_indv = diff(spikeTimestamp(i).timeStamp);                           % inter-spike interval for cluster i (along all recording conditions)       
    time_indiv = spikeTimestamp(i).timeStamp(1:end);                        % start of each ISI 
    aux = find(isi_indv < REFRACTORY_PERIOD);                               % identify the pairs of spikes closer than the refractory period.  
    n_tooClose = numel(aux);                                                % total number of spikes that violate the refractory period

    if n_tooClose>0  
        removedSpikes(i) = n_tooClose*100/spikeTimestamp(i).nUnits;         % percentage of echo spikes in cluster i
        if removedSpikes(i) > REF_PERIOD_THRESHOLD && n_tooClose > REF_PERIOD_MINSPIKES  % too many close spikes 
            echo_cluster(i) = 1;                                            % label the echo cluster (the cluster will be discarded)
        end
        keep = true(size(time_indiv));                                      % spikes to keep (revised October 2026: the second spikes are deleted all together at the end,
                                                                            % instead of one by one inside the loop; same result, since the pairs are processed in reverse order)
        for k=n_tooClose:-1:1                                               % across all close pair of spikes, in reverse order
            first = aux(k);                                                 % index of the first spike of the pair
            second = aux(k)+1;                                              % index of the second spike of the pair
            time_indiv(first) = (time_indiv(first)+time_indiv(second))/2;   % replace the first spike by the mid-point
            keep(second) = false;                                           % and delete the second spike
        end
        spikeTimestamp(i).timeStamp = time_indiv(keep);                     % replace the time stamp sequence, after removing the echo spikes
        spikeTimestamp(i).nUnits = spikeTimestamp(i).nUnits - n_tooClose;   % updated number of total spikes
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify the active clusters / discard the clusters with too low level of activity

active_cluster = zeros(n_cluster,1);                                        % clusters with enough level of activity      
silent_cluster = zeros(n_cluster,1);                                        % clusters with extremely low level of activity
ttx_cluster = zeros(n_cluster,1);                                           % clusters with response to TTX

cond_ttx = [];
ttx_n = 0;

for i=1:n_conditions
    if contains(rec_conditions{i},'ttx','IgnoreCase',true)
        ttx_n = ttx_n+1;
        cond_ttx(ttx_n) = i;                                              % indexes of the recording conditions that contain the TTX treatment
    end
end

if ttx_n==0
    fprintf('There is NO treatment with TTX \n');
else
    if ttx_n>1
        fprintf('The files that contain the treatment with TTX are: \n');
        for t=1:ttx_n
            fprintf('%s\n', rec_conditions{cond_ttx(t)});
        end
    else
        fprintf('The file that contains the treatment with TTX is: %s\n', rec_conditions{cond_ttx});
    end
end

for k=1:n_cluster                                                           % across all clusters
    freq = nan(n_conditions,1);
    for i=1:n_conditions                                                    % across all recording conditions
        freq(i) = double(Indiv_Cond{i}(k).nUnits)/double(Indiv_Cond{i}(k).recLength);       % rough spiking frequency for cluster k, in 1/sec
    end
    aux = find(freq > ACTIVITY_THRESHOLD);                                  % conditions in which the cluster k has enough activity level
    if sum(aux)>0                                                           % at least one condition with enough activity
       active_cluster(k) = 1;                                               % label the clusters with suprathreshold level of activity 
    end

    max_freq = max(freq);                                                   % maximum frequency reached by the cluster k (in any of the recording conditions)
    if ttx_n>0 && sum(max_freq*TTX_EFFECT > freq(cond_ttx))>0               % there is at least one TTX condition with activity level below TTX_EFFECT*max_freq
       ttx_cluster(k) = 1;                                                  % label the clusters with positive TTX effect
    end
    clear freq max_freq
end
silent_cluster(active_cluster==0) = 1;                                      % =1, "silent";  =0, active cluster

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify the channels with abnormal activity pattern                     % cluster with tonic firing and high level of activity 
                                                                            % at least the constrast between phasic vs tonic is super clear                                                                            
nUnits = nan(n_cluster,1);
cvUnits = zeros(n_cluster,1);
lvrUnits = zeros(n_cluster,1);

for i=1:n_cluster                                                           % across all clusters (active & "silent")
    nUnits(i) = spikeTimestamp(i).nUnits;                                   % total number of units in cluster i (NaN for "silent" cluster)                
    if nUnits(i)>1                                                          % redundant
       isi = diff(spikeTimestamp(i).timeStamp);                             % isi for all spikes in cluster i (across all recording conditions)
       cvUnits(i) = std(isi)/mean(isi);                                     % coefficient of variation for isi in cluster i (NaN for "silent" cluster)  
       lvrUnits(i) = lvr(isi);                                              % local variation ratio in cluster i (NaN for "silent" cluster)  
       clear isi
    end
end

uniform = find(cvUnits < CV_THRESHOLD & lvrUnits < LVR_THRESHOLD & silent_cluster < 1);          % clusters that displayed clear tonic/uniform firing (the CV and LVR are below threshold), among the active clusters
if numel(uniform)>0
    aux_nspikes = nUnits(uniform);                                          % total number of units for the tonic firing clusters
    outlier = isoutlier(aux_nspikes, "percentile",[0 PERCENTILE]);          % identify the clusters with high level of activity (only among tonic firing clusters)
    if sum(outlier)>0
        idx_outlier1 = uniform(find(outlier));                              % index of all outliers
        n_outlier = sum(outlier);                                           % total number of outliers
    else                                                                    % there is no cluster with high level of activity
        idx_outlier1 = [];
        n_outlier = 0;
    end
else                                                                        % the CV/LVR are above threshold for all active clusters
    idx_outlier1 = [];
    n_outlier = 0;
end

nonSilent = find(silent_cluster < 1);
outlier2 = [];                                                              % (revised October 2026) initialized, so that the code does not crash when all clusters are "silent"
idx_outlier = idx_outlier1;
if numel(nonSilent)>0
    aux_nspikes = nUnits(nonSilent);                                        % total number of units for the clusters with suprathreshold level of activity 
    outlier2_aux = isoutlier(aux_nspikes, "percentile",[0 PERCENTILE_TOP]); % identify the cluster with too much activity
    outlier2 = nonSilent(find(outlier2_aux));                               % index of the clusters with too much activity (cluster indexes, NOT a logical array)
    % outlier2 = isoutlier(nUnits, "percentile",[0 PERCENTILE_TOP]);              % identify the cluster with too much activity
    if sum(outlier2)>0
%        aux = vertcat(idx_outlier1, find(outlier2));                            % pool the clusters with: 1) tonic (high) activity & 2) too much activity
        aux = vertcat(idx_outlier1, outlier2);                            % pool the clusters with: 1) tonic (high) activity & 2) too much activity
        idx_outlier = unique(aux);                                              % id of clusters with either: 1) tonic (high) activity & 2) too much activity
        n_outlier = numel(idx_outlier);                                         % total number of outlier clusters
    else
        idx_outlier = idx_outlier1;                                              
        n_outlier = numel(idx_outlier);                                         
    end
end
% last check before discarding the outlier clusters
for i=n_outlier:-1:1                                                                        % across all outliers, in reverse order
    BurstDetected = BurstDetection(spikeTimestamp(idx_outlier(i)).timeStamp, BURST_METHOD); % call the function to detect bursts
    BurstFeature = BurstFeatures(BurstDetected);                                            % call the function to characterize the bursts
    isi_burst = BurstFeature.neuronMetric.isi.mean;                                         % mean inter-spike interval within bursts
    isi_nonburst = BurstFeature.neuronMetric.nonBursting.isi.mean;                          % mean inter-spike interval outside bursts    
    nSpikes_burst = BurstFeature.neuronMetric.nSpikes.norm;                                 % total number of spikes within burst normalized by the total number of spikes in the recording (ratio of spikes participating in bursting)
                                                                                            % it won't be discarded if it shows clear contrast between tonic and phasic activity
                                                                                            % delete the potential outlier from the black list when the firing rate during burst is ISI_RATIO larger than during tonic activity 
                                                                                            % or when most of the spikes are bursting
                                                                                            % the clusters with too much activity ("outlier2") can only be rescued if most of their spikes are bursting
                                                                                            % (revised October 2026: it was "find(outlier2)", which returned 1:numel(outlier2) instead of the cluster indexes,
                                                                                            %  since "outlier2" is no longer a logical array)
    if (isi_burst*ISI_RATIO < isi_nonburst && ~ismember(idx_outlier(i),outlier2)) || nSpikes_burst > N_SPIKES_BURST
       idx_outlier(i) = [];                                                                 % the contrast between tonic and phasic firing is clear, suggesting that it is a real neuron and not a noisy channel
    end
    clear BurstDetected BurstFeature isi_burst isi_nonburst nSpikes_burst
end

noisy_cluster = zeros(n_cluster,1);
noisy_cluster(idx_outlier) = 1;                                             % =1, noisy;  =0, normal cluster

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cropped the electrode map based on a picture

if option_picture                                                                       % if there is a picture to cropped the electrode map
   picture_cluster = zeros(n_cluster,1);                                                % clusters within the picture 
   ellipse = @(x,y,x0,y0,radius_x,radius_y) ((x-x0)/radius_x)^2 + ((y-y0)/radius_y)^2;  % function to define the size and location of the ellipse
   if isnan(ellipse_cropped)
       error('\n The ellipse for cropping has not been defined')
   else                                                                                 % parameters to define the size and location of the ellipse
       x0 = ellipse_cropped(1); 
       y0 = ellipse_cropped(2);
       radiusx = ellipse_cropped(3);
       radiusy = ellipse_cropped(4);
   end
   for i=1:n_cluster                                                                    % across all clusters
       pos_x = spikeTimestamp(i).positionXY(1)+1;                                       % position x of the cluster
       pos_y = spikeTimestamp(i).positionXY(2)+1;                                       % position y of the cluster
       if ellipse(pos_x,pos_y,x0,y0,radiusx*CROPPING_SAFETY,radiusy*CROPPING_SAFETY)<=1 % if the cluster is inside the ellipse    
           picture_cluster(i) = 1;                                                      % positive clusters
       end
   end
else
    picture_cluster = ones(n_cluster,1);                                                % when there is no picture, all the clusters are included
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discard the problematic clusters and Save the cleaned data

fprintf('\nTotal number of clusters originally identified by the spike sorting algorithm: %g\n', n_cluster);
fprintf('Clusters with too many echo spikes: %g\n', sum(echo_cluster));
fprintf('Clusters with low level of activity: %g\n', sum(silent_cluster));
fprintf('Clusters with high level of tonic activity: %g\n', sum(noisy_cluster));
fprintf('Clusters with TTX effect: %g\n', sum(ttx_cluster));

if ttx_n > 0                                                                                                % there is a TTX treatment   
    discarded_cluster = echo_cluster + silent_cluster + noisy_cluster + ~ttx_cluster + ~picture_cluster;    % pool all the clusters to be discarded      
else
    discarded_cluster = echo_cluster + silent_cluster + noisy_cluster + ~picture_cluster;                   % pool all the clusters to be discarded
end

idx_discarded = find(discarded_cluster>0);                                  % index of clusters to be discarded (=0 only for the clusters that checked all the QC)
n_discarded = numel(idx_discarded);                                         % total number of clusters to be discarded
fprintf('Total number of clusters to be discarded: %g (out of %g)\n', n_discarded, n_cluster);

discarded_signal = cell(n_discarded,1);                                     % (revised October 2026) column cell array (it was n_discarded x n_discarded),
                                                                            % initialized before the "if", so that it exists when no cluster is discarded (the "save" crashed)
if n_discarded > 0
   for k=n_discarded:-1:1                                                   % across all clusters to be discarded
       discarded_signal{k}=spikeTimestamp(idx_discarded(k)).timeStamp;      % time stamp of units from discarded clusters    
       spikeTimestamp(idx_discarded(k))=[];                                 % eliminate the cluster with abnormal pattern of activity
   end
end

n_cluster0 = n_cluster;
aux_file = split(file_data(i_file).name, '_raw');
file_output = aux_file{1};
file_save = fullfile(folder, sprintf('%s_cleaned.mat',file_output));              
save(file_save,'spikeTimestamp','spike_plot','rec_conditions','Indiv_Cond','removedSpikes','echo_cluster','silent_cluster','noisy_cluster','ttx_cluster','discarded_cluster','n_discarded','n_cluster0','picture_cluster','ellipse_cropped','discarded_signal', 'cond_ttx');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the result of the quality control

if option_plot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1)                                                              % plot number of spikes vs ISI CV, for all clusters
    
    plot(nUnits, cvUnits, 'ok', 'DisplayName','all clusters');    hold on

    if numel(find(echo_cluster)) > 0
        plot(nUnits(find(echo_cluster)), cvUnits(find(echo_cluster)), 'sg','DisplayName', sprintf('echo (%g)',numel(find(echo_cluster))));  % label the echo clusters that were eliminated
    else
        plot([], [], 'sg','DisplayName', sprintf('echo (%g)',numel(find(echo_cluster)))); 
    end
    if numel(find(silent_cluster)) > 0
        plot(nUnits(find(silent_cluster)), cvUnits(find(silent_cluster)), 'pm','DisplayName',sprintf('silent (%g)',numel(find(silent_cluster)))); % label the silent clusters that were eliminated
    else
        plot([], [], 'pm','DisplayName',sprintf('silent (%g)',numel(find(silent_cluster)))); 
    end
    if numel(find(noisy_cluster)) > 0
        plot(nUnits(find(noisy_cluster)), cvUnits(find(noisy_cluster)), '*r','DisplayName',sprintf('noisy (%g)',numel(find(noisy_cluster)))); % label the noisy clusters that were eliminated
    else
        plot([], [], '*r','DisplayName',sprintf('noisy (%g)',numel(find(noisy_cluster)))); 
    end
    if numel(find(ttx_cluster)) > 0
        plot(nUnits(find(ttx_cluster)), cvUnits(find(ttx_cluster)), '.c','DisplayName',sprintf('TTX effect (%g)',numel(find(ttx_cluster))));  % label the clusters with TTX effect
    else
        plot([], [], '*r','DisplayName',sprintf('TTX effect (%g)',numel(find(ttx_cluster)))); 
    end
    legend;     xlabel('Number of spikes');    ylabel('ISI coefficient of variation');
    title(sprintf('All neurons (%g discarded / %g total)', n_discarded, n_cluster));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2)                                                              % plot number of spikes vs ISI local variation ratio, for all clusters
    
    plot(nUnits, lvrUnits, 'ok', 'DisplayName','all clusters');    hold on

    if numel(find(echo_cluster)) > 0
        plot(nUnits(find(echo_cluster)), lvrUnits(find(echo_cluster)), 'sg','DisplayName', sprintf('echo (%g)',numel(find(echo_cluster))));  % label the echo clusters that were eliminated
    else
        plot([], [], 'sg','DisplayName', sprintf('echo (%g)',numel(find(echo_cluster)))); 
    end
    if numel(find(silent_cluster)) > 0
        plot(nUnits(find(silent_cluster)), lvrUnits(find(silent_cluster)), 'pm','DisplayName',sprintf('silent (%g)',numel(find(silent_cluster)))); % label the silent clusters that were eliminated
    else
        plot([], [], 'pm','DisplayName',sprintf('silent (%g)',numel(find(silent_cluster)))); 
    end
    if numel(find(noisy_cluster)) > 0
        plot(nUnits(find(noisy_cluster)), lvrUnits(find(noisy_cluster)), '*r','DisplayName',sprintf('noisy (%g)',numel(find(noisy_cluster)))); % label the noisy clusters that were eliminated
    else
        plot([], [], '*r','DisplayName',sprintf('noisy (%g)',numel(find(noisy_cluster)))); 
    end
    if numel(find(ttx_cluster)) > 0
        plot(nUnits(find(ttx_cluster)), lvrUnits(find(ttx_cluster)), '.c','DisplayName',sprintf('TTX effect (%g)',numel(find(ttx_cluster))));  % label the clusters with TTX effect
    else
        plot([], [], '*r','DisplayName',sprintf('TTX effect (%g)',numel(find(ttx_cluster)))); 
    end
    legend;     xlabel('Number of spikes');    ylabel('ISI local variation ratio');
    title(sprintf('All neurons (%g discarded / %g total)', n_discarded, n_cluster));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(3)                                                               % raster plot of representative clusters for each category

    subplot(3,1,1);                                                         % raster plot for "good" clusters
    n_good = min(N_PLOT,size(spikeTimestamp,2));                            % number of good clusters to plot
    plot_good = 1+fix((size(spikeTimestamp,2)-1)*rand(1,n_good));           % index of good clusters to plot (randomly generated)
    for i=1:n_good                                                         
        plot(spikeTimestamp(plot_good(i)).timeStamp*K_MSEC_SEC, ones(spikeTimestamp(plot_good(i)).nUnits,1)+(i-1)*SHIFT,'|',"MarkerEdgeColor",[0 0 0]); hold on
    end
    if ttx_n > 0                                                            % if there is TTX treatment
       ini=0;
       for cond = 1:n_conditions                                            % across all recording conditions
           aux = find(cond_ttx==cond);                                      % is the condition "cond" a TTX treatment?
           if numel(aux)>0                                                  % Yes, the current condition is a TTX treatment
               fin = ini + Indiv_Cond{cond}(1).recLength;                   % end of the recording condition "cond"
               plot([ini ini], [0.5 i*SHIFT+1.5],'--r'); hold on            % draw a line at the start of the TTX treatment
               plot([fin fin], [0.5 i*SHIFT+1.5],'--r'); hold on            % draw another line at the end of the TTX treatment
               text(ini,LABEL_POS,'TTX','Color','red')
           end
           ini = ini + Indiv_Cond{cond}(1).recLength;                       % accumulate the windows
       end
       xlim([0 ini]);
    else
        text(0,LABEL_POS,'No TTX treatment','Color','green')
    end
    ylim([0.5 i*SHIFT+1.5]);      xlabel('Time (sec)');    title('Good units')

    subplot(3,1,2);                                                         % raster plot for "silent" clusters
    idx_silent=find(silent_cluster);
    if numel(idx_silent)>0                     
       for k=1:min(N_PLOT,numel(idx_silent))
           plot(spike_plot(idx_silent(k)).timeStamp*K_MSEC_SEC, ones(spike_plot(idx_silent(k)).nUnits,1)+(k-1)*SHIFT,'|',"MarkerEdgeColor",[0 0 0]); hold on
       end
       ylim([0.5 k*SHIFT+1.5]);        xlabel('Time (sec)');        title('Discarded units - Low frequency')
    end

    subplot(3,1,3);                                                         % raster plot for the "noisy" clusters
    idx_noisy=find(noisy_cluster);
    if numel(idx_noisy)>0                     
       for k=1:min(N_PLOT,numel(idx_noisy))
           plot(spike_plot(idx_noisy(k)).timeStamp*K_MSEC_SEC, ones(spike_plot(idx_noisy(k)).nUnits,1)+(k-1)*SHIFT,'|',"MarkerEdgeColor",[0 0 0]); hold on
       end
       ylim([0.5 k*SHIFT+1.5]);        xlabel('Time (sec)');        title('Discarded units - High frequency tonic firing')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(4)                                                               % Electrode map showing the location of all "good" clusters

    map = nan(MATRIX_SIZE_X ,MATRIX_SIZE_Y ); 
    for i=1:n_cluster                                                       % across all clusters ("good" & discarded)
        pos_x = spike_plot(i).positionXY(1)+1;                              % position x of the cluster
        pos_y = spike_plot(i).positionXY(2)+1;                              % position y of the cluster
        label = 1;                                                          % label for the electrode (default = 0, "good" cluster)
        if discarded_cluster(i) > 0                                         % cluster to be discarded
           label = -1;                                                      % change the label (=-1, to be discarded)
        end
        if ttx_cluster(i) == 1                                              % cluster responsive to TTX
            label = 0;
            if discarded_cluster(i) <1
               label = 1;
            end
        end
        map(MATRIX_SIZE_Y - pos_y + 1,pos_x) = label;                       % label the electrode in the map
    end
    h = heatmap(map,'XDisplayLabels',NaN*ones(size(map,1),1),'YDisplayLabels',NaN*ones(size(map,1),1),'Colormap',cool);
    h.MissingDataColor = [0.5 0.5 0.5];
    title('Good electrodes (=1), TTX-sensitive electrodes (=1) & discarded electrodes (=2)')

    if  option_picture                                                  % if there is a picture to cropped the electrode map
        pos = get(gca, 'Position');
        xapf = @(x,pos,xl) pos(3)*(x-min(xl))/diff(xl)+pos(1);          % 'x' Annotation Position Function
        yapf = @(y,pos,yl) pos(4)*(y-min(yl))/diff(yl)+pos(2);          % 'y' Annotation Position Function
        pos = gca().Position;
        xx=xlim; xl(1) = str2double(xx(1)); xl(2) = str2double(xx(2));
        yy=ylim; yl(1) = str2double(yy(1)); yl(2) = str2double(yy(2));
        annotation('ellipse', [xapf(x0-radiusx,pos,xl) yapf(y0-radiusy,pos,yl) radiusx*2*pos(3)/diff(xl) radiusy*2*pos(4)/diff(yl)],"Color","yellow","LineWidth",2, "LineStyle","--");
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if option_pause                                                             % wait for the user before moving on to the next file
    fprintf("\n Let me know when you are ready to move on to the next file... \n\n")
    pause();
end
close all;
end                                                                         % end the for loop across files
end

