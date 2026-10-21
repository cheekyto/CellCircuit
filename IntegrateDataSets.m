function IntegrateDataSets(inputs)

% function to integrate the dataset coming from Kilosor spike sorting and 3Brain spike sorting.
% The datasets has been "shaped" and gone through "quality control"
% (optional) inputs:    structure to overwrite any variable of the sections INPUTS and PARAMETERS (e.g., inputs.folder, inputs.file_data_Kilosort); see "override_inputs"
%                       when it is not provided, the values defined in the code are used
% (revised October 2026: optional inputs; cell arrays n x 1 instead of n x n; no "clear" inside the loops; Kilosort + 3Brain no longer crashes;
%  3Brain units with bad waveform are always discarded; n_cluster_Kilosort / n_cluster_3Brain are the final numbers of clusters)

if ~parameter_set('active'), clc; end; close all;      % (revised October 2026) no "clc" while "run_pipeline" runs (it would erase its progress messages)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system = 3;                 % =1, Ubuntu; =2, Yang's PC; =3, Yang's MacBook; =0, do not modify the MATLAB path (e.g., test harness)
if nargin < 1, inputs = struct(); end                                       % no inputs provided -> use the values defined in the code
if isfield(inputs, 'system'), system = inputs.system; end                   % the system can also be provided as an input
switch system
    case 1      % Anlysis PCP - Ubuntu
        addpath(genpath('/home/cortex/Documents/Yang'));
    case 2      % Yang's PC - GitHub
        addpath(genpath('C:\Users\Yang\OneDrive - Harvard University\Documents\GitHub\CellCircuit\'));
    case 3      % Yang's MacBook
        addpath(genpath('/Users/sungminyang/Documents/GitHub/CellCircuit'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = '/path/to/recording_folder';
file_data_Kilosort = 'recording_name';
file_data_3Brain = []; 
file_output = 'recording_name';      % if you select a single dataset, it will keep the original filename

include_Kilosort = 1;
unitSelect_Kilosort = [1; 2; 3];                % classification of the unit type, performed by "bombcell", that will be considered for downstream analysis
                                                % =0, noise; =1, good; =2, MUA; =3, non-somatic
include_3Brain = 0;
unitSelect_3Brain = [1];                  % waveform-based classification of the unit type, that will be considered for downstream analysis
                                                % =0, unknown; =1, good; =2, weak; =-1, noise
option_criteria = 3;                        % =1, only spatial superposition
                                            % =2, only spike train similarity
                                            % =3, spatial superposition or spike train similarity
                                            % =4, spatial superposition and spike train similarity    
option_plot = 0;
option_print_list = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% matlab file sprintf('%s\\%s\\%s_integrated_QCmanual.mat', main_folder, project_folder, file_data), containing the following variables
%  1) spikeTimestamp, 2) PerCondition, 3) NetworkBurstDetected

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_SPIKES_MINIMUM = 20;                       % minimum number of hypersynchronized spikes to trigger assessment; it will apply to any category (bursting, nonbursting and total) (default: 20)

N_SPIKES_BURST_CUT = 50;                     % threshold for the percentage of hypersynchonized spikes within burst, in % (default: 50%)
N_SPIKES_NONBURST_CUT = 20;                  % threshold for the percentage of hypersynchonized spikes outside burst to make the first cut and proceed with the shuffling, in % (default: 20%)
N_SPIKES_TOTAL_CUT = 40;                     % threshold for the percentage of hypersynchonized spikes to make the first cut and proceed with the shuffling, in % (default: 40%)
N_SPIKES_3BRAIN_CUT = 50;                    % threshold for the percentage of hypersynchonized spikes in the 3Brain signal, in % (default: 60%)

K_MSEC_SEC = 1/1000;                        % factor to rescale from msec to sec
OPTION_METHOD_BURST = 1;                    % method for the bursting detection (=1, Maximal Interval)
REFRACTORY_PERIOD = 0;                      % refractory period for the shuffling, in msec (default: 0)
METHOD_SHUFFLING = 2;                       % method for the spike-train shuffling (=1, uniform; =2, biphasic <burst and non-bust>) 

apply_parameters('IntegrateDataSets');                                      % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")
override_inputs(inputs);                    % overwrite the INPUTS / PARAMETERS above with the ones provided as input (if any)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the datasets

file_Kilosort = fullfile(folder, sprintf('%s_QCmanual.mat',file_data_Kilosort));
file_3Brain = fullfile(folder, sprintf('%s_QCmanual.mat',file_data_3Brain));

if include_Kilosort
%    load(file_Kilosort,'spikeTimestamp','rec_conditions','n_clusterOriginal');  
    load(file_Kilosort,'spikeTimestamp','rec_conditions','good_cluster');  
    spikeTimestamp_Kilosort = spikeTimestamp;
    rec_conditions_Kilosort = rec_conditions;
%    n_clusterOriginal_Kilosort = n_clusterOriginal;
    n_clusterOriginal_Kilosort = good_cluster;
    n_cluster_Kilosort = size(spikeTimestamp_Kilosort,2);
    clear('spikeTimestamp','rec_conditions','n_clusterOriginal');

    metadata = readtable(fullfile(folder, 'metadata.csv'));
    unitType_list_Kilosort = metadata.unitType;                            % type of unit; classification made by "bombcell"
    clusterID_list_Kilosort = metadata.new_cluster_id;
else
    n_cluster_Kilosort = 0;
    bad_waveform_Kilosort = nan;
    similarCluster = nan; 
    overlapingCluster = nan;
end
if include_3Brain
    load(file_3Brain,'spikeTimestamp','rec_conditions','n_clusterOriginal','unitLabel');  
    spikeTimestamp_3Brain = spikeTimestamp;
    rec_conditions_3Brain = rec_conditions;
    n_clusterOriginal_3Brain = n_clusterOriginal;
    n_cluster_3Brain = size(spikeTimestamp_3Brain,2);
    clear('spikeTimestamp','rec_conditions','n_clusterOriginal');  
else
    n_cluster_3Brain = 0;
    bad_waveform_3Brain = nan;
    similarCluster = nan;
    overlapingCluster = nan;
    discard_3Brain = nan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare the recording conditions

if include_Kilosort
    fprintf("\n\n The recording conditions for Kilosort dataset are:")
    for i=1:size(rec_conditions_Kilosort,1)
        fprintf("\n %s", rec_conditions_Kilosort{i,1});
    end
end
if include_3Brain
    fprintf("\n\n The recording conditions for 3Brain dataset are:")
    for i=1:size(rec_conditions_3Brain,1)
        fprintf("\n %s", rec_conditions_3Brain{i,1});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select the desired units

bad_waveform_Kilosort = [];
if include_Kilosort
    for i=1:n_cluster_Kilosort
        cluster_id = spikeTimestamp_Kilosort(i).cluster;
        idx = clusterID_list_Kilosort == cluster_id;
        if ~ismember(unitType_list_Kilosort(idx), unitSelect_Kilosort)
            bad_waveform_Kilosort = vertcat(i, bad_waveform_Kilosort);
        end
    end
end

bad_waveform_3Brain = [];
if include_3Brain
    for i=1:n_cluster_3Brain
        if ~ismember(unitLabel(i, 2),unitSelect_3Brain)
            bad_waveform_3Brain = vertcat(i, bad_waveform_3Brain);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare the two datasets (position in the electrode grid)

overlapingCluster = [];
if include_Kilosort && include_3Brain
    if n_cluster_Kilosort>0
        for i=1:n_cluster_Kilosort
            if ~ismember(i, bad_waveform_Kilosort)
                cluster_idx = spikeTimestamp_Kilosort(i).position;   
                for k=1:n_cluster_3Brain 
                    if ~ismember(k, bad_waveform_3Brain) && cluster_idx == spikeTimestamp_3Brain(k).cluster
                        overlapingCluster = vertcat([i cluster_idx k spikeTimestamp_3Brain(k).cluster], overlapingCluster);   % [idx_kilosort cluster_kilosort idx_3BRain cluster_3Brain ]
                    end
                end
            end
        end
    end
end

      if option_print_list
         if size(overlapingCluster,1)>0
            fprintf("\n \n The overlaping clusters are: \n")
            for k=1:size(overlapingCluster,1)
                ik = overlapingCluster(k,1);
                i3 = overlapingCluster(k,3);
                fprintf("Kilosort: element %g, cluster %g (x: %g; y: %g) -- 3Brain: element %g, cluster %g (x: %g; y: %g) \n", ...
                 ik, overlapingCluster(k,2),spikeTimestamp_Kilosort(ik).positionXY(1),spikeTimestamp_Kilosort(ik).positionXY(2), ...
                    i3, overlapingCluster(k,4),spikeTimestamp_3Brain(i3).positionXY(1),spikeTimestamp_3Brain(i3).positionXY(2));
            end
         end
      end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assess the similarity between the two spike trains

if include_Kilosort && include_3Brain

t_step = spikeTimestamp_Kilosort(1).t_step;                                                        % temporal step, in msec
range = [spikeTimestamp_Kilosort(1).window(1) spikeTimestamp_Kilosort(1).window(2)]./K_MSEC_SEC;  % temporal window for analysis, in msec (entire recording)

similarCluster = [];
BurstDetected_Kilosort = cell(n_cluster_Kilosort,1); BurstFeature_Kilosort = cell(n_cluster_Kilosort,1);                                          % initialize the variables for burst structure (revised October 2026: n x 1, it was n x n)
for i=1:n_cluster_Kilosort
    if ~ismember(i, bad_waveform_Kilosort)
        BurstDetected_Kilosort{i} = BurstDetection(spikeTimestamp_Kilosort(i).timeStamp, OPTION_METHOD_BURST);  % Burst detection using spike train of cluster i  
        BurstFeature_Kilosort{i} = BurstFeatures(BurstDetected_Kilosort{i});                                              % Burst characterization for cluster i
    end
end

for k=1:n_cluster_3Brain
    
    if ~ismember(k, bad_waveform_3Brain) 
    
        BurstDetected_3Brain = BurstDetection(spikeTimestamp_3Brain(k).timeStamp, OPTION_METHOD_BURST);  % Burst detection using spike train of cluster ik  
       % BurstFeature_3Brain = BurstFeatures(BurstDetected_3Brain);                                        % Burst characterization for cluster ik
        signal1 = BurstDetected_3Brain.spikes(:,[1 3]);                             % time stamp for cluster k (column 1) and index of the burst where the spike belongs (column 2; =0 if outside of burst) 

        for i=1:n_cluster_Kilosort
            if ~ismember(i, bad_waveform_Kilosort)
                signal2 = BurstDetected_Kilosort{i}.spikes(:,[1 3]);                     % time stamp for cluster i (column 1) and index of the burst where the spike belongs (column 2; =0 if outside of burst)
                QC = HyperSynTrain(signal1,signal2,t_step,range,0,METHOD_SHUFFLING,REFRACTORY_PERIOD,BurstDetected_Kilosort{i},BurstFeature_Kilosort{i});  % analyze the hypersynchronicity for cluster i vs j (no iterations) 
         
                cluster_idx = spikeTimestamp_Kilosort(i).position; 
                cond_burst = QC.burst.nspikes_hs*100 > N_SPIKES_BURST_CUT && QC.burst.nspikes_total*QC.burst.nspikes_hs > N_SPIKES_MINIMUM;
                cond_nonburst = QC.nonburst.nspikes_hs*100 > N_SPIKES_NONBURST_CUT && QC.nonburst.nspikes_total*QC.nonburst.nspikes_hs > N_SPIKES_MINIMUM;
                cond_total = QC.total.nspikes_hs*100 > N_SPIKES_TOTAL_CUT && QC.total.nspikes_total*QC.total.nspikes_hs > N_SPIKES_MINIMUM;
                cond_3Brain = QC.total.ch_1.nspikes_hs*100 > N_SPIKES_3BRAIN_CUT;
                if cond_burst || cond_nonburst || cond_total || cond_3Brain
                   similarCluster = vertcat([i cluster_idx k spikeTimestamp_3Brain(k).cluster], similarCluster);   % [idx_kilosort cluster_kilosort idx_3BRain cluster_3Brain ]                    
                end
            end
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discard the overlaping clusters and integrate the datasets

discard_Kilosort = [];
discard_3Brain = [];

if include_Kilosort 
    if include_3Brain           % Kilosort + 3Brain
        discard_Kilosort = bad_waveform_Kilosort;
        discard_3Brain = sort(bad_waveform_3Brain);                         % (revised October 2026) the 3Brain units with bad waveform are always discarded
                                                                            % (before, they were kept when no overlaping/similar clusters were found)
        switch option_criteria
            case 1      % only spatial superposition
                if numel(overlapingCluster)>0
                    discard_3Brain = union(overlapingCluster(:,3), bad_waveform_3Brain);
                end
            case 2      % only spike train similarity
                if numel(similarCluster)>0
                    discard_3Brain = union(similarCluster(:,3), bad_waveform_3Brain);
                end
            case 3      % spatial superposition or spike train similarity
                if numel(overlapingCluster)>0
                    aux = union(overlapingCluster(:,3), bad_waveform_3Brain);
                    if numel(similarCluster)>0
                        discard_3Brain = union(aux,similarCluster(:,3));
                    else
                        discard_3Brain = aux;
                    end
                elseif numel(similarCluster)>0
                    discard_3Brain = union(bad_waveform_3Brain, similarCluster(:,3));
                end

            case 4      % spatial superposition and spike train similarity
                if numel(similarCluster)>0 && numel(overlapingCluster)>0
                    discard_3Brain = union(intersect(overlapingCluster(:,3),similarCluster(:,3)), bad_waveform_3Brain);
                end
        end
    else                        % Kilosort only   
        discard_Kilosort = bad_waveform_Kilosort;
    end
    rec_conditions = rec_conditions_Kilosort;

elseif include_3Brain           % 3Brain only
      discard_3Brain = bad_waveform_3Brain;
      rec_conditions = rec_conditions_3Brain;
else                      % Neither Kilosort or 3Brain
      error('Please select at least one of the dataset')
end

if include_Kilosort
    for i=n_cluster_Kilosort:-1:1
        if ismember(i,discard_Kilosort)
            spikeTimestamp_Kilosort(i) = [];
        else
            spikeTimestamp_Kilosort(i).idx_toc = nan;
        end
    end
    spikeTimestamp = spikeTimestamp_Kilosort;
end

if include_3Brain
    spikeTimestamp_3Brain(discard_3Brain) = [];                             % eliminate the discarded 3Brain clusters (all at once)
    if include_Kilosort                                                     % (revised October 2026) the Kilosort clusters were already discarded above; before, they were
       spikeTimestamp = cat(2,spikeTimestamp_Kilosort, spikeTimestamp_3Brain);  % deleted a second time here, with indexes that had already shifted (crash, or wrong clusters deleted)
    else
       spikeTimestamp = spikeTimestamp_3Brain;
    end
end

        % (revised October 2026) number of clusters that remain from each dataset after the integration (it is reported by "Compile_MEAproject"
        % as the number of good clusters from Kilosort / 3Brain); before, it was the number of clusters before discarding any of them
n_cluster_Kilosort_original = n_cluster_Kilosort;                           % number of clusters in each dataset before the integration
n_cluster_3Brain_original = n_cluster_3Brain;
if include_Kilosort,    n_cluster_Kilosort = size(spikeTimestamp_Kilosort,2);    end
if include_3Brain,      n_cluster_3Brain = size(spikeTimestamp_3Brain,2);        end

if include_Kilosort
    if include_3Brain
        n_clusterOriginal = nan;
        file_save = fullfile(folder,sprintf('%s_integrated_QCmanual.mat', file_output));
    else
        file_output = file_data_Kilosort;
        n_clusterOriginal = n_clusterOriginal_Kilosort;
        file_save = fullfile(folder,sprintf('%s_curated_QCmanual.mat', file_output));
    end
else
    if include_3Brain
        file_output = file_data_3Brain;
        n_clusterOriginal = n_clusterOriginal_3Brain;
        file_save = fullfile(folder,sprintf('%s_curated_QCmanual.mat', file_output));
    else
        error('\n Not known dataset was selected!!! \n');
    end    
end
save(file_save, 'spikeTimestamp','rec_conditions','n_cluster_Kilosort','n_cluster_3Brain','bad_waveform_Kilosort','bad_waveform_3Brain','similarCluster','overlapingCluster','discard_3Brain','n_clusterOriginal','n_cluster_Kilosort_original','n_cluster_3Brain_original');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display a summary of the results (and parameters)

if include_Kilosort
    fprintf("\n\n You have included a dataset sorted by Kilosort.");
    fprintf("\n Based on their waveform, you have included units that are: ");
    for i=1:numel(unitSelect_Kilosort)
        switch unitSelect_Kilosort(i)
            case 0
                fprintf(" - noise");
            case 1
                fprintf(" - good");
            case 2
                fprintf(" - MUA");
            case 3
                fprintf(" - non-somatic");
        end
    end 
    n_Kilosort = size(spikeTimestamp_Kilosort,2);
else
    n_Kilosort = 0;
end
if include_3Brain
    fprintf("\n\n You have included a dataset sorted by 3Brain.");
    fprintf("\n Based on their waveform, you have included units that are: ");
    for i=1:numel(unitSelect_3Brain)
        switch unitSelect_3Brain(i)
            case 0
                fprintf(" - unknown");
            case 1
                fprintf(" - good");
            case 2
                fprintf(" - weak");
            case -1
                fprintf(" - noise");
        end
    end  
    n_3Brain = size(spikeTimestamp_3Brain,2);
else
    n_3Brain = 0;
end

fprintf("\n\n Final number of clusters after integration: %g units", size(spikeTimestamp,2));
fprintf("\n %g units are from the Kilosort dataset & %g units are from the 3Brain dataset", n_Kilosort, n_3Brain);

if include_Kilosort
   fprintf("\n\n Kilosort dataset originally had: %g units - ", n_cluster_Kilosort_original);
   fprintf(" %g units were discarded because of their waveform", numel(bad_waveform_Kilosort));
 
   if include_3Brain
      fprintf("\n\n 3Brain dataset originally had: %g units - ", n_cluster_3Brain_original);
      fprintf(" %g units were discarded because of their waveform", numel(bad_waveform_3Brain));

      fprintf("\n\n The criteria to assess the overlapping between the two datasets is: ");
      switch option_criteria
          case 1
              criteria = 'Only spatial superposition (same electrode)';
          case 2
              criteria = 'Only spike train similarity (overlaping spike time stamps)';
          case 3
              criteria = 'Spatial superposition (same electrode) OR spike train similarity (overlaping spike time stamps)';
          case 4
              criteria = 'Spatial superposition (same electrode) AND spike train similarity (overlaping spike time stamps)';
       end
      fprintf(" %s \n", criteria);
      fprintf(" From that assessment, %g units from the 3Brain dataaset were discarded because of their overlaping with Kilosort dataset (after removing bad waveform) \n", numel(discard_3Brain)-numel(bad_waveform_3Brain));
      switch option_criteria
          case 3
              fprintf("  (%g cluster pairs were found with spatial superposition and %g with similar spike trains)", size(overlapingCluster,1),size(similarCluster,1))
          case 4
              fprintf("  (%g clusters pairs were found with spatial superposition and %g with similar spike trains)", size(overlapingCluster,1),size(similarCluster,1))
      end
   end
end

fprintf("\n\n The integrated has been saved \n");

%%% save the final data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot

if option_plot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure()        % raster plot of integrated units

    for i=1:size(spikeTimestamp,2)
        data = spikeTimestamp(i);
        if isnan(data.idx_toc)         % from Kilosort dataset
           plot(data.timeStamp(:)*K_MSEC_SEC, ones(data.nUnits,1)+i,'|',"MarkerEdgeColor",'red'); hold on
        else                                        % from 3Brain dataset
           plot(data.timeStamp(:)*K_MSEC_SEC, ones(data.nUnits,1)+i,'|',"MarkerEdgeColor",'blue'); hold on
        end
    end
    ylim([0 i+1]);        xlim([data.window(1) data.window(2)]);             xlabel('Time (sec)');            ylabel('Cluster');
    str = {'Intergrated Dataset - GOOD cllusters (\color{red}Kilosort - \color{blue}3Brain)'};
    title(str)

end
end