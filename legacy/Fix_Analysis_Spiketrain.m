function Fix_Analysis_Spiketrain()

clc; close all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system = 3;                 % =1, Ubuntu; =2, Yang's PC; =3, Yang's MacBook
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

folderPath = '/path/to/project_folder';
group = 'groupA';

BURST_DETECT_METHOD = 1;                            % algorithm to be used for individual burst detection
                                                    %  =1       Maximal Interval 
NETWORK_BURST_DETECT_METHOD = 3;                    % algorithm to be used for newtork burst detection
                                                    %  =1       Population firing rate
                                                    %  =2       looking for peaks, combining the regular population firing rate and the normalized population firing; and using a progressive narrowing of the binSize
                                                    %  =3       option 2 + clustering
NETWORK_BURST_ANALYSIS_METHOD = 2;                  % algorithm to be used for network burst analysis
                                                    %  =1       single family of network bursting, with a single set of participating channels
                                                    %  =2       cluster analysis to identify different families of network bursts

                                                    
                                                    RATIO_MATCH = 0.60;                         % minimum matching between the treatment labels and the filename
TREATMENT{1} = 'CGP';
TREATMENT{2} = 'Gabazine';
TREATMENT{3} = 'DNQX';
TREATMENT{4} = 'D-AP5';
TREATMENT{5} = 'TTX';
n_treatment = size(TREATMENT,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% matlab file sprintf('%s\\%s\\%s_Analyzed.mat', main_folder, project_folder, file_data), containing the following variables
%  1) spikeTimestamp, 2) PerCondition, 3) NetworkBurstDetected
% PerCondition(cond) has the following fields (per each recording condition): 
%   1) spikeTrain, 2) SpikeFeature, 3) BurstFeature, 4) BurstDetected, 5) BurstFeature_Population, 6) TonicFeature_Population, 7) NetworkBurstFeature, 8) NetworkBurstStructure, 9) condition, 10) window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                        % factor to rescale from msec to sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify all the files

folder = fullfile(folderPath, group);
listing = dir(folder);

file_data = zeros(size(listing,1),1);
n_recordings=0;
for i=1:size(listing,1)
    if contains(listing(i).name,'Analyzed','IgnoreCase',true)               % look for all the files in the folder that contain "Analyzed" in their names
        file_data(i) = 1;                                                   % label the files that are going to be compiled
        n_recordings = n_recordings + 1;
        Recording{n_recordings,1} = extractBefore(listing(i).name,'_Analyzed');   % cell array with the name of all recording files found in the folder
    end
end

for i=1:n_recordings                                                        % across all experiments
    file_input = fullfile(folderPath, group, sprintf('%s_Analyzed.mat', Recording{i}));
    load(file_input, 'spikeTimestamp','PerCondition','n_conditions');

    file_manual = fullfile(folderPath, group, sprintf('%s_QCmanual.mat', Recording{i}));
    load(file_manual, 'rec_conditions');

    n_clusters = size(spikeTimestamp,2);
    dose = strings(1,n_conditions);
    for cond=1:n_conditions
        t_window = PerCondition(cond).window;
        spikeTrain = Reshape_spikeTimestamp(spikeTimestamp, t_window);          % reshape the structure spikeTimestamp for the specified time window
        NetworkBurstDetected = PerCondition(cond).NetworkBurstDetected;

        nBurst = nan(n_clusters,1);       BurstDetected = cell(n_clusters,1);   BurstFeature = cell(n_clusters,1);
    
        for j=1:n_clusters                                                          % across all clusters
            BurstDetected{j} = BurstDetection(spikeTrain(j).timeStamp, BURST_DETECT_METHOD);    % detect the bursts for cluster i
            BurstFeature{j} = BurstFeatures(BurstDetected{j});                      % characterize the bursts for cluster i
            nBurst(j) = BurstFeature{j}.nBurst;                                     % total number of burst events identified for cluster i
        end


        aux_cluster = split(rec_conditions{cond,1},".b");                               % assuming the raw data file was a .brw or .bxr file
        condition_file(cond) = aux_cluster(1);                                          % file name for the condition  
        pat = caseInsensitivePattern("_" + digitsPattern + ("u"|"m"|"n") + "M");
        if cond==1                                                              % assuming that the first condition is the baseline
        else
            aux_condition0 = rec_conditions{1,1};                               % assuming that the first condition is the baseline
            aux0 = split(aux_condition0,".b");
            condition0 = aux0(1);                                               % main name for the experiment
            aux_condition = extractAfter(rec_conditions{cond,1},condition0);
            aux_label = split(aux_condition,".b");                              % specific label for the condition, contained in the file name
            condition_label = [];
            for k=1:n_treatment
                if contains(convertCharsToStrings(aux_label{1}),lavelVAR(TREATMENT{k},RATIO_MATCH),'IgnoreCase',true)
                    condition_label = horzcat(condition_label,'+');
                    condition_label = horzcat(condition_label,TREATMENT{k});
                end
            end
            if size(condition_label)==0
            else
                aux = extractAfter(extract(aux_condition, pat),"_");
                if size(aux,1)>1
                    dose_aux = [];
                    for j=1:size(aux,1)
                        dose_aux = strcat(dose_aux, string(aux{j}));
                    end
                    dose(cond) = dose_aux; 
                else
                if size(aux,1)>0
                    dose(cond) = string(aux);
                end
            end
        end
        clear condition_label 
        end

        PerCondition(cond).dose = dose(cond);
        PerCondition(cond).file = rec_conditions{cond,1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyze the network bursting
 
nNetBurst1 = NetworkBurstDetected.nBurst;
nNetBurst2 = NetworkBurstDetected.Secondary.nBurst;
if ~isnan(nNetBurst1) && nNetBurst1 > 0

    switch NETWORK_BURST_ANALYSIS_METHOD
        
        case 1     % single family of network bursting, with a single set of participating channels
            netburst_member = ones(nNetBurst1,1);
            channels_part = NetworkBurstDetected.Channels(:,1);
            NetworkBurstFeature{1} = NetworkBurstFeatures(NetworkBurstDetected, spikeTrain, BurstDetected, t_window/K_MSEC_SEC, netburst_member, channels_part, 0);
            if nNetBurst2 > 0 
                netburst_member2 = ones(nNetBurst2,1);
                channels_part2 = NetworkBurstDetected.Secondary.Channels(:,1);
                NetworkBurstFeature{2} = NetworkBurstFeatures(NetworkBurstDetected.Secondary, spikeTrain, BurstDetected, t_window/K_MSEC_SEC, netburst_member2, channels_part2, 1);
            end

        case 2     % cluster analysis to identify different families of network bursts
            
            channels_part = NetworkBurstDetected.Channels(:,1);
            cluster1 = NetworkBurstDetected.clustering.label(1:nNetBurst1,1);
            aux_cluster = unique(cluster1);
            nFamilies = numel(aux_cluster);
            for k=1:nFamilies
                netburst_member_k = (cluster1 == aux_cluster(k));
                NetworkBurstFeature{k} = NetworkBurstFeatures(NetworkBurstDetected, spikeTrain, BurstDetected, t_window/K_MSEC_SEC, netburst_member_k, channels_part, nNetBurst1, 0);
            end

            if nNetBurst2 > 0 
                netburst_member2 = ones(nNetBurst2,1);
                channels_part2 = NetworkBurstDetected.Secondary.Channels(:,1);
                NetworkBurstFeature{k+1} = NetworkBurstFeatures(NetworkBurstDetected.Secondary, spikeTrain, BurstDetected, t_window/K_MSEC_SEC, netburst_member2, channels_part2, nNetBurst1, 1);
            end
    end
else
   netburst_member = [];
   channels_part = ones(n_clusters,1);
   NetworkBurstFeature{1} = NetworkBurstFeatures(NetworkBurstDetected, spikeTrain, BurstDetected, t_window/K_MSEC_SEC, netburst_member,channels_part, nNetBurst1, 0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare the output

    PerCondition(cond).NetworkBurstFeature = NetworkBurstFeature;

    clear NetworkBurstFeature 
    
    end
    
    save(file_input, 'PerCondition', '-append');
    clear('spikeTimestamp','PerCondition','n_conditions');
end

end