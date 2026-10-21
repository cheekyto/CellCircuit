function Compare_Kilosort_vs_3Brain()

clc; close all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system = 3;                 % =1, Ubuntu; =2, Yang's PC; =3, Yang's MacBook
switch system
    case 1      % Anlysis PCP - Ubuntu
        addpath(genpath('/home/cortex/Documents/Yang'));
        delimiter = '/';
    case 2      % Yang's PC - GitHub
        addpath(genpath('C:\Users\Yang\OneDrive - Harvard University\Documents\GitHub\MEA-analysis\'));
        delimiter = '\\';
    case 3      % Yang's MacBook
        addpath(genpath('/Users/sungminyang/Documents/GitHub/CellCircuit'));
        delimiter = '/';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folderPath = '/path/to/recording_folder';             % Folder where to find the original data (xxx.bxr files)
file_data_3Brain = 'recording_name_3Brain';
file_data_Kilosort = 'recording_name_AllUnits';
file_output = 'recording_name kilosort vs 3Brain';  % (revised October 2026) file in folderPath where options 4 and 5 save their results ("unitLabel", "units_3Brain", 
                                                % "units_Kilosort"), as in "Brain3_processedData"; before, option 4 appended them to "xxx_3Brain_raw.mat", which is 
                                                % overwritten whenever "ShapeSpikeTrainData_3Brain" runs again, and option 5 read them from "xxx.mat", which no code creates

option_plot = 2;                            % =1, raster plot of all clusters   
                                            % =2, waveform and raster plot of the most active electrodes, one-by-one
                                            % =3, compare result from Kilosort vs 3Brain (using kilosort as reference)    
                                            % =4, compare result from Kilosort vs 3Brain (using 3Brain as reference)
                                            % =5, summarize the comparison Kilosort vs 3Brain        

QC_3Brain = 1;                              % =1, the 3Brain data has gone through QC (there is a "xxx_QCmanual.mat" file)
                                            % =0, the 3Brain data is raw, no QC has been performed (there is only a "xxx_raw.mat" file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                        % factor to rescale from msec to sec
MATRIX_SIZE = 64;
NEW_CLUSTER_ID = 4;
UNIT_TYPE = 2;
PEAK_CHANNEL = 6;
WAVEFORM_TEMPLATE = 5;
N_CLUSTER_PLOT = 10;
N_NEIGHBOR = 1;                             % to determine the matrix of neighbors to the electrode of interest (move as much as N_NEIGHBOR electrodes to each side)
BIN = 1000;                                 % bin size in msec
PERCENTILE = 0.975;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data

if QC_3Brain
    file_3Brain1 = sprintf('%s%s%s_3Brain_raw.mat', folderPath, delimiter, file_data_3Brain);
    file_3Brain2 = sprintf('%s%s%s_3Brain_QCmanual.mat', folderPath, delimiter, file_data_3Brain);
    load(file_3Brain1, 'Indiv_Cond','waveformLengths','samplingRates','dacFactors','offsetValues');
    load(file_3Brain2, 'spikeTimestamp','rec_conditions','good_cluster');
else
    file_3Brain = sprintf('%s%s%s_3Brain_raw.mat', folderPath, delimiter, file_data_3Brain);
    load(file_3Brain, 'spikeTimestamp','rec_conditions', 'Indiv_Cond','waveformLengths','samplingRates','dacFactors','offsetValues');
end
spikeTimestamp_3Brain = spikeTimestamp;
Indiv_Cond_3Brain = Indiv_Cond;
n_clusters_3Brain = size(spikeTimestamp,2);
waveformLength = waveformLengths(1);
samplingRate = samplingRates(1);
clear spikeTimestamp Indiv_Cond

file_Kilosort = sprintf('%s%s%s_raw.mat', folderPath, delimiter, file_data_Kilosort);
load(file_Kilosort, 'spikeTimestamp');
spikeTimestamp_Kilosort = spikeTimestamp;
file_name = sprintf('%s%swaveforms.mat',folderPath,delimiter);
load(file_name, 'waveforms'); 
waveforms_Kilosort = waveforms;

clear spikeTimestamp file_name waveforms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the raster dataset with spike times and convert spike times from frame position to seconds for each channel

cont_spikes = zeros(n_clusters_3Brain,1);
for i=1:n_clusters_3Brain
    cont_spikes(i) = spikeTimestamp_3Brain(i).nUnits;
end

if QC_3Brain
    threshold_activity = 0;
else
    threshold_activity = quantile(cont_spikes,PERCENTILE);
end
activeClusters = find(cont_spikes>threshold_activity);

map = nan(MATRIX_SIZE,MATRIX_SIZE); 
for idx=1:numel(activeClusters)    
    chIdx = activeClusters(idx)-1;                                          % identifier of channel to be analysed
    pos_x = spikeTimestamp_3Brain(chIdx+1).positionXY(1)+1;
    pos_y = spikeTimestamp_3Brain(chIdx+1).positionXY(2)+1;
    map(pos_y,pos_x)=-1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the population firing rate data

spike_data = spikeTimestamp_3Brain(activeClusters(1)).timeStamp;
tRange = spikeTimestamp_3Brain(1).window;
aux_bin = bin_spike_time(spike_data, BIN, tRange);

spike_bin_indiv = nan(numel(aux_bin.count), numel(activeClusters));
spikes_bin = nan(numel(aux_bin.count),2);
spikes_bin(:,2) = aux_bin.binMidPoint(:);

for idx=1:numel(activeClusters)    
    spike_data = spikeTimestamp_3Brain(activeClusters(idx)).timeStamp;
    aux = bin_spike_time(spike_data, BIN, tRange);
    spike_bin_indiv(:,idx) = aux.count(:);
    
    clear spike_data aux
end
spikes_bin(:,1) = mean(spike_bin_indiv,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting option

switch option_plot                      
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1        % raster plot of all clusters 
        
        figure(1)
        for k=1:n_clusters_3Brain
            data = spikeTimestamp_3Brain(k);
            plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+k,'|',"MarkerEdgeColor",[0 0 0]);      hold on;      clear data;
        end 
        ylim([0 1+k]);        xlabel('Time (sec)');        ylabel('Neuron');    title('ALL electrodes - 3Brain')

        figure(2)
        histogram(cont_spikes)

        figure(3)
        for k=1:numel(activeClusters)
            data = spikeTimestamp_3Brain(activeClusters(k));
            plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+k,'|',"MarkerEdgeColor",[0 0 0]);      hold on;        clear data;
        end 
        ylim([0 1+k]);        xlabel('Time (sec)');        ylabel('Neuron');    title('Most active electrodes - 3Brain')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2       % waveform and raster plot of the most active electrodes, one-by-one

        x = double((1:waveformLength)) / samplingRate;
        contador=0;
        for idx=1:numel(activeClusters)
            chIdx = activeClusters(idx);                                    % identifier of channel to be analysed
            fprintf(sprintf('\n Electrode number   %g', spikeTimestamp_3Brain(chIdx).cluster));
            fprintf(sprintf('\n It has %g units \n', spikeTimestamp_3Brain(chIdx).nUnits));
            
            cond = 1;
            waveformData = Indiv_Cond_3Brain{cond,2};
            dacFactor = dacFactors(cond);
            offsetValue = offsetValues(cond);

            traces = waveformData{spikeTimestamp_3Brain(chIdx).cluster+1};
            y = offsetValue + dacFactor * mean(traces,1);
            dy = std(offsetValue + dacFactor*traces,0,1);

            fig1=figure();
            subplot(2,2,[1,2])
            for k=1:numel(activeClusters)
                data = spikeTimestamp_3Brain(activeClusters(k));
                plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+k,'|',"MarkerEdgeColor",[0 0 0]);
                hold on;                clear data
            end 
            data = spikeTimestamp_3Brain(chIdx);
            plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+idx,'|',"MarkerEdgeColor","r");
            clear data;            ylim([0 1+k]);            xlabel('Time (sec)');            ylabel('Neuron')

            subplot(2,2,3)
            for i=1:size(traces ,1)
                yy = offsetValue + dacFactor * traces(i,:);
                plot(x,yy);                hold on
            end

            subplot(2,2,4)
            fill([x';flipud(x')],[(y-dy)';flipud((y+dy)')],[0 1 1],'linestyle','none')
            line(x,y);

            prompt = "Save the waveform? (y=yes / n=no):  ";
            answer = input(prompt,"s"); 
            switch answer 
                case 'y'
                    contador = contador + 1;
                    file_wave = sprintf('%s%s%s_3Brain_waveform%g.fig', folderPath, delimiter, file_data_3Brain, contador);
                    savefig(fig1,file_wave);
                otherwise
            
            end

            fig2 = figure(2);
            pos_x = spikeTimestamp_3Brain(chIdx).positionXY(1)+1;
            pos_y = spikeTimestamp_3Brain(chIdx).positionXY(2)+1;
            map(pos_y,pos_x)=1;
            heatmap(map,'XDisplayLabels',NaN*ones(size(map,1),1),'YDisplayLabels',NaN*ones(size(map,1),1),'Colormap',cool);
            
            pause()
            map(pos_y,pos_x)=0;     close(fig1);            close(fig2);
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3              % compare result from Kilosort vs 3Brain (using kilosort as reference)     

        file_excel = sprintf('%s%smetadata.csv', folderPath, delimiter);
        T = readtable(file_excel);         
        for i=1:size(spikeTimestamp_Kilosort,2)
            cluster_idx = spikeTimestamp_Kilosort(i).position;                              
            extra_cluster = findClusterID(spikeTimestamp_Kilosort, cluster_idx, N_NEIGHBOR, 2);
            idx_3Brain_neighbors = findClusterID(spikeTimestamp_3Brain, cluster_idx, N_NEIGHBOR, 3);
            idx_3Brain  = idx_3Brain_neighbors(ceil(numel(idx_3Brain_neighbors)/2));

            if isnan(idx_3Brain)
                fprintf('\n The cluster could NOT be found in the 3Brain dataset \n')
            else
            
                % plot raster plot (3Brain vs Kilosort)
                fig1 = figure(1);
                n_clusters_plot = min(size(spikeTimestamp_Kilosort,2),N_CLUSTER_PLOT);
                for k=1:n_clusters_plot
                    data = spikeTimestamp_Kilosort(k);
                    plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+k,'|',"MarkerEdgeColor",[0 0 0]);       hold on;                clear data
                end 
                shift_x0 = k+n_clusters_plot/2;
                data = spikeTimestamp_Kilosort(i);
                plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+shift_x0,'|',"MarkerEdgeColor",'r');            hold on
                shift_x = shift_x0 + 1;
                plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+shift_x,'|',"MarkerEdgeColor",'r');            hold on

                shift_x = shift_x + 1;
                if numel(extra_cluster)>1
                    for k=1:numel(extra_cluster)
                        if extra_cluster(k) ~= i
                            data = spikeTimestamp_Kilosort(extra_cluster(k));
                            plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+shift_x,'|');                        hold on
                            plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+shift_x0,'|',"MarkerEdgeColor",'m');                        hold on
                            shift_x = shift_x + 1;
                        end
                    end
                end

                data_3Brain = spikeTimestamp_3Brain(idx_3Brain);
                plot(data_3Brain.timeStamp*K_MSEC_SEC, ones(data_3Brain.nUnits,1)+shift_x,'|',"MarkerEdgeColor",'b');            hold on
                plot(data_3Brain.timeStamp*K_MSEC_SEC, ones(data_3Brain.nUnits,1)+shift_x0,'|',"MarkerEdgeColor",'b');            ylim([0 shift_x+4]);            xlabel('Time (sec)');            ylabel('Cluster');
                str = {'\color{black}Kilosort - \color{red}Kilosort - \color{magenta}Kilosort (extra) - \color{blue}3Brain'};       title(str)

                % plot raster plot (3Brain; selected vs neighboring electrodes)
                fig2 = figure(2);
                data = spikeTimestamp_Kilosort(i);
                plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1),'|',"MarkerEdgeColor",'r');            hold on

                shift_x = 2;
                data_3Brain = spikeTimestamp_3Brain(idx_3Brain);
                plot(data_3Brain.timeStamp*K_MSEC_SEC, ones(data_3Brain.nUnits,1)+shift_x,'|',"MarkerEdgeColor",'b');            hold on
                shift_x = shift_x + 1;
                for j=1:power(N_NEIGHBOR*2+1,2)
                    if idx_3Brain_neighbors(j)<power(MATRIX_SIZE,2)+1 && idx_3Brain_neighbors(j)>0
                        data_3Brain = spikeTimestamp_3Brain(idx_3Brain_neighbors(j));
                    else
                        data_3Brain.timeStamp = 0;
                        data_3Brain.nUnits = 1;
                    end
                    if j==ceil(power(N_NEIGHBOR*2+1,2)/2)
                        plot(data_3Brain.timeStamp*K_MSEC_SEC, ones(data_3Brain.nUnits,1)+shift_x,'|',"MarkerEdgeColor",'b');
                    else
                        plot(data_3Brain.timeStamp*K_MSEC_SEC, ones(data_3Brain.nUnits,1)+shift_x,'|',"MarkerEdgeColor",[0 0 0]);
                        plot(data_3Brain.timeStamp*K_MSEC_SEC, ones(data_3Brain.nUnits,1)+2,'|',"MarkerEdgeColor",'m');
                    end
                    shift_x = shift_x + 1;                hold on
                end
                ylim([0 shift_x+1]);            xlabel('Time (sec)');            ylabel('Cluster');
                str = {'\color{black}3Brain - \color{magenta}3Brain (neighbors) - \color{blue}3Brain (target) - \color{red}Kilosort'};      title(str)
        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                idx_waveform = 0;            k=1;
                while idx_waveform<1
                    if T{k,NEW_CLUSTER_ID}==spikeTimestamp_Kilosort(i).cluster
                        idx_waveform = k;
                    else
                        k=k+1;
                    end
                end
                if spikeTimestamp_Kilosort(i).position ~= T{idx_waveform,PEAK_CHANNEL}
                    error('Error finding the cluster - the electrode position does not match')
                end

                % plot waveforms (3Brain vs Kilosort)
                fig3 = figure(3);        
                subplot(1,2,1)
                score_aux = T{idx_waveform,UNIT_TYPE};
                switch score_aux
                    case 0
                        score_unit = 'noise';
                    case 1
                        score_unit = 'good';
                    case 2
                        score_unit = 'MUA';
                    case 3
                        score_unit = 'non-somatic';
                end
                xx = double((1:numel(waveforms_Kilosort{idx_waveform}))') / samplingRate; 
                plot(xx,waveforms_Kilosort{idx_waveform}, 'r');         title(sprintf('Kilosort (%s)', score_unit));        
                if numel(extra_cluster)>1
                    for k=1:numel(extra_cluster)
                        if extra_cluster(k) ~= i
                            idx_waveform_extra = 0;
                            j=1;
                            while idx_waveform_extra<1
                                if T{j,NEW_CLUSTER_ID}==spikeTimestamp_Kilosort(extra_cluster(k)).cluster
                                    idx_waveform_extra = j;
                                else
                                    j=j+1;
                                end
                            end
                            hold on; plot(xx,waveforms_Kilosort{idx_waveform_extra}, 'r--');
                        end
                    end
                end

                subplot(1,2,2)
                cond = 1;
                waveformData = Indiv_Cond_3Brain{cond,2};
                dacFactor = dacFactors(cond);
                offsetValue = offsetValues(cond);
                x = double((1:waveformLength)) / samplingRate;
                traces = waveformData{spikeTimestamp_3Brain(idx_3Brain).cluster+1};
                y = offsetValue + dacFactor * mean(traces,1);
                dy = std(offsetValue + dacFactor*traces,0,1);
                fill([x';flipud(x')],[(y-dy)';flipud((y+dy)')],[0 1 1],'linestyle','none')
                line(x,y);          title('3Brain')

                % plot waveforms (3Brain; selected vs neighboring electrodes)
                fig4 = figure(4); 
                for j=1:power(N_NEIGHBOR*2+1,2)
                    subplot(3,3,j)
                    idx_ch = idx_3Brain_neighbors(j);
                    if ~isnan(idx_ch) && idx_ch>0
                        traces = waveformData{spikeTimestamp_3Brain(idx_ch).cluster+1};
                        if numel(traces)>0
                            y = offsetValue + dacFactor * mean(traces,1);
                            dy = std(offsetValue + dacFactor*traces,0,1);
                            fill([x';flipud(x')],[(y-dy)';flipud((y+dy)')],[0 1 1],'linestyle','none')    
                        else
                            y = x*0;
                        end
                        pos_x = spikeTimestamp_3Brain(idx_3Brain_neighbors(j)).positionXY(1);
                        pos_y = spikeTimestamp_3Brain(idx_3Brain_neighbors(j)).positionXY(2);
                    else
                        y = x*0;
                        pos_x = nan;
                        pos_y = nan;
                    end
                    line(x,y);      legend(sprintf('x=%g, y=%g', pos_x, pos_y))
                end
                
                display('Press a button to continue...')
                pause()
                close(fig1, fig2, fig3, fig4)
            end
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 4              % compare result from Kilosort vs 3Brain (using 3Brain as reference)    

        file_excel = sprintf('%s%smetadata.csv', folderPath, delimiter);
        T = readtable(file_excel);         
        idx_set = activeClusters;
        unitLabel = nan(numel(idx_set),2);
        Kilosort_matrix = struct('ch',cell(N_NEIGHBOR*2+1,N_NEIGHBOR*2+1),'clusters',cell(N_NEIGHBOR*2+1,N_NEIGHBOR*2+1),'scores',cell(N_NEIGHBOR*2+1,N_NEIGHBOR*2+1),'waveIDX',cell(N_NEIGHBOR*2+1,N_NEIGHBOR*2+1),'position',cell(N_NEIGHBOR*2+1,N_NEIGHBOR*2+1));

        % Load the information for the neighboring electrodes from the Kilosort data
        for i=1:numel(idx_set)                                              % across all the interesting units identified by 3Brain
            idx_target = spikeTimestamp_3Brain(idx_set(i)).cluster;         % label of the target cluster i (0-4095)         
            for j=-1*N_NEIGHBOR:1:N_NEIGHBOR                                % look at neighboring electrodes
                for l=-1*N_NEIGHBOR:1:N_NEIGHBOR                            % look at neighboring electrodes
                    idx = idx_target + MATRIX_SIZE*j + l;                   % index of the neighboring electrode             
                    jj = N_NEIGHBOR+1+j;    ll = N_NEIGHBOR+1+l;            % (revised October 2026) position in "Kilosort_matrix" (before, "2+j" and "2+l": only valid for N_NEIGHBOR = 1)
                    Kilosort_matrix(jj,ll).ch = idx;                        % square around idx_target
                    cluster_pos = [];  position = [];  scores = []; wave_idx = [];
                    m = 1;                                                  % to count the number of different units identified by Kilosort in the electrode "idx"
                    for k=1:size(spikeTimestamp_Kilosort,2)                 % across all units identified by Kilosort
                        if idx == spikeTimestamp_Kilosort(k).position       % kilosort has identified a unit in the neighboring electrode "idx"
                            cluster_pos = vertcat(cluster_pos, k);          % position in spikeTimestamp of the clusters from Kilosort in the channel "idx"  
                            position = vertcat(position, spikeTimestamp_Kilosort(k).positionXY);
                            idx_waveform = 0;            count=1;
                            while idx_waveform<1
                                if T{count,NEW_CLUSTER_ID}==spikeTimestamp_Kilosort(k).cluster
                                    idx_waveform = count;
                                else
                                    count=count+1;
                                end
                            end
                            wave_idx = vertcat(wave_idx,idx_waveform);
                            if spikeTimestamp_Kilosort(k).position ~= T{idx_waveform,PEAK_CHANNEL}
                                error('Error finding the cluster - the electrode position does not match')
                            end
                            
                            score_aux = T{idx_waveform,UNIT_TYPE};
                            switch score_aux
                                case 0
                                    score_unit = 'noise';
                                case 1
                                    score_unit = 'good';
                                case 2
                                    score_unit = 'MUA';
                                case 3
                                    score_unit = 'non-somatic';
                            end
                            scores{m} = score_unit;
                            m = m +1;
                        end
                    end
                    Kilosort_matrix(jj,ll).clusters = cluster_pos;
                    Kilosort_matrix(jj,ll).position = position;
                    Kilosort_matrix(jj,ll).scores = scores;
                    Kilosort_matrix(jj,ll).waveIDX = wave_idx;
                    clear cluster_pos position scores wave_idx
                 end
            end

            % Plot the electrode map from 3Brain data
            fig1 = figure(1);
            pos_x = spikeTimestamp_3Brain(idx_set(i)).positionXY(1)+1;     % (revised October 2026: +1, positionXY is 0-based -as in options 1 and 2-; before, the 
            pos_y = spikeTimestamp_3Brain(idx_set(i)).positionXY(2)+1;     %  electrodes of the first row/column crashed, and the others were shifted by one)
            map(pos_y,pos_x)=1;
            heatmap(map,'XDisplayLabels',NaN*ones(size(map,1),1),'YDisplayLabels',NaN*ones(size(map,1),1),'Colormap',cool);     title('3Brain')

            % plot raster plot (3Brain vs Kilosort)
            fprintf(sprintf('\n Unit number(cluster id)   %g', idx_target));
            fprintf(sprintf('\n It has %g spikes \n', spikeTimestamp_3Brain(idx_set(i)).nUnits));

            fig2 = figure(2);
            subplot(4,1,4)
            bar(spikes_bin(:,2)*K_MSEC_SEC,spikes_bin(:,1)/(BIN*K_MSEC_SEC));      ylabel('Pop Firing - 3Brain');      xlabel('Time (sec)', 'FontSize', 18);       xlim([0 spikeTimestamp_3Brain(1).window(2)]);
            
            subplot(4,1,[1,2,3])
            data = spikeTimestamp_3Brain(idx_set(i));
            plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1),'|',"MarkerEdgeColor",'r');            hold on
            plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+1,'|',"MarkerEdgeColor",'r');            hold on

            shift_x = N_NEIGHBOR*2+1;
            for j=1:N_NEIGHBOR*2+1
                for l=1:N_NEIGHBOR*2+1
                    if numel(Kilosort_matrix(j,l).clusters) > 0 && j==N_NEIGHBOR+1 && l==N_NEIGHBOR+1               % electrode of interest                     
                        for k=1:numel(Kilosort_matrix(j,l).clusters)
                            data = spikeTimestamp_Kilosort(Kilosort_matrix(j,l).clusters(k));
                            plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+shift_x,'|',"MarkerEdgeColor",'b');            hold on
                            plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+1,'|',"MarkerEdgeColor",'b');            hold on
                            shift_x = shift_x + 1;
                        end
                    else                                                                                            % neighboring electrodes
                        for k=1:numel(Kilosort_matrix(j,l).clusters)
                            data = spikeTimestamp_Kilosort(Kilosort_matrix(j,l).clusters(k));
                            plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+shift_x,'|',"MarkerEdgeColor",'k');            hold on
                            shift_x = shift_x + 1;
                        end
                    end
                end
            end
            ylim([0 shift_x+1]);        xlim([0 spikeTimestamp_3Brain(1).window(2)]);             xlabel('Time (sec)');            ylabel('Cluster');
            str = {'\color{red}3Brain - \color{blue}Kilosort (target) - \color{black}Kilosort (neighbors)'};      title(str);

            % plot waveforms (3Brain vs Kilosort)
            fig3 = figure(3);
            cond = 1;
            waveformData = Indiv_Cond_3Brain{cond,2};
            dacFactor = dacFactors(cond);
            offsetValue = offsetValues(cond);

            traces = waveformData{spikeTimestamp_3Brain(idx_set(i)).cluster+1};
            x = double((1:waveformLength)) / samplingRate;
            y = offsetValue + dacFactor * mean(traces,1);
            dy = std(offsetValue + dacFactor*traces,0,1);
            
            subplot(1,2,1)
            fill([x';flipud(x')],[(y-dy)';flipud((y+dy)')],[0 1 1],'linestyle','none');     line(x,y);      title('3Brain')

            subplot(1,2,2)
            center = N_NEIGHBOR+1;
            if numel(Kilosort_matrix(center,center).clusters)>0
                for k=1:numel(Kilosort_matrix(center,center).clusters)
                    idx=Kilosort_matrix(center,center).waveIDX(k); 
                    xx = double((1:numel(waveforms_Kilosort{idx}))') / samplingRate; 
                    plot(xx,waveforms_Kilosort{idx},'DisplayName',Kilosort_matrix(center,center).scores{k});       hold on;
                    legend('-DynamicLegend');       legend('show');
                end
            end
            title(sprintf('Kilosort'));

            % plot waveforms (Kilosort)
            fig4 = figure(4);
            neighbor_n =  N_NEIGHBOR*2+1; 
            for j=1:neighbor_n
                for l=1:neighbor_n
                    subplot(neighbor_n,neighbor_n,(j-1)*neighbor_n+l)
                    if numel(Kilosort_matrix(l,j).clusters)>0
                        for k=1:numel(Kilosort_matrix(l,j).clusters)
                            idx=Kilosort_matrix(l,j).waveIDX(k); 
                            xx = double((1:numel(waveforms_Kilosort{idx}))') / samplingRate; 
                            plot(xx,waveforms_Kilosort{idx});       hold on;
                        end
                        legend(sprintf('x=%g, y=%g', Kilosort_matrix(l,j).position(1,1),Kilosort_matrix(l,j).position(1,2)))
                    end
                end
            end            
            map(pos_y,pos_x)=0;

            % Manual classification of 3Brain units
            prompt = "Classification? (G=good / W=weak / N=noise):  ";
            label = input(prompt,"s"); 
            unitLabel(i,1) = idx_target;
            switch label 
                case 'g'
                    unitLabel(i,2) = 1;
                case 'w'
                    unitLabel(i,2) = 2;
                case 'n'
                    unitLabel(i,2) = -1;
                otherwise
                    unitLabel(i,2) = 0;
            end

            file_save = sprintf('%s%s%s.mat', folderPath, delimiter, file_output);   % (revised October 2026: separate output file, see INPUTS)
            if isfile(file_save)
                save(file_save,'unitLabel','-append');
            else
                save(file_save,'unitLabel');                                % (the first time, the file is created)
            end
            pause();        close(fig1, fig2, fig3, fig4)
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 5              % compare result from Kilosort vs 3Brain (using 3Brain as reference)    

        file_save = sprintf('%s%s%s.mat', folderPath, delimiter, file_output);   % (revised October 2026: the file saved by option 4; before, "xxx.mat", which no code creates)
        load(file_save,'unitLabel');
        file_excel = sprintf('%s%smetadata.csv', folderPath,delimiter);
        T = readtable(file_excel);
        
        aux = find(unitLabel(:,2)==1);
        units_3Brain.good = numel(aux);
        aux = find(unitLabel(:,2)==2);
        units_3Brain.weak = numel(aux);
        aux = find(unitLabel(:,2)==-1);
        units_3Brain.noise = numel(aux);
        aux = find(unitLabel(:,2)==0);
        units_3Brain.undetermined = numel(aux);
        units_3Brain.total = size(unitLabel,1) - sum(isnan(unitLabel(:,2)));


        % aux_kilosort = nan(size(spikeTimestamp_Kilosort,2),1);
        % for k=1:size(spikeTimestamp_Kilosort,2)
        %     idx_waveform = 0;            count=1;
        %     while idx_waveform<1
        %           if T{count,NEW_CLUSTER_ID}==spikeTimestamp_Kilosort(k).cluster
        %              idx_waveform = count;
        %           else
        %              count=count+1;
        %           end
        %     end
        %     score_aux = T{idx_waveform,UNIT_TYPE};
        %     aux_kilosort(k) = score_aux;
        % end
        aux_kilosort = nan(size(T,1),1);
        for k=1:size(T,1)
            aux_kilosort(k) = T{k,UNIT_TYPE};
        end

        aux = find(aux_kilosort==1);
        units_Kilosort.good = numel(aux);
        aux = find(aux_kilosort==2);
        units_Kilosort.mua = numel(aux);
        aux = find(aux_kilosort==3);
        units_Kilosort.nonsomatic = numel(aux);
        aux = find(aux_kilosort==0);
        units_Kilosort.noise = numel(aux);
        units_Kilosort.total = k;

        a=1;
        save(file_save,'units_3Brain','units_Kilosort','-append')

    end
end

