function Brain3_processedData()

% function to handle processed data from 3Brain system
% (October 2026: moved to the folder "legacy"; superseded by "ShapeSpikeTrainData_3Brain" + "Compare_Kilosort_vs_3Brain". It only reads BrainWave v5 files,
%  option 4 starts at the unit 3501 and option 3 has a plot call that fails; it was not revised)

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

folder_3Brain = '/path/to/recording_folder';
file_data_3Brain = 'recording_name_3Brain';
folder_Kilosort = '/path/to/recording_folder';
file_data_Kilosort = 'recording_name';
file_output = 'recording_name kilosort vs 3Brain';

dataDurationSec = 900;                      % duration in seconds of the portion to extract and plot
dataStartPositionSec = 0;                   % start position in seconds of the portion to extract and plot

option_plot = 2;                            % =0, it's the first time we run it; it will processed the data (loading the waveform from 3Brain is an heavy process) and save the results in a file
                                            % =1, raster plot of all clusters   
                                            % =2, waveform and raster plot of the most active electrodes, one-by-one
                                            % =3, compare result from Kilosort vs 3Brain (using kilosort as reference)    
                                            % =4, compare result from Kilosort vs 3Brain (using 3Brain as reference)
                                            % =5, summarize the comparison Kilosort vs 3Brain        

K_MSEC_SEC = 1/1000;                        % factor to rescale from msec to sec
MATRIX_SIZE = 64;
NEW_CLUSTER_ID = 4;
UNIT_TYPE = 2;
PEAK_CHANNEL = 6;
WAVEFORM_TEMPLATE = 5;
N_CLUSTER_PLOT = 10;
BIN = 0.25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data 

file_input = sprintf('%s%s%s.bxr', folder_3Brain, delimiter, file_data_3Brain);
samplingRate = h5readatt(file_input, '/','SamplingRate');
chIdxs = h5read(file_input, '/Well_A1/StoredChIdxs');
spikeTOC = h5read(file_input, '/Well_A1/SpikeTOC');
toc = h5read(file_input, '/TOC');

info = h5info(file_input);
spikeSortingPerformed = 0;
for i=1:size(info.Groups.Datasets,1)
    if strcmp(info.Groups.Datasets(i).Name,'SpikeUnits')
        spikeSortingPerformed = 1;
    end
end
 
startFrame = round(dataStartPositionSec * samplingRate);
nFrames = round(dataDurationSec * samplingRate);

minDigitalValue = h5readatt(file_input, '/','MinDigitalValue');
maxDigitalValue = h5readatt(file_input, '/','MaxDigitalValue');
minAnalogValue = h5readatt(file_input, '/','MinAnalogValue');
maxAnalogValue = h5readatt(file_input, '/','MaxAnalogValue');
dacFactor = (maxAnalogValue - minAnalogValue) / (maxDigitalValue - minDigitalValue);
offsetValue = minAnalogValue - dacFactor * minDigitalValue;

file_name = sprintf('%s%swaveforms.mat',folder_Kilosort,delimiter);
load(file_name, 'waveforms');
waveforms_Kilosort = waveforms;

if option_plot==0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% From the given start position and duration (in frames), find the corresponding range of spike positions using the TOC

aux1 = find(toc(1,:)<startFrame);
if numel(aux1)>0
    tocStartIdx = aux1(end);
else
    tocStartIdx = 1;
end

aux2 = find(toc(1,:)>startFrame+nFrames);
if numel(aux2)>0
    tocEndIdx = aux2(1);
else
    tocEndIdx = size(toc,2);
end

spikeStartPosition = spikeTOC(tocStartIdx);
spikeStartPosition = max(spikeStartPosition,1);
spikeEndPosition = spikeTOC(tocEndIdx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Collect the required spike data and waveform data

spikeData = h5read(file_input, '/Well_A1/SpikeTimes');
spikeDataTimestamps = spikeData(spikeStartPosition:spikeEndPosition);
spikeDataCh = h5read(file_input, '/Well_A1/SpikeChIdxs');
spikeDataChIdxs = spikeDataCh(spikeStartPosition:spikeEndPosition);
dataLength = spikeEndPosition - spikeStartPosition;

waveformLength = int64(h5readatt(file_input,'/Well_A1/SpikeForms','Wavelength'));
wave = h5read(file_input, '/Well_A1/SpikeForms');
spikeDataWaveforms = wave(((spikeStartPosition-1)*waveformLength+1):spikeEndPosition*waveformLength);

if spikeSortingPerformed
    spikeChUnit = h5read(file_input, '/Well_A1/SpikeUnits');
    spikeDataChUnits = spikeChUnit(spikeStartPosition:spikeEndPosition);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the raster dataset with spike times and convert spike times from frame position to seconds for each channel
% Create the structure "spikeTimestamp"

if spikeSortingPerformed
    n_cluster = double(max(spikeDataChUnits));
else
    n_cluster = double(numel(chIdxs));
end

if option_plot<4 

    spikeTimestamp = struct('idx_toc',cell(1,n_cluster),'timeStamp',cell(1,n_cluster),'nUnits',cell(1,n_cluster),'recLength',cell(1,n_cluster),'window',cell(1,n_cluster),'cluster',cell(1,n_cluster),'position',cell(1,n_cluster),'t_step',cell(1,n_cluster),'positionXY',cell(1,n_cluster));

for i=1:n_cluster                                                           % across all different clusters
    if spikeSortingPerformed
        aux_idx = find(spikeDataChUnits==i);                                % index of all units in cluster i
        aux = spikeDataChIdxs(aux_idx);
        if numel(aux)>0
            cluster_id = double(mode(aux));
            pos_x = fix(cluster_id/MATRIX_SIZE)+1;
            pos_y = rem(cluster_id, MATRIX_SIZE)+1;
         else
            cluster_id = nan;
            pos_x = nan;
            pos_y = nan;
        end
    else
        aux_idx = find(spikeDataChIdxs==chIdxs(i));                         % index of all units in channel i
        cluster_id = chIdxs(i);
        pos_x = fix((i-1)/MATRIX_SIZE)+1;
        pos_y = rem(i-1, MATRIX_SIZE)+1;
    end
    n_units = numel(aux_idx);                                               % number of total units in cluster i
 
    if n_units>0                                                            % there is at least one unit from cluser i 
       spikes = spikeDataTimestamps(aux_idx); 
       spikes = double(spikes)/samplingRate;         
       spikeTimestamp(i).timeStamp = spikes;                                % spike time stamp for cluster i, in sec
       spikeTimestamp(i).nUnits = n_units;                                  % number of units in the window of interest
       spikeTimestamp(i).idx_toc = aux_idx;                                 % index of all units in channel/cluster i (to be used as index for "spikeDataChIdxs")
    else                    
        spikeTimestamp(i).timeStamp = [];                                   % no units
        spikeTimestamp(i).nUnits = 0;
        spikeTimestamp(i).idx_toc = [];
    end
   
    spikeTimestamp(i).recLength = (spikeDataTimestamps(end)-spikeDataTimestamps(1))/samplingRate;       % length of the window, in sec 
    spikeTimestamp(i).window = [dataStartPositionSec dataStartPositionSec+dataDurationSec];             % initial and last point of the window, in sec
    spikeTimestamp(i).cluster = cluster_id;                                 % label of the cluster i (0-4095)
    spikeTimestamp(i).position = nan;                                       % position of the peak amplitude signal for the cluster i
    spikeTimestamp(i).t_step = 1/samplingRate/K_MSEC_SEC;                   % temporal resolution, in msec
    spikeTimestamp(i).positionXY = [pos_x; pos_y];                          % coordinates (X,Y) of the peak amplitude signal for the cluster i

    clear spikes n_units aux_idx aux cluster_id
end
spikeTimestamp_3Brain = spikeTimestamp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the waveform dataset with individual traces for each electrode/unit

for idx=1:n_cluster
    aux_idx = spikeTimestamp(idx).idx_toc;                                  % index of all units in channel/cluster i (to be used as index for "spikeDataChIdxs")
    k=1;
    waveform1 = [];
    waveform2 = [];
    for i=1:dataLength
        % if spikeDataChIdxs(i) == chIdx && startFrame <= spikeDataTimestamps(i) && spikeDataTimestamps(i) < startFrame + nFrames
        if ismember(i,aux_idx) && startFrame <= spikeDataTimestamps(i) && spikeDataTimestamps(i) < startFrame + nFrames
           % waveformData = horzcat(waveformData, double(spikeDataWaveforms(i*waveformLength:(i+1)*waveformLength))); 
           waveform1(k,:) = double(spikeDataWaveforms(i*waveformLength:(i+1)*waveformLength-1));
           waveform2(k,:) = double(spikeDataWaveforms(((i-1)*waveformLength+1):i*waveformLength));
           k=k+1;
        end
    end
    waveformData1{idx} = waveform1;
    waveformData2{idx} = waveform2;
    clear k waveform1 waveform2 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the raster dataset with spike times and convert spike times from frame position to seconds for each channel

cont_spikes = zeros(n_cluster,1);
for i=1:n_cluster
    cont_spikes(i) = spikeTimestamp(i).nUnits;
end

threshold_activity = quantile(cont_spikes,0.975);
activeClusters = find(cont_spikes>threshold_activity);

map = nan(MATRIX_SIZE,MATRIX_SIZE); 
for idx=1:numel(activeClusters)
    
    chIdx = activeClusters(idx)-1;                            % identifier of channel to be analysed
    pos_x = spikeTimestamp(chIdx+1).positionXY(1);
    pos_y = spikeTimestamp(chIdx+1).positionXY(2);
    map(pos_y,pos_x)=-1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the population firing rate data

length = dataDurationSec;                                                   % length of the recording
n_bins = ceil(length/(BIN));                                                % total number of bins to process the data
spikes_bin = nan(n_bins,4);

for i=1:n_bins                                                              % across all the bins
    aux_spikes = 0;
    spikes_bin(i,2) = dataStartPositionSec+(i-1)*BIN;                       % left side of the bin i, in sec
    spikes_bin(i,3) = spikes_bin(i,2) + BIN;                                % right side of the bin i, in sec
    spikes_bin(i,4) = (spikes_bin(i,2)+spikes_bin(i,3))/2;                  % midpoint of the bin i, in sec
    for ch=1:n_cluster                                                      % across all clusters
        aux1 = find(spikeTimestamp_3Brain(ch).timeStamp >= spikes_bin(i,2));
        aux2 = find(spikeTimestamp_3Brain(ch).timeStamp < spikes_bin(i,3));
        aux = intersect(aux1, aux2);                                        % index of all spikes of cluster ch within the bin i
        if numel(aux) > 0
           aux_spikes = aux_spikes + numel(aux);                            % total number of spikes in the bin i (cumulative across all clusters)
        end
        clear aux aux1 aux2
    end
    spikes_bin(i,1) = aux_spikes/n_cluster;                                 % total number of spikes in bin i, per channel
    clear aux_spikes
end

end
else
        file_save = sprintf('%s%s%s.mat', folder_Kilosort, delimiter, file_output);
        load(file_save,'spikeTimestamp_3Brain','spikeTimestamp_Kilosort','waveformData2','spikes_bin','cont_spikes','activeClusters','waveformLength','samplingRate','n_cluster','map');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting option

switch option_plot                      
    
    case 0      % Process the data and save the result (do not plot anything)
        
        clear spikeTimestamp
        file_name = sprintf('%s%s%s_cleaned.mat', folder_Kilosort, delimiter,file_data_Kilosort);
%        file_name = sprintf('%s%s%s_QCmanual.mat', folder_Kilosort, delimiter,file_data_Kilosort);
        load(file_name, 'spikeTimestamp');
        spikeTimestamp_Kilosort = spikeTimestamp; clear spikeTimestamp;
        file_save = sprintf('%s%s%s.mat', folder_Kilosort, delimiter, file_output);
        save(file_save,'spikeTimestamp_3Brain','spikeTimestamp_Kilosort','waveformData2','spikes_bin','cont_spikes','activeClusters','waveformLength','samplingRate','n_cluster','map');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1        % raster plot of all clusters 
        
        figure(1)
        for k=1:n_cluster
            data = spikeTimestamp_3Brain(k);
            plot(data.timeStamp, ones(data.nUnits,1)+k,'|',"MarkerEdgeColor",[0 0 0]);
            hold on;            clear data
        end 
        ylim([0 1+k]);        xlabel('Time (sec)');        ylabel('Neuron');    title('ALL electrodes')

        figure(2)
        histogram(cont_spikes)

        figure(3)
        for k=1:numel(activeClusters)
            data = spikeTimestamp_3Brain(activeClusters(k));
            plot(data.timeStamp, ones(data.nUnits,1)+k,'|',"MarkerEdgeColor",[0 0 0]);
            hold on;            clear data
        end 
        ylim([0 1+k]);        xlabel('Time (sec)');        ylabel('Neuron');    title('Most active electrodes')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2       % waveform and raster plot of the most active electrodes, one-by-one

        x = double((1:waveformLength)) / samplingRate;
        spikeTimestamp=spikeTimestamp_3Brain;
        for idx=1:numel(activeClusters)

            chIdx = activeClusters(idx)-1;                                  % identifier of channel to be analysed
            fprintf(sprintf('\n Electrode number   %g', spikeTimestamp(chIdx+1).cluster));
            fprintf(sprintf('\n It has %g units \n', spikeTimestamp(chIdx+1).nUnits));

            traces = waveformData2{activeClusters(idx)};
            y = offsetValue + dacFactor * mean(traces,1);
            dy = std(offsetValue + dacFactor*traces,0,1);

            fig1=figure();
            subplot(2,2,[1,2])
            for k=1:numel(activeClusters)
                data = spikeTimestamp(activeClusters(k));
                plot(data.timeStamp, ones(data.nUnits,1)+k,'|',"MarkerEdgeColor",[0 0 0]);
                hold on;                clear data
            end 
            data = spikeTimestamp(activeClusters(idx));
            plot(data.timeStamp, ones(data.nUnits,1)+idx,'|',"MarkerEdgeColor","r");
            clear data;            ylim([0 1+k]);            xlabel('Time (sec)');            ylabel('Neuron')

            subplot(2,2,3)
            for i=1:size(traces ,1)
                yy = offsetValue + dacFactor * traces(i,:);
                plot(x,yy);                hold on
            end

            subplot(2,2,4)
            fill([x';flipud(x')],[(y-dy)';flipud((y+dy)')],[0 1 1],'linestyle','none')
            line(x,y);

            fig2 = figure(2);
            pos_x = spikeTimestamp(chIdx+1).positionXY(1);
            pos_y = spikeTimestamp(chIdx+1).positionXY(2);
            map(pos_y,pos_x)=1;
            heatmap(map,'XDisplayLabels',NaN*ones(size(map,1),1),'YDisplayLabels',NaN*ones(size(map,1),1),'Colormap',cool);
            
            pause()
            map(pos_y,pos_x)=0;
            close(fig1);            close(fig2);
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3              % compare result from Kilosort vs 3Brain (using kilosort as reference)     

        % file_name = sprintf('%s\\%s_QCmanual.mat', folder_Kilosort, file_data_Kilosort);
        % file_name = sprintf('%s\\%s_cleaned.mat', folder_Kilosort, file_data_Kilosort);
        % load(file_name, 'spikeTimestamp');
        % spikeTimestamp_Kilosort = spikeTimestamp;
        file_excel = sprintf('%s%smetadata.csv', folder_Kilosort, delimiter);
        T = readtable(file_excel);
         
        for i=1:size(spikeTimestamp_Kilosort,2)

            cluster_idx = spikeTimestamp_Kilosort(i).position;
            extra_cluster = [];
            for k=1:size(spikeTimestamp_Kilosort,2)
                if spikeTimestamp_Kilosort(k).position == cluster_idx
                    extra_cluster = vertcat(extra_cluster,k);
                end
            end

            idx_3Brain = nan;
            idx_3Brain_neighbors = nan(9,1);
            for k=1:size(spikeTimestamp_3Brain,2)
                if spikeTimestamp_3Brain(k).cluster == cluster_idx
                    idx_3Brain = k;
                    idx_3Brain_neighbors(1) = k-MATRIX_SIZE-1;
                    idx_3Brain_neighbors(2) = k-MATRIX_SIZE;
                    idx_3Brain_neighbors(3) = k-MATRIX_SIZE+1;
                    idx_3Brain_neighbors(4) = k-1;
                    idx_3Brain_neighbors(5) = k;                % selected
                    idx_3Brain_neighbors(6) = k+1;
                    idx_3Brain_neighbors(7) = k+MATRIX_SIZE-1;
                    idx_3Brain_neighbors(8) = k+MATRIX_SIZE;
                    idx_3Brain_neighbors(9) = k+MATRIX_SIZE+1;

                    aux = find(idx_3Brain_neighbors<0);
                    idx_3Brain_neighbors(aux) = nan;
                    aux = find(idx_3Brain_neighbors>4095);
                    idx_3Brain_neighbors(aux) = nan;
                end
            end
            if isnan(idx_3Brain)
                error('\n The cluster could NOT be found in the 3Brain dataset \n')
            end

            % plot raster plot (3Brain vs Kilosort)
            fig1 = figure(1);
            n_clusters_plot = min(size(spikeTimestamp_Kilosort,2),N_CLUSTER_PLOT);
            for k=1:n_clusters_plot
                data = spikeTimestamp_Kilosort(k);
                plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+k,'|',"MarkerEdgeColor",[0 0 0]);
                hold on;                clear data
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
                        plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+shift_x, 'Kilosort (extra)','|');                        hold on
                        plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+shift_x0,'|',"MarkerEdgeColor",'m');                        hold on
                        shift_x = shift_x + 1;
                    end
                end
            end

            data_3Brain = spikeTimestamp_3Brain(idx_3Brain);
            plot(data_3Brain.timeStamp, ones(data_3Brain.nUnits,1)+shift_x,'|',"MarkerEdgeColor",'b');            hold on
            plot(data_3Brain.timeStamp, ones(data_3Brain.nUnits,1)+shift_x0,'|',"MarkerEdgeColor",'b');            ylim([0 shift_x+4]);            xlabel('Time (sec)');            ylabel('Cluster');
            str = {'\color{black}Kilosort - \color{red}Kilosort - \color{magenta}Kilosort (extra) - \color{blue}3Brain'};
            title(str)

            % plot raster plot (3Brain; selected vs neighboring electrodes)
            fig2 = figure(2);
            data = spikeTimestamp_Kilosort(i);
            plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1),'|',"MarkerEdgeColor",'r');            hold on

            shift_x = 2;
            data_3Brain = spikeTimestamp_3Brain(idx_3Brain);
            plot(data_3Brain.timeStamp, ones(data_3Brain.nUnits,1)+shift_x,'|',"MarkerEdgeColor",'b');            hold on
            shift_x = shift_x + 1;
            for j=1:9
                if idx_3Brain_neighbors(j)<4097 && idx_3Brain_neighbors(j)>0
                    data_3Brain = spikeTimestamp_3Brain(idx_3Brain_neighbors(j));
                else
                    data_3Brain.timeStamp = 0;
                    data_3Brain.nUnits = 1;
                end
                if j==5
                    plot(data_3Brain.timeStamp, ones(data_3Brain.nUnits,1)+shift_x,'|',"MarkerEdgeColor",'b');
                else
                    plot(data_3Brain.timeStamp, ones(data_3Brain.nUnits,1)+shift_x,'|',"MarkerEdgeColor",[0 0 0]);
                    plot(data_3Brain.timeStamp, ones(data_3Brain.nUnits,1)+2,'|',"MarkerEdgeColor",'m');
                end
                shift_x = shift_x + 1;                hold on
            end
            ylim([0 shift_x+1]);            xlabel('Time (sec)');            ylabel('Cluster');
            str = {'\color{black}3Brain - \color{magenta}3Brain (neighbors) - \color{blue}3Brain (target) - \color{red}Kilosort'};
            title(str)
        
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
            plot(xx,waveforms_Kilosort{idx_waveform}, 'r');
            title(sprintf('Kilosort (%s)', score_unit));
        
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
            x = double((1:waveformLength)) / samplingRate;
            traces = waveformData2{idx_3Brain};
            y = offsetValue + dacFactor * mean(traces,1);
            dy = std(offsetValue + dacFactor*traces,0,1);
            fill([x';flipud(x')],[(y-dy)';flipud((y+dy)')],[0 1 1],'linestyle','none')
            line(x,y);
            title('3Brain')

            % plot waveforms (3Brain; selected vs neighboring electrodes)
            fig4 = figure(4); 

            for j=1:9
                subplot(3,3,j)
                idx_ch = idx_3Brain_neighbors(j);
                if ~isnan(idx_ch)
                    traces = waveformData2{idx_ch};
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
                line(x,y);
                legend(sprintf('x=%g, y=%g', pos_x, pos_y))
            end
            display('Press a button to continue...')
            pause()
            close(fig1, fig2, fig3, fig4)
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 4              % compare result from Kilosort vs 3Brain (using 3Brain as reference)    

        % file_name = sprintf('%s\\%s_QCmanual.mat', folder_Kilosort, file_data_Kilosort);
        % file_name = sprintf('%s\\%s_cleaned.mat', folder_Kilosort, file_data_Kilosort);
        % load(file_name, 'spikeTimestamp');
        % spikeTimestamp_Kilosort = spikeTimestamp;
        file_excel = sprintf('%s%smetadata.csv', folder_Kilosort, delimiter);
        file_save = sprintf('%s%s%s.mat', folder_Kilosort, delimiter, file_output);
        T = readtable(file_excel);
        % load(file_save,'spikeTimestamp_3Brain','spikeTimestamp_Kilosort','waveformData2','spikes_bin','unitLabel')
        % save(file_save,'spikeTimestamp_3Brain','spikeTimestamp_Kilosort','waveformData2','spikes_bin')

        option_set = 1;
        min_freq = 0.025;
        map = nan(MATRIX_SIZE,MATRIX_SIZE); 
        switch option_set
            case 1      % all clusters
                n_set = n_cluster;
                idx_set = 1:n_cluster;

                for k=1:n_cluster
                    if spikeTimestamp_3Brain(k).nUnits>round(min_freq*dataDurationSec)
                        map(spikeTimestamp_3Brain(k).positionXY(2),spikeTimestamp_3Brain(k).positionXY(1))=-1;
                    end
                end
    
            case 2      % only most active electrodes
                n_set = numel(activeClusters);
                idx_set = activeClusters;
        end
      
        unitLabel = nan(n_set,2);
        for i=3501:n_set
        % for i=1:n_set
            if spikeTimestamp_3Brain(idx_set(i)).nUnits>round(min_freq*dataDurationSec)

            idx_target = spikeTimestamp_3Brain(idx_set(i)).cluster;        % label of the target cluster i (0-4095)
            
            for j=-1:1:1
                for l=-1:1:1
                    idx = idx_target + MATRIX_SIZE*j + l;
                    Kilosort_matrix(2+j,2+l).ch = idx;                     % square around idx_target
                    cluster_pos = [];  position = [];  scores = []; wave_idx = [];
                    m = 1;
                    for k=1:size(spikeTimestamp_Kilosort,2)
                        if idx == spikeTimestamp_Kilosort(k).position
                            cluster_pos = vertcat(cluster_pos, k);         % position in spikeTimestamp of the clusters from Kilosort in the channel "idx"  
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
                    Kilosort_matrix(2+j,2+l).clusters = cluster_pos;
                    Kilosort_matrix(2+j,2+l).position = position;
                    Kilosort_matrix(2+j,2+l).scores = scores;
                    Kilosort_matrix(2+j,2+l).waveIDX = wave_idx;
                    clear cluster_pos position scores
                 end
            end

            % Plot the electrode map
            fig1 = figure(1);
            pos_x = spikeTimestamp_3Brain(idx_set(i)).positionXY(1);
            pos_y = spikeTimestamp_3Brain(idx_set(i)).positionXY(2);
            map(pos_y,pos_x)=1;
            heatmap(map,'XDisplayLabels',NaN*ones(size(map,1),1),'YDisplayLabels',NaN*ones(size(map,1),1),'Colormap',cool);
            title('3Brain')

            % plot raster plot 
            fprintf(sprintf('\n Unit number   %g', idx_target));
            fprintf(sprintf('\n It has %g units \n', spikeTimestamp_3Brain(idx_set(i)).nUnits));

            fig2 = figure(2);
            subplot(4,1,4)
            bar(spikes_bin(:,2),spikes_bin(:,1)/BIN);
            xlabel('Time (sec)', 'FontSize', 18)
            xlim([0 dataDurationSec]);
            
            subplot(4,1,[1,2,3])
            data = spikeTimestamp_3Brain(idx_set(i));
            plot(data.timeStamp, ones(data.nUnits,1),'|',"MarkerEdgeColor",'r');            hold on
            plot(data.timeStamp, ones(data.nUnits,1)+1,'|',"MarkerEdgeColor",'r');            hold on

            shift_x = 3;
            for j=1:3
                for l=1:3
                    if numel(Kilosort_matrix(j,l).clusters) && j==2 && l==2 >0
                        for k=1:numel(Kilosort_matrix(j,l).clusters)
                            data = spikeTimestamp_Kilosort(Kilosort_matrix(j,l).clusters(k));
                            plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+shift_x,'|',"MarkerEdgeColor",'b');            hold on
                            plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+1,'|',"MarkerEdgeColor",'b');            hold on
                            shift_x = shift_x + 1;
                        end
                    else
                        for k=1:numel(Kilosort_matrix(j,l).clusters)
                            data = spikeTimestamp_Kilosort(Kilosort_matrix(j,l).clusters(k));
                            plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+shift_x,'|',"MarkerEdgeColor",'k');            hold on
                            shift_x = shift_x + 1;
                        end
                    end
                end
            end
            ylim([0 shift_x+1]);        xlim([0 dataDurationSec]);             xlabel('Time (sec)');            ylabel('Cluster');
            str = {'\color{red}3Brain - \color{blue}Kilosort (target) - \color{black}3Brain (neighbors)'};
            title(str)

            % plot waveforms (3Brain vs Kilosort)
            fig3 = figure(3);
            traces = waveformData2{idx_set(i)};
            x = double((1:waveformLength)) / samplingRate;
            y = offsetValue + dacFactor * mean(traces,1);
            dy = std(offsetValue + dacFactor*traces,0,1);

            subplot(1,2,1)
            fill([x';flipud(x')],[(y-dy)';flipud((y+dy)')],[0 1 1],'linestyle','none')
            line(x,y);
            title('3Brain')

            subplot(1,2,2)
            if numel(Kilosort_matrix(2,2).clusters)>0
                for k=1:numel(Kilosort_matrix(2,2).clusters)
                    idx=Kilosort_matrix(2,2).waveIDX(k); 
                    xx = double((1:numel(waveforms_Kilosort{idx}))') / samplingRate; 
                    plot(xx,waveforms_Kilosort{idx},'DisplayName',Kilosort_matrix(2,2).scores{k});       hold on;
                    legend('-DynamicLegend');
                    legend('show');
                end
            end
            title(sprintf('Kilosort'));

            % plot waveforms (Kilosort)
            fig4 = figure(4);
            for j=1:3
                for l=1:3
                    subplot(3,3,(j-1)*3+l)
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
            fig1.Position(1:2) = [10 10];
            fig2.Position(1:2) = [10 350];
            fig3.Position(1:2) = [500 350];
            fig4.Position(1:2) = [1000 10];            
            map(pos_y,pos_x)=0;

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
            save(file_save,'unitLabel','-append');
            % pause()
            close(fig1, fig2, fig3, fig4)
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 5              % compare result from Kilosort vs 3Brain (using 3Brain as reference)    

        file_save = sprintf('%s%s%s.mat', folder_Kilosort, delimiter,file_output);
        load(file_save,'unitLabel');
        file_excel = sprintf('%s%smetadata.csv', folder_Kilosort,delimiter);
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



