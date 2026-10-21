function plotMEArecording()

% Function to plot data from MEA recording

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

folder = '/path/to/recording_folder';
file_data = 'recording_name';

option = 1;             % =0, raster plot of all clusters   (TO REVIEW)
                        % =1, raster plot of selected clusters, for selected time window + population firing rate   (TO REVIEW)
                        % =2, raster plot + burst identification (e.g., to test the burst detection algorithm)
                        % =3, raster plot (single cluster) + burst identification + shuffling (e.g., to test the shuffling algorithm)
                        % =4, raster plot + hypersynchronized units + shuffling     (TO REVIEW)
                        % =5, electrode map (color-coded mean firing rate)      (TO REVIEW)
                        % =6, waveforms     (TO REVIEW)
                        % =7, autocorrelation       (TO REVIEW)
                        % =8, raw raster plot       (TO REVIEW)

switch option
   
    case 0                      % raster plot of all clusters 

        option_save = 1;
        option_condition = 1;   % =1, show the label for each recording condition
        LABEL_POS = 0;

    case 1                      % raster plot of selected clusters, for selected time window + population firing rate

        option_explore = 0;
        channels = -1;
        range = [300 900];     % range, in sec
        option_id = 0;
        save_png = 1;

    case 2                      % raster plot + burst identification (e.g., to test the burst detection algorithm) 

        channels = 1:10;        % list of channels to be plotted (=-1, all channels)
        range = [100 600];      % time window for the plot, in sec
        n_colors = 5;           % numbe of different colors to label to different elements in the burst
        color{1} = [0 1 1];        color{2} = [1 0 1];        color{3} = [0 1 0];        color{4} = [0 0 1];        color{5} = [1 0 0];
        METHOD = 1;             % burst detection method (=1      Maximal Interval)

    case 3                      % raster plot (single cluster) + burst identification + shuffling (e.g., to test the shuffling algorithm)

        channel = 5;            % cluster to be analyzed and plotted (single cluster)
        range = [1 300];        % time window for the plot, in sec
        n_iteration = 5;        % total number of iterations for the shuffling
        METHOD = 1;             % burst detection method (=1      Maximal Interval)
        methodSF = 3;           % spike train shuffling method

    case 4                      % raster plot + hypersynchronized units + shuffling     

        channel = [113 114];
        range = [1 1200];     % range, in sec
        n_iteration = 5;
        option_method = 1;

    case 5                      % electrode map (color-coded mean firing rate)      

        channels = -1;
        range = [1 900];    % range, in sec    
        MAX_PLOT = 0.01;
        MATRIX_SIZE = 64;

        option_metric = 2;          % =1, total # units / recording length
                                    % =2, firing rate using median(isi)
                                    % =3, maximal projection    
        QUANT = 0.1;

        option_circle = 0;
        if  option_circle
            x0 = 25;
            y0 = 37;
            radiusx = 25;
            radiusy = 25;
        end

    case 6                      % waveforms     

        channels = 360:390;

    case 7                      % autocorrelation       

        channels = 30:35;

    case 8                      % raw raster plot       

        option_condition = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                        % factor to rescale from msec to sec
SHIFT = 0.1;                                % factor to move the rows in the raster plot
BIN = 0.20;                                 % bin size (0.30)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the data

file_name = sprintf('%s%s%s_QCmanual.mat', folder, delimiter, file_data);
%file_name = sprintf('%s%s%s_cleaned.mat', folder, delimiter, file_data);
load(file_name, 'spikeTimestamp','rec_conditions');

switch option
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0              %  raster plot of all clusters         
        
        figure(1)
        n_clusters = size(spikeTimestamp,2);
        for k=1:n_clusters
            data = spikeTimestamp(k);
            plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+k,'|',"MarkerEdgeColor",[0 0 0]);   hold on;    clear data
        end 
        ylim([0 1+k]);      xlabel('Time (sec)');   ylabel('Neuron');   title(sprintf('Raster Plot (%g neurons)', n_clusters));

        if option_condition
            n_conditions = size(rec_conditions,1);
            for cond=1:n_conditions
                if cond==1                                                              % assuming that the first condition is the baseline
                    aux_condition = rec_conditions{cond,1};
                    aux = split(aux_condition,".b");
                    condition = aux(1);                                                 % name of the condition (file main name)    
                    t_ini = rec_conditions{cond,2}(1);
                else
                    aux_condition0 = rec_conditions{1,1};
                    aux0 = split(aux_condition0,".b");
                    condition0 = aux0(1);                                               % main name
                    aux_condition = extractAfter(rec_conditions{cond,1},condition0);
                    aux = split(aux_condition,".b");
                    condition = aux(1);                                                 % name of the recording condition, without the main name
                    t_ini = rec_conditions{cond,2}(1);
                    plot([t_ini t_ini], [0 n_clusters], '--r','LineWidth',2);       hold on
                end
                text(t_ini, LABEL_POS, condition,'Color','red');        hold on
            end
        end

        if option_save
           file_output = sprintf('%s%s%s_ALLclusters.fig', folder, delimiter, file_data);
           saveas(gcf, file_output)
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1              %  raster plot of selected clusters, for selected time window + population firing rate

        n_clusters = size(spikeTimestamp,2);
        if channels == -1
            channels = 1:size(spikeTimestamp,2);
        end

        if option_explore 
          
            figure(1)
            for k=1:n_clusters
                data = spikeTimestamp(k);
                plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+k,'|',"MarkerEdgeColor",[0 0 0]);
                text(data.timeStamp(end)*K_MSEC_SEC, k+0.5,sprintf('# %g', k));     hold on;        clear data
            end 
            ylim([0 1+k]);      xlabel('Time (sec)');       ylabel('Neuron')
        else
            
            n_channel = numel(channels);
            figure(1)
            subplot(3,1,[1 2])
            for ch=1:n_channel
                electr = channels(ch);
                spikes = spikeTimestamp(electr).timeStamp;                     % time stamp for the spikes in electrode "electr", in msec
                aux1 = find(spikes>range(1)/K_MSEC_SEC);
                aux2 = find(spikes<range(2)/K_MSEC_SEC);
                aux = intersect(aux1, aux2);
                if numel(aux)>0
                    plot(spikes(aux)*K_MSEC_SEC, ones(numel(aux),1)+ch,'|',"MarkerEdgeColor",[0 0 0]);
                    if option_id
                        text(spikes(aux(end))*K_MSEC_SEC, ch+0.5,sprintf('# %g', electr));
                    end
                end
                hold on
                clear data
            end
            ylim([1 2+ch])
%            xlabel('Time (sec)', 'FontSize', 18)
            ylabel('Neuron', 'FontSize', 18)

            %% Calculate the population activity (entire recording)
    
            length = range(2)-range(1);                                     % length of the recording
            n_bins = ceil(length/(BIN));                                    % total number of bins to process the data
            spikes_bin = nan(n_bins,4);

            for i=1:n_bins                                                  % across all the bins
                aux_spikes = 0;
                spikes_bin(i,2) = range(1)/K_MSEC_SEC+(i-1)*BIN/K_MSEC_SEC; % left side of the bin i, in msec
                spikes_bin(i,3) = spikes_bin(i,2) + BIN/K_MSEC_SEC;         % right side of the bin i, in msec
                spikes_bin(i,4) = (spikes_bin(i,2)+spikes_bin(i,3))/2;      % midpoint of the bin i, in msec
                for ch=1:n_channel                                          % across all clusters
                    electr = channels(ch);
                    aux1 = find(spikeTimestamp(electr).timeStamp >= spikes_bin(i,2));
                    aux2 = find(spikeTimestamp(electr).timeStamp < spikes_bin(i,3));
                    aux = intersect(aux1, aux2);                            % index of all spikes of cluster ch within the bin i
                    if numel(aux) > 0
                        aux_spikes = aux_spikes + numel(aux);               % total number of spikes in the bin i (cumulative across all clusters)
                    end
                    clear aux aux1 aux2
                end
                spikes_bin(i,1) = aux_spikes/n_channel;                     % total number of spikes in bin i, per channel
                % spikes_bin(i,2) = spikes_bin(i,2) + 800/K_MSEC_SEC; % borrar
                clear aux_spikes
            end

            %figure(2)       % population firing
            subplot(3,1,3)
            bar(spikes_bin(:,2)*K_MSEC_SEC,spikes_bin(:,1)/BIN);
            xlabel('Time (sec)', 'FontSize', 18)
            ylabel('Population Firing Rate (Hz)', 'FontSize', 15)

            if save_png
                file_png = fullfile(folder, sprintf('%s_plot.png',file_data));
                exportgraphics(gcf, file_png);
        
            end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2              %  raster plot + burst identification
        
        if channels == -1
            channels = 1:size(spikeTimestamp,2);                            % all clusters will be plotted
        end
        n_channel = numel(channels);                                        % total number of clusters to be plotted
        
        figure(1)                                                           % initiate the figure #1
        list_electr=[];                                                     % vector to load the list of clusters (as characters) 
        for ch=1:n_channel                                                  % across all clusters
            electr = channels(ch);                                          % number of cluster "ch"
            spikes = spikeTimestamp(electr).timeStamp;                      % time stamp data for the cluster ch
            burst_plot = [];                                                % to load the number of burst
            aux1 = find(spikes>range(1)/K_MSEC_SEC);                        
            aux2 = find(spikes<range(2)/K_MSEC_SEC);
            aux = intersect(aux1, aux2);                                    % index of spikes within the time window

            if numel(aux)>0                                                             % there are spikes within the window of interest
                BurstDetected = BurstDetection(spikes(aux), METHOD);                    % run the burst detection algorithm
                for idx=1:numel(aux)                                                    % all spikes within the plot window
                   if BurstDetected.spikes(idx,3)==0                                    % the spike "idx" of the cluster "electr" is not within a burst
                      aux_color = [0 0 0];                                              % assign black as color (non-bursting)
                   else                                                                 % the spike "idx" of the cluster "electr" is withing a burst
                      aux_color = color{rem(BurstDetected.spikes(idx,3),n_colors)+1};   % assign one out of n_colors
                      burst_plot = [burst_plot BurstDetected.spikes(idx,3)];            % update the list of bursts
                   end
                   plot(BurstDetected.spikes(idx,1),ch,'|',"MarkerSize",10,"MarkerEdgeColor",aux_color);  hold on     % plot the spike "idx" with the appropiate color
               end
            end
        
            if numel(burst_plot)>0                                                              % there is at least one burst
               for k=1:numel(burst_plot)                                                        % across all the burst detected in the cluster "ch"
                   plot(BurstDetected.IndivBurst(burst_plot(k)).start, 0.3+ch, 'r*'); hold on   % label the start (first spike) of the burst k
                   plot(BurstDetected.IndivBurst(burst_plot(k)).last, 0.3+ch, 'b+'); hold on    % label the end (last spike) of the burst k
                   plot(BurstDetected.IndivBurst(burst_plot(k)).peak, 0.3+ch, 'ko'); hold on    % label the peak (smaller isi) of the burst k
                   plot(BurstDetected.IndivBurst(burst_plot(k)).median, 0.3+ch, 'mo'); hold on  % label the mid-point of the burst k
               end
            end
           list_electr = sprintf('%s #%g', list_electr, electr);            % update the list of clusters in the plot
        end

        title('Raster Plot',sprintf('Electrodes:%s', list_electr));         % Title of the plot, indicating the list of clusters in the plot
        xlabel('Time (msec)'); ylabel('Neuron'); ylim([0 1+n_channel]);
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3              % raster plot + burst identification + shuffling

        spikes = spikeTimestamp(channel).timeStamp;                         % time stamp data for the cluster 
        aux1 = find(spikes>range(1)/K_MSEC_SEC);                        
        aux2 = find(spikes<range(2)/K_MSEC_SEC);
        aux = intersect(aux1, aux2);                                        % index of spikes within the time window
        if numel(aux)>0                                                     % there is at least one spike 
            signal(:,1) = spikes(aux);                                      % time spike of all spikes within the window
            n_spikes = numel(signal);                                       % total number of spikes (the original number, before deleting the last spike)
            BurstDetected0 = BurstDetection(signal, METHOD);                % run the burst detection algorithm for the original signal (it only takes the first column)
            BurstFeature0 = BurstFeatures(BurstDetected0);                  % characterize the bursting activity, for the original signal
            n_burst=size(BurstDetected0.IndivBurst,2);                      % total number of burst in the original signal
  
            figure(1)                                                                   % initiate figure #1
            plot(signal*K_MSEC_SEC, ones(n_spikes,1),'|',"MarkerEdgeColor",[0 0 0]);    % plot the raster plot (original signal, in black)
            for i=1:n_burst                                                             % across all bursts
                hold on; plot(BurstDetected0.IndivBurst(i).start*K_MSEC_SEC, 1, '*r');  % label the start of the burst
                hold on; plot(BurstDetected0.IndivBurst(i).last*K_MSEC_SEC, 1, '*m');   % label the end of the burst
            end

            % to test the shuffling method
             
            % f_burst = nan(n_iteration+1,1);
            % f_nonburst = nan(n_iteration+1,1);
            % nSpikes_burst = nan(n_iteration+1,1);
            % nSpikes_total = nan(n_iteration+1,1);
            % nSpikes_nonburst = nan(n_iteration+1,1);
            % f_burst(1) = 1/(BurstFeature0.neuronMetric.isi.mean*K_MSEC_SEC);
            % f_nonburst(1) = 1/(BurstFeature0.neuronMetric.nonBursting.isi.mean*K_MSEC_SEC);
            % nSpikes_burst(1) = n_burst*BurstFeature0.neuronMetric.nSpikes.mean;
            % nSpikes_total(1) = n_spikes;
            % nSpikes_nonburst(1) = BurstFeature0.neuronMetric.nonBursting.nSpikes_norm*nSpikes_total(1);

            if n_iteration>0
               refractPeriod = 0;
               ShuffleTrain = ShuffleSpikeTrain(signal,range./K_MSEC_SEC,spikeTimestamp(channel).t_step,methodSF,n_iteration,refractPeriod,BurstDetected0,BurstFeature0);
               for k=1:n_iteration                                          % across all iterations
                   hold on; plot(ShuffleTrain{k}(:,1)*K_MSEC_SEC, ones(size(ShuffleTrain{k},1))+k*SHIFT,'|',"MarkerEdgeColor",[0.7 0.7 0.7]);
                   % BurstDetected = BurstDetection(ShuffleTrain{k}(:,1), 1);                     
                   % BurstFeature = BurstFeatures(BurstDetected);   
                   % f_burst(k+1) = 1/(BurstFeature.neuronMetric.isi.mean*K_MSEC_SEC);
                   % f_nonburst(k+1) = 1/(BurstFeature.neuronMetric.nonBursting.isi.mean*K_MSEC_SEC);
                   % nSpikes_burst(k+1) = BurstFeature.nBurst*BurstFeature.neuronMetric.nSpikes.mean;
                   % nSpikes_total(k+1) = size(ShuffleTrain{k}(:,1),1);
                   % nSpikes_nonburst(k+1) = BurstFeature.neuronMetric.nonBursting.nSpikes_norm*nSpikes_total(1);
                   % clear BurstDetected BurstFeature
               end
            end
            for i=1:n_burst                                                 % across all bursts (from the original data)
                hold on;  aux1 = BurstDetected0.IndivBurst(i).start*K_MSEC_SEC;
                plot([aux1; aux1], [1; 1+SHIFT*(n_iteration+1)], 'r:');     % draw a line at the start of the burst i
                hold on;  aux2 = BurstDetected0.IndivBurst(i).last*K_MSEC_SEC;
                plot([aux2; aux2], [1; 1+SHIFT*(n_iteration+1)], 'm:');     % draw a line at the end of the burst i
            end
            xlabel('Time (sec)');   ylim(range./K_MSEC_SEC);  ylim([0.5 1.5+SHIFT*n_iteration]);
        else
            error('There is no spikes in the selected time window')
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 4              % raster plot + hypersynchronized units + shuffling
        
        ini = nan(2,1);
        n_unit = nan(2,1);
        n_burst = nan(2,1);
        n_hs = nan(2+n_iteration,2);
        for ch=1:2
            aux_ini = find(spikeTimestamp(channel(ch)).timeStamp(:)<range(1)/K_MSEC_SEC);
            if numel(aux_ini)==0 || aux_ini<1
                ini(ch)=1;
            end
            aux_fin = find(spikeTimestamp(channel(ch)).timeStamp(:)>range(2)/K_MSEC_SEC);
            signal{ch} = spikeTimestamp(channel(ch)).timeStamp(ini(ch):aux_fin(1));
            signalSpike = signal{ch};
            signalSpike(end) = [];
            signalSpike(:,2) = diff(signal{ch});
            n_hs(ch,2) = numel(signal{ch});
            BurstDetected{ch} = BurstDetection(signalSpike, option_method);
            BurstFeature{ch} = BurstFeatures(BurstDetected{ch});
            n_burst(ch)=BurstDetected{ch}.nBurst;
                    
            burst{ch} = [BurstDetected{ch}.spikes(:,3); 0];
                        
            hold on
            clear signalSpike aux_fin aux_ini
        end
        signal1 = horzcat(signal{1},burst{1});
        signal2 = horzcat(signal{2},burst{2});
        HS_spikes0 = HyperSyncPair(signal1, signal2, t_step);
        
        aux0=find(HS_spikes0.data{1}(:,2)==0);
        aux1=find(HS_spikes0.data{1}(:,2)>0);
        figure(1)
        plot(HS_spikes0.data{1}(aux0,1)*K_MSEC_SEC,ones(numel(aux0),1),'|',"MarkerEdgeColor",[0 0 0]);
        hold on
        plot(HS_spikes0.data{1}(aux1,1)*K_MSEC_SEC,ones(numel(aux1),1),'*',"MarkerEdgeColor",[1 0 0]);
        hold on
        n_hs(1,1)=numel(aux1);
        clear aux0 aux1

        aux0=find(HS_spikes0.data{2}(:,2)==0);
        aux1=find(HS_spikes0.data{2}(:,2)>0);
        plot(HS_spikes0.data{2}(aux0,1)*K_MSEC_SEC,ones(numel(aux0),1)*1.1,'|',"MarkerEdgeColor",[0 0 0]);
        hold on
        plot(HS_spikes0.data{2}(aux1,1)*K_MSEC_SEC,ones(numel(aux1),1)*1.1,'*',"MarkerEdgeColor",[1 0 0]);
        hold on
        n_hs(2,1)=numel(aux1);
 
        xlabel('Time (sec)');
        ylim([0.5 1.5])

        % to test the hypersynchronicity on shuffled spike trains

        ShuffleTrain = ShuffleSpikeTrain(signal{2}, range./K_MSEC_SEC, t_step, 2, n_iteration, 0, BurstDetected{2}, BurstFeature{2});
        for it=1:n_iteration
            HS_spikes = HyperSyncPair(signal1, ShuffleTrain{it}, t_step);

            aux0=find(HS_spikes.data{2}(:,2)==0);
            aux1=find(HS_spikes.data{2}(:,2)>0);
            plot(HS_spikes.data{2}(aux0,1)*K_MSEC_SEC,ones(numel(aux0),1)+SHIFT*(it+1),'|',"MarkerEdgeColor",[0 0 0]);
            hold on
            plot(HS_spikes.data{2}(aux1,1)*K_MSEC_SEC,ones(numel(aux1),1)+SHIFT*(it+1),'*',"MarkerEdgeColor",[1 0 0]);
            hold on
            n_hs(it+2,1)=numel(aux1);
            n_hs(it+2,2)=size(ShuffleTrain{it},1);

        end
        xlabel('Time (sec)');
        ylim([0.5 1.7+n_iteration*SHIFT])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 5              %  electrode map (color-coded mean firing rate)        

        file_name = sprintf('%s%s%s_raw.mat', folder, delimiter, file_data);
        load(file_name, 'spikeTimestamp');
 
        if channels == -1
            channels = 1:size(spikeTimestamp,2);
        end
        n_channel = numel(channels);
        map = nan(MATRIX_SIZE,MATRIX_SIZE); 

        for i=1:n_channel
            pos_x = spikeTimestamp(i).positionXY(1)+1;
            pos_y = spikeTimestamp(i).positionXY(2)+1;

            spikes = spikeTimestamp(i).timeStamp;        
            aux1 = find(spikes>range(1)/K_MSEC_SEC);
            aux2 = find(spikes<range(2)/K_MSEC_SEC);
            aux = intersect(aux1, aux2);
            if numel(aux)>2
                isi = diff(spikes(aux));
            else 
                firing = 0;
            end

            switch option_metric
                case 1  % =1, total # units / recording length
                    firing = numel(aux) / (range(2)-range(1));
        
                case 2  % =2, firing rate using median(isi)
                    firing = 1/(median(isi)*K_MSEC_SEC);

                case 3  % =3, maximal projection    
                    firing = 1/(quantile(isi,QUANT)*K_MSEC_SEC);
            end

            map(pos_y,pos_x)=firing;
        end


       elipse = @(x,y,x0,y0,radius_x,radius_y) ((x-x0)/radius_x)^2 + ((y-y0)/radius_y)^2;
       map_plot = nan(MATRIX_SIZE,MATRIX_SIZE); 
       map_cropped = nan(MATRIX_SIZE,MATRIX_SIZE); 
       for i=1:MATRIX_SIZE 
           y_displayed = MATRIX_SIZE-i+1;                                  % because the heatmap starts the axis at the top-left; y=5 is y_displayed = 60
           for j=1:MATRIX_SIZE
               map_plot(i,j) = map(y_displayed,j); 
               map_cropped(i,j) = map_plot(i,j);
               if option_circle && elipse(j,y_displayed,x0,y0,radiusx*1.1,radiusy*1.1)>1 
                  map_cropped(i,j) = nan;
               end
           end
       end

       figure(1);
       heatmap(map_plot,'XDisplayLabels',NaN*ones(size(map,1),1),'YDisplayLabels',NaN*ones(size(map,1),1),'Colormap',parula,'ColorLimits',[0 MAX_PLOT]);

       if  option_circle
            pos = get(gca, 'Position');
            xapf = @(x,pos,xl) pos(3)*(x-min(xl))/diff(xl)+pos(1);                % 'x' Annotation Position Function
            yapf = @(y,pos,yl) pos(4)*(y-min(yl))/diff(yl)+pos(2);                % 'y' Annotation Position Function
            pos = gca().Position;
            xx=xlim;
            xl(1) = str2double(xx(1)); xl(2) = str2double(xx(2));
            yy=ylim;
            yl(1) = str2double(yy(1)); yl(2) = str2double(yy(2));
            annotation('ellipse', [xapf(x0-radiusx,pos,xl) yapf(y0-radiusy,pos,yl) radiusx*2*pos(3)/diff(xl) radiusy*2*pos(4)/diff(yl)],"Color","white","LineWidth",1, "LineStyle","--");
       end

       figure(2)
       heatmap(map_cropped,'XDisplayLabels',NaN*ones(size(map,1),1),'YDisplayLabels',NaN*ones(size(map,1),1),'Colormap',jet,'ColorLimits',[0 MAX_PLOT]);
       annotation('ellipse', [xapf(x0-radiusx,pos,xl) yapf(y0-radiusy,pos,yl) radiusx*2*pos(3)/diff(xl) radiusy*2*pos(4)/diff(yl)],"Color","white","LineWidth",1, "LineStyle","--");

       ellipse_cropped = [x0 y0 radiusx radiusy];
       % save(file_name,'ellipse_cropped','-append');   
       % figure(3)
       % heatmap(map,'XDisplayLabels',NaN*ones(size(map,1),1),'YDisplayLabels',NaN*ones(size(map,1),1),'Colormap',jet,'ColorLimits',[0 MAX_PLOT]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 6              %  waveforms     

        % file_name = sprintf('%s\\%s\\%s_QCmanual.mat', main_folder, project_folder, file_data);
        % load(file_name, 'spikeTimestamp');

        file_name = sprintf('%s\\%s\\waveforms.mat', main_folder, project_folder);
        load(file_name, 'waveforms');

        n_channel = numel(channels);

        for i=1:n_channel
            figure(i)
            plot(waveforms{channels(i)});
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 7              %  autocorrelation

        file_name = sprintf('%s\\%s\\%s_QCmanual.mat', main_folder, project_folder, file_data);
        load(file_name, 'spikeTimestamp');
        n_channel = numel(channels);

        for i=1:n_channel
            spikes = spikeTimestamp(channels(i)).timeStamp;
            isi = diff(spikes);
            figure(1)
            histogram(isi);
            [acf,lags] = autocorr(spikes);
            figure(2)
            plot(lags,acf);
            a=1;
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 8              %  raster plot of all clusters from raw data

        file_name = sprintf('%s\\%s\\%s_raw.mat', main_folder, project_folder, file_data);
        load(file_name, 'spikeTimestamp','rec_conditions');

        n_conditions = size(rec_conditions,1);

        figure(1)
        n_clusters = size(spikeTimestamp,2);
        for k=1:n_clusters
            data = spikeTimestamp(k);
            plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+k,'|',"MarkerEdgeColor",[0 0 0]);
            hold on
            clear data
        end 
        ylim([0 1+k])
        xlabel('Time (sec)')
        ylabel('Neuron')
        title(sprintf('Raster Plot (%g neurons)', n_clusters));


end
end