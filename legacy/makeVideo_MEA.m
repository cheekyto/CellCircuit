function makeVideo_MEA = makeVideo_MEA()

% function to make a video of MEA recordings, showing the firing rate on a
% display od the electrodes as well as the raster plot

clc; clear; close all; 
addpath(genpath('C:\Users\Yang\OneDrive - Harvard University\Desktop\Research Associate\Analysis\'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

option_edges = 1;
channels = -1;                  % =-1, all channels    

t_duration = 70;                % length of the data (s) you want to extract
BIN = 400;
GAP_SIZE = 0.04;
SCALE_UPPLIM = 30;

MATRIX_SIZE = 64;

option_number = 1;              % =1, single recording
                                % =2, paired recording

main_folder = '/path/to/main_folder';
project_folder = 'project_folder';
file_data = 'recording_name';
outputFile = sprintf('%s\\%s\\%s.gif', main_folder, project_folder, file_data);

switch option_number

    case 1      % single recording
        project_folder = 'project_folder';
        file_data = 'recording_name';
        t_start1 = 820;                 % time from which to start data extraction (s)

    case 2      % paired recording
        
        project_folder1 = 'project_folder';
        file_data1 = 'recording_name';
        t_start1 = 562;                 % time from which to start data extraction (s)

        project_folder2 = 'project_folder';
        file_data2 = 'recording_name';
        t_start2 = 520;                 % time from which to start data extraction (s)

        % project_folder1 = 'project_folder';
        % file_data1 = 'recording_name';
        % t_start1 = 208;                 % time from which to start data extraction (s)
        % 
        % project_folder2 = 'project_folder';
        % file_data2 = 'recording_name';
        % t_start2 = 520;                 % time from which to start data extraction (s)
end        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MATRIX_SIZE = 64;

K_SECtoMSEC = 1000;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Load data

switch option_number

    case 1      % single recording

        file_name = sprintf('%s\\%s\\%s_QCmanual.mat', main_folder, project_folder, file_data);
        load(file_name, 'spikeTimestamp');
        spikeTimestamp1 = spikeTimestamp;
        n_clusters1 = size(spikeTimestamp,2);
        if channels == -1
            channels1 = 1:size(spikeTimestamp,2);
        end

    case 2      % paired recording
        
        file_name = sprintf('%s\\%s\\%s_QCmanual.mat', main_folder, project_folder1, file_data1);
        load(file_name, 'spikeTimestamp');
        spikeTimestamp1 = spikeTimestamp;
        n_clusters1 = size(spikeTimestamp1,2);
        if channels == -1
            channels1 = 1:size(spikeTimestamp,2);
        end
        clear 'spikeTimestamp'

        file_name = sprintf('%s\\%s\\%s_QCmanual.mat', main_folder, project_folder2, file_data2);
        load(file_name, 'spikeTimestamp');
        spikeTimestamp2 = spikeTimestamp;
        n_clusters2 = size(spikeTimestamp2,2);
        if channels == -1
            channels2 = 1:size(spikeTimestamp,2);
        end

end        

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Bin the spike rate

nbin_total = round(t_duration*K_SECtoMSEC/BIN);
spiking1 = zeros(n_clusters1, nbin_total);

for i=1:nbin_total                                                          % across all bins 

    for ch=1:n_clusters1 
        ini = t_start1*K_SECtoMSEC + (i-1)*BIN;                                  % left side of the bin i, in msec
        fin = ini + BIN;                                                        % right side of the bin i, in msec
        aux1 = find(spikeTimestamp1(ch).timeStamp >= ini);
        aux2 = find(spikeTimestamp1(ch).timeStamp < fin);
        aux = intersect(aux1, aux2);                                        % index of all spikes of cluster ch within the bin i
        if numel(aux) > 0
           spiking1(ch,i) = numel(aux)/(BIN/K_SECtoMSEC);      
        end
        clear aux aux1 aux2       
    end
end
if option_number == 2
    spiking2 = zeros(n_clusters2, nbin_total);
    for ch=1:n_clusters2 
        ini = t_start2*K_SECtoMSEC + (i-1)*BIN;                                  % left side of the bin i, in msec
        fin = ini + BIN;                                                        % right side of the bin i, in msec
        aux1 = find(spikeTimestamp2(ch).timeStamp >= ini);
        aux2 = find(spikeTimestamp2(ch).timeStamp < fin);
        aux = intersect(aux1, aux2);                                        % index of all spikes of cluster ch within the bin i
        if numel(aux) > 0
           spiking2(ch,i) = numel(aux)/(BIN/K_SECtoMSEC);      
        end
        clear aux aux1 aux2       
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setting up

map1 = nan(MATRIX_SIZE,MATRIX_SIZE,nbin_total); 
for i=1:n_clusters1
    pos_x = spikeTimestamp1(i).positionXY(1)+1;
    pos_y = spikeTimestamp1(i).positionXY(2)+1;
    for t=1:nbin_total
        map1(pos_x,pos_y,t)=spiking1(i,t);
    end
end
map2 = nan(MATRIX_SIZE,MATRIX_SIZE,nbin_total); 
if option_number == 2
    for i=1:n_clusters2
        pos_x = spikeTimestamp2(i).positionXY(1)+1;
        pos_y = spikeTimestamp2(i).positionXY(2)+1;
        for t=1:nbin_total
            map2(pos_x,pos_y,t)=spiking2(i,t);
        end
    end
end

if option_edges
    gap = round(MATRIX_SIZE*GAP_SIZE);
    for x=1:gap
        map1(x,:,:) = nan;
        map2(x,:,:) = nan;
    end
    for x=(MATRIX_SIZE-gap):MATRIX_SIZE
        map1(x,:,:) = nan;
        map2(x,:,:) = nan;
    end
    for y=1:gap
        map1(:,y,:) = nan;
        map2(:,y,:) = nan;
    end
    for y=(MATRIX_SIZE-gap):MATRIX_SIZE
        map1(:,y,:) = nan;
        map2(:,y,:) = nan;
    end
end
% heatmap(map,'XDisplayLabels',NaN*ones(size(map,1),1),'YDisplayLabels',NaN*ones(size(map,1),1),'Colormap',jet,'ColorLimits',[0 MAX_PLOT]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

fig_num = figure('Position', [100 100 400 200]);

for index=1:nbin_total
 
    %%% Plot firing rate map
    switch option_number

        case 1      % single recording

            platform = nan(MATRIX_SIZE,MATRIX_SIZE);
            platform(:,:) = map1(:,:,index);
            imAlpha = ones(size(platform));
            imAlpha(isnan(platform)) = 0;
            im = imagesc(platform,[0 SCALE_UPPLIM]);
            im.AlphaData = imAlpha;
            set(gca,'color',0*[1 1 1]);
            colormap(jet);
            colorbar('color','white'); 
            axis off;axis image;
            set(gcf,'color','black')
            % text(20, 70, 'APM', 'Color', 'w');

        case 2      % paired recording

            t=tiledlayout(1, 2, 'tileSpacing', 'none','Padding','compact');

            ax = nexttile(1,[1 1]);
            platform = nan(MATRIX_SIZE,MATRIX_SIZE);
            platform(:,:) = map1(:,:,index);
            imAlpha = ones(size(platform));
            imAlpha(isnan(platform)) = 0;
            im = imagesc(platform,[0 SCALE_UPPLIM]);
            im.AlphaData = imAlpha;
            set(gca,'color',0*[1 1 1]);
            % colormap(jet);
            % colorbar('color','white'); 
            axis off;axis image;
            set(gcf,'color','black')
            text(20, 70, 'APM', 'Color', 'w');

            ax2 = nexttile(2,[1 1]);
            platform = nan(MATRIX_SIZE,MATRIX_SIZE);
            platform(:,:) = map2(:,:,index);
            imAlpha = ones(size(platform));
            imAlpha(isnan(platform)) = 0;
            im = imagesc(platform,[0 SCALE_UPPLIM]);
            im.AlphaData = imAlpha;
            set(gca,'color',0*[1 1 1]);
            colormap(jet);
            colorbar('color','white'); 
            axis off;axis image;
            set(gcf,'color','black')
            text(25, 70, 'CDM4', 'Color', 'w');

    end

    %%% Save the video
    drawnow;  
    
    time = index*BIN/K_SECtoMSEC;
    sgtitle(sprintf('Time: %.1f  sec', time),'color','red') 
    del=0.1;

    frame = getframe(fig_num);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    index     
    if index == 1
        imwrite(imind,cm,outputFile,'gif','Loopcount',inf,'DelayTime',del);
    else
        imwrite(imind,cm,outputFile,'gif','WriteMode','append','DelayTime',del);
    end
    
end

end

