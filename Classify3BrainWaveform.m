function unitLabel = Classify3BrainWaveform(fileRaw, spikeTimestamp, unitLabel, newRun)

% Function to manually classify the units identify by 3Brain, based on their waveform
% (revised October 2026: corrected the description of "newRun"; the raster plot and the waveforms are drawn with a single plot call each;
%  the answer to the prompt is accepted in upper or lower case)
if 0
delimiter = '/';
folderPath = '/path/to/recording_folder';             % Folder where to find the original data (xxx.bxr files)
file_data_3Brain = 'recording_name_3Brain';
fileRaw = sprintf('%s%s%s_3Brain_raw.mat', folderPath, delimiter, file_data_3Brain);    
file_3Brain = sprintf('%s%s%s_3Brain_cleaned.mat', folderPath, delimiter, file_data_3Brain);
load(file_3Brain, 'spikeTimestamp');

unitLabel = nan(size(spikeTimestamp,2),2);
unitLabel(1,1) = 1; 
unitLabel(1,2) = 1; 
unitLabel(2,1) = 628;
newRun=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fileRaw:          path to the "xxx_raw.mat" file  
% spikeTimestamp:   data curated after "QualityControl_Manual" 
% unitLabel:        variable with cluster id and waveform-based classification (in case the experiment was partially run before) 
% newRun:           =1, start the classification from the first unit (first time running the experiment)
%                   =0, continue a partially analyzed file, from the first unit without classification (NaN in unitLabel(:,2))
%                   (revised October 2026: the description was inverted with respect to the code and to "QualityControl_Manual")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_MSEC_SEC = 1/1000;                        % factor to rescale from msec to sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   unitLabel(i,1):         % cluster id of the unit i
%            (i,2):         % waveform-based classification of the unit i
                            % =1, good
                            % =2, weak
                            % =-1, noise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fileRaw, 'Indiv_Cond','waveformLengths','samplingRates','dacFactors','offsetValues');

n_clusters = size(spikeTimestamp,2);
waveformLength = waveformLengths(1);
samplingRate = samplingRates(1);
x = double((1:waveformLength)) / samplingRate;

if size(unitLabel,1) ~= n_clusters
    error('\n The number of clusters in spikeTimestamp does not match the number of units in unitLabel!!!')
end
if newRun==1
    startCluster = 1;
else
    aux = find(isnan(unitLabel(:,2))>0);
    if numel(aux)>0
       startCluster = aux(1);
       if unitLabel(startCluster,1) ~= spikeTimestamp(startCluster).cluster
           error('\n The cluster id in spikeTimestamp does not match the cluster id in unitLabel!!!')
       end
    else
       startCluster = 1;
    end
end

for idx=startCluster:1:n_clusters
    if idx==1
        fprintf("\n The waveform-based classification is starting from the begining \n")
    elseif startCluster>1 && idx==startCluster
        fprintf(sprintf('\n The waveform-based classification is starting from the middle of the list (%gth unit) \n', idx));
    end

    fprintf(sprintf('\n Cluster number   %g', spikeTimestamp(idx).cluster));
    fprintf(sprintf('\n It has %g units \n', spikeTimestamp(idx).nUnits));

    cond = 1;
    dacFactor = dacFactors(cond);
    offsetValue = offsetValues(cond);
    traces = [];
    for cond = 1:size(Indiv_Cond,1)
        waveformData = Indiv_Cond{cond,2};
        traces = vertcat(traces, waveformData{spikeTimestamp(idx).cluster+1});
        clear waveformData
    end
    y = offsetValue + dacFactor * mean(traces,1);
    dy = std(offsetValue + dacFactor*traces,0,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig1=figure();
    subplot(2,2,[1,2])
    all_times = vertcat(spikeTimestamp.timeStamp);                          % raster plot of all units with a single plot call (revised October 2026: it was one call per unit, for every classified unit)
    all_rows = repelem((1:n_clusters)', cellfun(@numel, {spikeTimestamp.timeStamp})') + 1;     % row of each spike in the raster plot (same as "ones(data.nUnits,1)+k")
    plot(all_times*K_MSEC_SEC, all_rows,'|',"MarkerEdgeColor",[0 0 0]);
    hold on;        k = n_clusters;
    data = spikeTimestamp(idx);
    plot(data.timeStamp*K_MSEC_SEC, ones(data.nUnits,1)+idx,'|',"MarkerEdgeColor","r","MarkerSize",12);
    clear data;            ylim([0 1+k]);            xlabel('Time (sec)');            ylabel('Neuron')

    subplot(2,2,3)
    plot(x, (offsetValue + dacFactor * traces)');                           % all the waveforms of the unit with a single plot call (one line per waveform, as before)
    hold on

    subplot(2,2,4)
    fill([x';flipud(x')],[(y-dy)';flipud((y+dy)')],[0 1 1],'linestyle','none')
    line(x,y);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Manual classification of 3Brain units
    prompt = "Classification? (G=good / W=weak / N=noise / Q=quit(&save)):  ";
    label = input(prompt,"s"); 
    unitLabel(idx,1) = spikeTimestamp(idx).cluster;
            switch lower(label)                                     % (revised October 2026) upper case is also accepted (before, 'G' was saved as 0 and 'Q' did not quit)
                case 'g'
                    unitLabel(idx,2) = 1;
                case 'w'
                    unitLabel(idx,2) = 2;
                case 'n'
                    unitLabel(idx,2) = -1;
                case 'q'
                    break
                otherwise
                    unitLabel(idx,2) = 0;
            end
    close(fig1);           
end           

end