function event = PopulationSilence(spikeTrain, bins, min_gap, option_plot)

% function to dectect periods of population silence

addpath(genpath('C:\Users\Yang\OneDrive - Harvard University\Desktop\Research Associate\Analysis\'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FACTOR_FILTER = 8;
N_PEAKS = 5;                                % (5.0)
THRESHOLD = 0.5;                            % (0.5)
K_MSEC_SEC = 1/1000;                        % factor to rescale from msec to sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare data

n_clusters = size(spikeTrain,2);
% t_step = spikeTrain(1).t_step;
range = spikeTrain(1).window ./ K_MSEC_SEC;
length = range(2)-range(1);                                                 % length of the recording
nSpikes_total = 0;
        
for ch=1:n_clusters
    nSpikes_total = nSpikes_total + spikeTrain(ch).nUnits;
end
meanRate = nSpikes_total/(length*K_MSEC_SEC)/n_clusters;         % mean spike rate, 1/sec

idx = 0;
n_sweeps = numel(bins);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bin the spike train
for k=1:n_sweeps
    bin = bins(k);
    idx = idx+1;
    n_bins = ceil(length*K_MSEC_SEC/(bin));                                 % total number of bins to process the data
    spikes = zeros(n_bins,n_clusters);
    
    for i=1:n_bins                                                          % across all the bins
        ini = range(1)+(i-1)*bin/K_MSEC_SEC;                                % left side of the bin i, in msec
        fin = ini + bin/K_MSEC_SEC;                                         % right side of the bin i, in msec
        for ch=1:n_clusters                              % across all clusters
            aux1 = find(spikeTrain(ch).timeStamp >= ini);
            aux2 = find(spikeTrain(ch).timeStamp < fin);
            aux = intersect(aux1, aux2);                            % index of all spikes of cluster ch within the bin i
            if numel(aux) > 0
               spikes(i,ch) = 1;                                   % total number of spikes in the bin i (cumulative across all clusters)
            elseif numel(aux) > 1
                   spikes(i,ch) = 2; 
            end
            clear aux aux1 aux2
        end
     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify silent bins
        
     silentbin = zeros(n_bins,2);
     for i=1:n_bins
         aux = find(spikes(i,:)==0);
         silentbin(i,2) = numel(aux)/n_clusters;
         silentbin(i,1) = (i-0.5)*bin;
         clear aux
     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify population silence
           
    silent_prob = exp(-1*meanRate*bin);
    expected_silentCh = silent_prob;

    filtered_silent = FiltroPasabajo(silentbin(:,2),n_bins,bin, 1/(bin*FACTOR_FILTER));
    filtered_time = FiltroPasabajo(silentbin(:,1),n_bins,bin, 1/(bin*FACTOR_FILTER));

    aux = sort(filtered_silent,'descend');
    peak_mean = mean(aux(1:N_PEAKS));
    threshold = expected_silentCh + (peak_mean-expected_silentCh)*THRESHOLD;
    aux_peak = find(filtered_silent>threshold);                     % index of bins with more silents neurons than threshold
    aux_diff = diff(aux_peak);                                      % distance (in number of bins) between peaks
    aux_gap = find(aux_diff>round(min_gap/bin));                    % index of bins that face a gap until the following supra-threshold bin
    if numel(aux_peak)>0                                            % if there is at least one supra-threshold bin
       peak = nan(numel(aux_gap)+1,1);
       edges = nan(numel(aux_gap)+1,2);
  
       for j=1:numel(aux_gap)+1                                   % across all the peaks that sufficiently separated one from each other
           
           if j==1                                                 % if it is the first event
              ini = aux_peak(1);
           else                                                    % for the rest
              ini = aux_peak(aux_gap(j-1)+1);
           end
           
           if j==numel(aux_gap)+1                                  % if it is the last event
              fin = aux_peak(end);
           else                                                    % for the rest
              fin = aux_peak(aux_gap(j));
           end
           edges(j,1) = ini;                                       % index for the first bin of event k
           edges(j,2) = fin;                                       % index for the last bin of event k
           [~, peak(j)] = max(filtered_silent(ini:fin));  
           peak(j) = ini+peak(j)-1;                                % index for the bin with the peak/max # silent clusters in event k 
        end
        n_events= numel(peak);                                       % total number of potential wave events
    else                                                                % there is no network burst
        n_events = 0;
        peak=[];
    end
    event(idx).bin = bin;
    event(idx).peakTime = filtered_time(peak);
    event(idx).nTotal = n_events;
    event(idx).width = max(edges(:,2)-edges(:,1))*bin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the results

if option_plot
                
    figure(k)
    bar(silentbin(:,1),silentbin(:,2))
    hold on
    plot([silentbin(1,1) silentbin(end,1)], [expected_silentCh expected_silentCh], ':r','LineWidth',2)
    hold on
    plot(filtered_time,filtered_silent,'m');
    hold on
    plot_peak = nan(n_events,2);
    for i=1:n_events
        plot_peak(i,1) = filtered_time(peak(i));
        plot_peak(i,2) = filtered_silent(peak(i))+threshold-expected_silentCh/2;
    end
    plot(plot_peak(:,1), plot_peak(:,2),'g*');
    legend('Binned signal', 'Level Expected by chance', 'Filtered signal', 'Suprathreshold events')
    xlabel('Time (sec)')
    ylabel('Proportion of silent clusters')
    title(sprintf('Bin size: %g sec (%g events, width = %g sec)', bin, n_events, event(idx).width))
end

end
end