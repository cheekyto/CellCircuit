function DynamicPopActDetect = DynamicPopActDetect(units, spikeTimestamp, BurstDetected, range, bin_min, bin_max, bin_center, threshold, opt_distance, opt_binopti, opt_norm, opt_plot) 

% function to detect population activation in a spike train data, using dynamic bin size
% It calls the functions "PopulationActivityDetection" and "PopActivityPeak_Integration"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   structure array spikeTimestamp (its length is the number of recording channels), with the following fields
% 	.recLength:	length of the recording, in seconds.
% 	.nUnits:	total number of identified units in the channel
% 	.timeStamp:	time stamp for all units (nUnits) in the channel, in msec

%   units(j):   =0, if channel j should be ignored;     =1, if channel j should be incoporated for population activity
%   range:      range of the recording (range(1)=start; range(2)=end), in sec
%   bin_min:    minimum size of the temporal bin to calculate the population firing rate, in sec        
%   bin_max:    maximum size of the temporal bin to calculate the population firing rate, in sec
%   bin_center: optimal size of the temporal bin, obtained from a previous round of screening, in sec
%   threshold:  in terms of the standard deviation 
%   opt_distance:   option for "PopActivityPeak_Integration"
%   opt_binopti     =0, do not run the bin optimization for the secondary peaks; =1, run the bin optimization
%                   =value, run the detection with that bin value
%   opt_plot:   =0, no plot; =1: plot the data and the identified network bursts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DynamicPopActDetect.spikesBin(i,1):           total number of spikes in bin i
% DynamicPopActDetect.time(i):                  midpoint of the bin i, in sec 
% DynamicPopActDetect.nBurst:                   total number of network bursts 
% DynamicPopActDetect.IndivBurst(j).peak:       peak point of network burst j, in sec       
% DynamicPopActDetect.IndivBurst(j).start:      start point of network burst j, in sec
% DynamicPopActDetect.IndivBurst(j).last:       last point of network burst j, in sec
% DynamicPopActDetect.IndivBurst(j).peakIdx:    bin index of the peak point of network burst j
% DynamicPopActDetect.bin:                      size of the bin to produce the final output, in sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SCALING = 2.0;                                  % factor to scale up the bin size to run the dynamic search (default: 2.0)
SCALING_ADJ = 0.75;                             % factor to progressively change the bin size (default: 0.75)

N_SHORTLIST = 3;                                % number of pre-candidates to build the shortlist of optimal bin sizes 
CHANGE_THRESHOLD = 0.40;                        % threshold for the metric: change of number of burst events (default: 0.40)

OPTION_LEVEL = 1;                               % =1, just analyze the level1 peaks to determine the optimal bin size and the secondary peaks
                                                % =0, analyze all peaks (level0 and level1) to determine the optimal bin size
MAX_IT = 100;
apply_parameters('DynamicPopActDetect');                                    % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")
opt_test = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Population Firing rate and metrics

if isnan(bin_center)    % standard progression
    bin_size = bin_min;
    contador=0;
    while bin_size(end)<bin_max & contador < MAX_IT
        bin_size = vertcat(bin_size, bin_size(end)*SCALING);               
        contador=contador+1;
    end
    n_bin = numel(bin_size);
else                    % refined progression around bin_center
    bin_size = bin_center/SCALING/SCALING;
    contador = 0;
    while bin_size(end) < (bin_max+bin_center)/2 & contador < MAX_IT
        bin_size = vertcat(bin_size, bin_size(end)*SCALING*SCALING_ADJ);
        contador = contador + 1;
    end
    n_bin = numel(bin_size);
end

PopulationActivityDetected = cell(1,n_bin);
n_burst = nan(n_bin,1);
ibi = nan(n_bin,1);
%duration = nan(n_bin,1);           % it does not seem to add anything    
cv_ibi = nan(n_bin,1);
change_Nburst = nan(n_bin,1);       % change in the number of detected burst events, between consecutive bin sizes
for k=1:n_bin                       % across all bin sizes
    bin = bin_size(k);
    PopulationActivityDetected{k} = PopulationActivityDetection(units, spikeTimestamp, range, bin, threshold, 0, opt_norm, 0, 0);
    
    switch OPTION_LEVEL
        
        case 0  % all peaks (large and small)
                    
            n_burst(k) = PopulationActivityDetected{k}.nBurst;
            if n_burst(k)>1
                aux_ibi = nan(n_burst(k),1);
%               aux_duration = nan(n_burst(k),1);
%               aux_duration(1) = PopulationActivityDetected{k}.IndivBurst(1).last-PopulationActivityDetected{k}.IndivBurst(1).start;
                for i=1:n_burst(k)
                    aux_ibi(i) = PopulationActivityDetected{k}.IndivBurst(i).peak;
%                    aux_ibi(i)=PopulationActivityDetected{k}.IndivBurst(i).peak - PopulationActivityDetected{k}.IndivBurst(i-1).peak;
%                   aux_duration(i) = PopulationActivityDetected{k}.IndivBurst(i).last-PopulationActivityDetected{k}.IndivBurst(i).start;
                end
                ibi(k) = mean(diff(aux_ibi), 'omitnan');
                cv_ibi(k) = std(diff(aux_ibi), 'omitnan')/ibi(k);
%               duration(k) = mean(aux_duration, 'omitnan');
                clear aux_ibi aux_duration
            end

        case 1  % only level1 peaks (large peaks)
            
            n_burst(k) = PopulationActivityDetected{k}.nBurst_level1;
            n_burst_total = PopulationActivityDetected{k}.nBurst;
           if n_burst(k)>1
                aux_ibi = nan(n_burst(k),1);
                idx=0;
                for i=1:n_burst_total
                    if PopulationActivityDetected{k}.IndivBurst(i).level == 1
                        idx = idx + 1;
                        aux_ibi(idx) = PopulationActivityDetected{k}.IndivBurst(i).peak;
                    end
                end
                ibi(k) = mean(diff(aux_ibi), 'omitnan');
                cv_ibi(k) = std(diff(aux_ibi), 'omitnan')/ibi(k);
                clear aux_ibi aux_duration
            end
    end
    if k>1
       change_Nburst(k-1) = (n_burst(k) - n_burst(k-1))/n_burst(k-1); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Identify the ideal bin size

[~, sort_idx] = sort(abs(change_Nburst), 'ascend');
% (revised October 2026) the shortlist of low CV is built only with the bin sizes that have a CV (at least 2 level-1 events); before, 
% "sort(cv_ibi)" put the NaN values at the end and the first N_SHORTLIST were taken, so a bin size without CV could enter the shortlist 
% (when fewer than N_SHORTLIST bin sizes have a CV), and it crashed with fewer than N_SHORTLIST bin sizes. Same shortlist in any other case
valid_cv = find(~isnan(cv_ibi));                                            % bin sizes with at least 2 level-1 events (the CV of the IBI is defined)
[~, order_cv] = sort(cv_ibi(valid_cv), 'ascend');
cv_shortlist = valid_cv(order_cv(1:min(N_SHORTLIST, numel(valid_cv))));    % (up to) N_SHORTLIST bin sizes with the lowest CV
idx_best = nan;
for i=min(N_SHORTLIST,n_bin):-1:1 
    idx = sort_idx(i);
    if abs(change_Nburst(idx)) < CHANGE_THRESHOLD      % subthreshold change in the number of burst events 
        if ismember(idx, cv_shortlist)                  % low variability of the inter-burst-interval
            idx_best = sort_idx(i);
        end
    end    
end
if n_bin==1
    idx_best=1;
end
if ~isnan(idx_best)         % the optimal bin size was identified, then look for secondary peaks
%    opt_plot=1;
    DynamicPopActDetect = PopActivityPeak_Integration(units, spikeTimestamp, BurstDetected, range, bin_size(idx_best), threshold, opt_norm, opt_distance, opt_binopti, OPTION_LEVEL, opt_plot);    
    DynamicPopActDetect.bin = bin_size(idx_best);
else
    DynamicPopActDetect.nBurst = 0;
    DynamicPopActDetect.bin = nan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if opt_test && opt_plot

    if ~isnan(idx_best)
        if idx_best>1
            PopulationActivityDetected{k} = PopulationActivityDetection(units, spikeTimestamp, range, bin_size(idx_best-1), threshold, 0, 1,0, 1);
        end
        PopulationActivityDetected{k} = PopulationActivityDetection(units, spikeTimestamp, range, bin_size(idx_best), threshold, 0, 1, 0,1);
        title('Ideal bin size')
        if idx_best<numel(bin_size)
            PopulationActivityDetected{k} = PopulationActivityDetection(units, spikeTimestamp, range, bin_size(idx_best+1), threshold, 0, 1,0, 1);
        end
    else
        idx_best=1;
    end

    figure()

    subplot(2,2,1)
    plot(bin_size, n_burst,'--o');  hold on;
    plot(bin_size(idx_best), n_burst(idx_best), 'xr')
    xlabel('Bin size (sec)');       ylabel('# Burst');      yscale('log')

    subplot(2,2,2)
    plot(bin_size, change_Nburst,'--o');    hold on;
    plot(bin_size(idx_best), change_Nburst(idx_best), 'xr')
    xlabel('Bin size (sec)');       ylabel('Change in #burst')

    subplot(2,2,3)
    plot(bin_size, ibi,'--o');  hold on
    plot(bin_size(idx_best), ibi(idx_best), 'xr')
    xlabel('Bin size (sec)');       ylabel('I.B.I.')

    subplot(2,2,4)
    plot(bin_size, cv_ibi,'--o');       hold on;
    plot(bin_size(idx_best), cv_ibi(idx_best), 'xr')
    xlabel('Bin size (sec)');       ylabel('CV of IBI')
end

end