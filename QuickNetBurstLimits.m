function QuickNetBurstLimits = QuickNetBurstLimits(pks, locs, proms, smoothedSpikes, START_LEVEL, END_LEVEL)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Identifity the starting and ending point of each potential network burst 

contador = numel(pks);                                                      % total number of peaks
idx_peaks = locs; 
start_Nburst = nan(contador,1);
end_Nburst = nan(contador,1);

for k=1:contador
        
    threshold_ini = pks(k)-(1-START_LEVEL)*proms(k);                        % threshold to identify the start of the burst event k    
    threshold_fin = pks(k)-(1-END_LEVEL)*proms(k);                          % threshold to identify the end of the burst event k
            
    if k==1                                                                 % if it is the first event
        aux = find(smoothedSpikes(1:idx_peaks(k),1)<threshold_ini);         % index of all bins below the threshold
        [~, idx_min] = min(smoothedSpikes(1:idx_peaks(k),1));               % lowest point in between the peak k and the beginning
        shift = 0;                                                          % to adjust the index
    else                                                                    % for the rest (not the first event)
        aux = find(smoothedSpikes(idx_peaks(k-1):idx_peaks(k),1)<threshold_ini);% index of all bins (between peakl k and k-1) below threshold
        [~, idx_min] = min(smoothedSpikes(idx_peaks(k-1):idx_peaks(k),1));  % lowest point in between the peak k and k-1 
        shift = idx_peaks(k-1)-1;                                           % to adjust the index
    end

    if numel(aux)>0                                                         % at least one point crossing the threshold level
        start_Nburst(k) = aux(end)+shift;                                   % first point below the threshold, to the left of the peak k
    else 
        start_Nburst(k) = idx_min+shift;                                    % lowest point to the left of the peak k
    end
    clear aux shift idx_min
        
    if k==contador                                                          % if it is the last event
        aux = find(smoothedSpikes(idx_peaks(k):end,1)<threshold_fin);       % index of all bins below the threshold    
        [~, idx_min] = min(smoothedSpikes(idx_peaks(k):numel(smoothedSpikes)));  % lowest point in between the peak k and the end
    else                                                                        % for the rest (not the last one)                      
        aux = find(smoothedSpikes(idx_peaks(k):idx_peaks(k+1),1)<threshold_fin); % index of all bins (between peak k and k+1) below threshold
        [~, idx_min] = min(smoothedSpikes(idx_peaks(k):idx_peaks(k+1)));    % lowest point in between the peak k and k+1
    end
    shift = idx_peaks(k)-1;
            
    if numel(aux)>0                                                         % at least one point crossing the threshold level
        end_Nburst(k) = aux(1)+shift;                                       % first point below the threshold, to the right of the peak
    else
        end_Nburst(k) = idx_min + shift;                                    % lowest point to the right of the peak k
    end
    clear aux shift idx_min    
end

QuickNetBurstLimits.end = end_Nburst;
QuickNetBurstLimits.start = start_Nburst;

end