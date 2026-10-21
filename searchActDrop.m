function searchActDrop = searchActDrop(activity, time_binned, activityMatrix, threshold, noise, search_start)

% function to detect the first significant drop in the general level of activity,
% using a binned population firing rate data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   activity:       array with binned population firing rate
%   time_binned:    time array, in sec
%   threshold:      absolute value of the threshold for detection (average activity level) 
%   search_start:   temporal point to start the search (use it also as the size the moving averaging)  
%   noise:          std of the activity level

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   searchActDrop.time
%   searchActDrop.activity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MOViNG_STEP = 0.25;                             % steps to move the averaging window, in term of the size of the window
THRESHOLD_BIS = 10.0;                            % (default: 10.0)
apply_parameters('searchActDrop');                                          % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Search for the primary peaks

binSize = median(diff(time_binned));
window_bin = round(search_start/binSize);
movingstep_bin = round(search_start*MOViNG_STEP/binSize);
limit = ceil((time_binned(end)-search_start)/(search_start*MOViNG_STEP));

searchActDrop.activity = nan; 
searchActDrop.time = search_start;
flag = 0;

for i=0:limit
    ini = i*movingstep_bin+1;
    fin = min(ini + window_bin, numel(activity));
    auc = mean(activity(ini:fin),'omitnan');
    if auc < threshold
        
        searchActDrop.activity = auc/threshold; 
        searchActDrop.time = time_binned(fin);
        flag = 1;
        break
    end
end

if flag     % look for a pronounced peak after the drop in activity
    for i=1:size(activityMatrix,1)
        endGap = find(~isnan(activityMatrix(i,:)), 1, 'last');            % last bin with data in the gap i (revised October 2026: also the longest gap, which 
        if numel(endGap) > 0                                                % has no NaN and was skipped with "aux = find(isnan(...)); if numel(aux) > 0"; for 
                                                                            % the other gaps, the same bin as before: the one before the first NaN)
            halfGap = round(endGap/2);
            if halfGap > fin
                signal = activityMatrix(i, 1:halfGap);
                aux2 = find(signal > threshold + noise*THRESHOLD_BIS);
                if numel(aux2)>0
                    if aux2(end) > fin
                        location = aux2(end)+window_bin;
                        time2 = nan;
                        if numel(time_binned) > location 
                            time2 = time_binned(location);
                        else
                            time2 = time_binned(end) - time_binned(endGap-window_bin); 
                        end
                        searchActDrop.time = max(time2, searchActDrop.time);
                    end
                end
            end
        end
    end
end

end