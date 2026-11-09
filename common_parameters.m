function C = common_parameters()

% THE FEW PARAMETERS THAT MAY NEED ADJUSTING (new, November 2026)
% The window "CellCircuit" shows these parameters first, with a plain description ("Common adjustments"); the other parameters of
% "pipeline_parameters.m" stay under "All parameters (advanced)". They are the values that were changed for particular recordings (very
% few units, slow events, weak or false small events, secondary events). This is the DEFAULT list (at most 15 parameters are shown):
% in the window, each user can move parameters from the full list to this list and back, and the window keeps that user's list (MATLAB
% preferences; the button "Default list" goes back to this one). To change the default list, add or remove a row: the name in
% "pipeline_parameters.m" (Function.PARAMETER), a short name for the window, and what it does. The values themselves are not set here:
% they come from the parameter file of the run.
%
% C: table with the columns key (Function.PARAMETER), topic, setting (short name), help (what it does) and display: '' (the value as in
%    the parameter file) or 'period' (a rate in events per second, shown and typed in the window as its period, 1/rate, in seconds;
%    e.g., Compile_MEAproject.MIN_NETBURST = 1/210 is shown as 210 s)

rows = {
% key                                               topic               setting                                                 help    display
'Analyze_spikeTrain.window',                        'Analysis',         'Analysis window (s)',                                  'Seconds analyzed at the end of each condition (e.g., 600 = the last 10 min); -1 = the whole condition.', ''
'Analyze_spikeTrain.MIN_UNITS',                     'Analysis',         'Minimum number of units',                              'Recordings with fewer good units are not analyzed; the Excel tables report them with 0 active neurons and 0 Hz.', ''
'QualityControl_General.ACTIVITY_THRESHOLD',        'Quality control',  'Minimum firing rate, first check (Hz)',                'Units firing less are discarded by the first quality control.', ''
'QualityControl_Manual.minimum_activity',           'Quality control',  'Minimum firing rate, final check (Hz)',                'Units firing less in every condition (except TTX) are discarded by the last quality control.', ''
'Analyze_spikeTrain.NETWORK_BURST_DETECT_METHOD',   'Network bursts',   'Detection method (3, 4 or 5)',                         '3: population firing rate. 4: fraction of units firing, with thresholds from shuffled data. 5 (default): method 3, plus the large events that only method 4 finds.', ''
'NetworkBurstDetection.method3.THRESHOLD',          'Network bursts',   'Threshold of the large events (methods 3 and 5)',     'In standard deviations of the population firing rate (default 5; it must be above 3.5).', ''
'NetworkBurstDetection.MIN_CHAN',                   'Network bursts',   'Minimum units taking part in network bursts',          'With fewer participating units, the recording has no network bursting (lower it only for recordings with very few units).', ''
'PopulationActivitySize.MIN_CHAN',                  'Network bursts',   'Minimum active units in one event',                    'Events in which fewer units fire are not counted.', ''
'PopulationActivitySize.MIN_MULTIPLE',              'Network bursts',   'Minimum units firing 3 or more spikes in one event',   'Events with fewer such units are not counted.', ''
'NetworkBurstDetection.RECRUIT_REL_MIN',            'Network bursts',   'Minimum recruitment of the small events (0 to 1)',     'Units firing at the peak of a small event, above its surroundings, as a share of a typical large event of the recording; higher removes weaker small events (methods 3 and 5).', ''
'PopulationActivityDetection.PERCENTILE_SMALL',     'Network bursts',   'Threshold of the small events (method 4, percentile)', 'Peak above this percentile of shuffled data (default 99); 100 = no small events.', ''
'NetworkBurstDetection.GAP_WITH_PEAK',              'Network bursts',   'Minimum share of large events with secondary events',  'Below this share (0 to 1) of the large events followed by a secondary event, the secondary events are not kept (they stay primary events).', ''
'cluster_neuronal_events.MIN_FAMILY_EVENTS',        'Families',         'Minimum events in a family',                           'Smaller groups of similar events belong to no family.', ''
'cluster_neuronal_events.MIN_SIMILARITY',           'Families',         'Minimum similarity within a family (0 to 1)',          'The events of a family are, on average, at least this similar; higher is stricter.', ''
'Compile_MEAproject.MIN_NETBURST',                  'Excel tables',     'Longest average time between events to measure a family (s)', 'A family (or all the events together) is measured only if its events come, on average, more often than once every this many seconds; otherwise only its number of events is reported. Default 210 s: 3 events in a 600-s window are measured. (The pipeline uses the rate, 1/seconds.)', 'period'
};
C = cell2table(rows, 'VariableNames', {'key', 'topic', 'setting', 'help', 'display'});
C = convertvars(C, 1:5, 'string');

end
