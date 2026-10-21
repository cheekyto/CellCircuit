function evaluatePrimaryBursts = evaluatePrimaryBursts(PopulationActivityDetected)

N_BURST_MIN = 8;                % (default: 8)
MIN_DISTANCE = 4.0;             % (default: 4.0)
apply_parameters('evaluatePrimaryBursts');                                  % (revised October 2026) values of the parameter set of "run_pipeline", if any (see "pipeline_parameters")


nPeaks = PopulationActivityDetected.nBurst;
nPeaks0 = nPeaks; 
if PopulationActivityDetected.nBurst > N_BURST_MIN
    distance = nan(nPeaks-1,1);
    for i=2:nPeaks
        distance(i-1) = PopulationActivityDetected.IndivBurst(i).peak - PopulationActivityDetected.IndivBurst(i-1).peak;        
    end
    bimodality = test_bimodality(distance);
    distance1 = min(bimodality.gm2.mu);
    distance2 = max(bimodality.gm2.mu);
    
    if bimodality.final && distance2 > distance1*MIN_DISTANCE
        for i=nPeaks:-1:2
            if abs(distance(i-1)-distance1) < abs(distance(i-1)-distance2)
                PopulationActivityDetected.IndivBurst(i) = [];
                nPeaks = nPeaks - 1;
            end
        end
        PopulationActivityDetected.nBurst = nPeaks;
        if isfield(PopulationActivityDetected.IndivBurst, 'level')           % (revised October 2026) also update the number of level-1 events after discarding peaks;
            PopulationActivityDetected.nBurst_level1 = sum([PopulationActivityDetected.IndivBurst.level] == 1);  % before, it kept the old value, and the windows
        end                                                                  % between level-1 events got NaN rows (crash with opt_distance = 1)
        evaluatePrimaryBursts = PopulationActivityDetected;
    else
        evaluatePrimaryBursts = PopulationActivityDetected;
    end
    fprintf("\n *evaluatePrimaryBursts* has discarded %g peaks because of the distance \n", nPeaks0-nPeaks);
else
    evaluatePrimaryBursts = PopulationActivityDetected;
end

end