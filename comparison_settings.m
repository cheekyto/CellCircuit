% Settings of "CompareNeurons" (new, November 2026): copy this file, edit the copy, and run
%   CompareNeurons("<path of the copy>")
% Only TABLES and COMPARISONS must be given; every other setting can be left as it is (the defaults below) or removed (then the default is
% used too). The groups are those of the column "group" of the recordings tables, compiled by the pipeline in <project>/<group>/ and
% <project>/<group>_Compiled.mat (every "_Analyzed.mat" file of a group folder is an organoid of that group, as for "Compile_MEAproject")

% Recordings
C.TABLES = ["/path/to/recordings_A.xlsx", "/path/to/recordings_B.csv"];    % recordings tables of "run_pipeline" (one or more)
C.PROJECTS = [];                    % project folder of each table (the folder with the group folders and the compiled files); [] = the folder
                                    % of each table, as "run_pipeline"
C.CONDITION = "baseline";           % recording condition of the activity measures: "baseline" (identified by its file name, as the pipeline)
                                    % or the label of a condition (column "Condition" of the compiled tables)

% Axes of the neuron scores: always fitted again, on the organoids of these groups (equal weight per organoid; at least 20 organoids of at
% least 2 groups, otherwise nothing is run). Groups that only help the fitting can be listed here without being compared
C.FIT_GROUPS = [];                  % [] = every group of the tables

% Comparisons: only these are run; each one is a list of 2 or more groups (the first one is the reference)
C.COMPARISONS = {["groupA", "groupB"], ["groupA", "groupC", "groupD"]};

% Activity of each organoid (question: does the composition explain the activity?). One row per measure:
%   name | column of a compiled table ("<table>.<column>"; tables: Firing, Bursting, NETbursting0 (pooled), NETbursting1 (primary),
%          NETbursting2 (all clusters), cNEs) | transformation ("none", "log" = log10, "yes/no" = value > 0, with a logistic model) |
%   organoids ("all", or "bursting" = only the organoids with primary network bursts)
C.ACTIVITY = [
    "Bursting",                     "NETbursting1.Number_NETburst_Events",  "yes/no",   "all"
    "Network-burst rate",           "NETbursting1.Rate_NETburst",           "log",      "bursting"
    "Network-burst participation",  "NETbursting1.Participation_NETburst",  "none",     "bursting"
    "Network-burst duration",       "NETbursting1.Duration_NETburst",       "log",      "bursting"
    "Firing rate",                  "Firing.Firing_Rate",                   "log",      "all"
    ];

% Models of the activity: activity ~ group + composition + covariates (+ random intercepts), fitted with and without the composition
C.COMPOSITION = ["MeanWaveformScore1", "MeanWaveformScore2", "MeanFiringPatternScore"];    % mean scores of the organoid
C.COMPOSITION_MODE = "separate";    % "separate": one model per mean score; "joint": all of them in the same model
C.COVARIATES = "MeanNoiseSimilarity";   % built in: MeanNoiseSimilarity (mean similarity of the scored units to the QC noise shape; always added
                                    % with MeanWaveformScore2, user's choice), LogUnits (log10 of the number of units), PrimaryBursts; or a column
                                    % of a compiled table ("<table>.<column>"), or a column of the recordings tables (numbers, or text = categories)
C.RANDOM = [];                      % columns of the recordings tables used as random intercepts (e.g., "batch"); [] = none (ordinary regression)

% Roles of the units in the network bursts of the baseline (question: do the scores relate to who leads and who takes part?)
C.ROLES = ["onset", "participation"];   % onset: percentile of the onset rank among the robust participants (0 = first); participation:
                                    % fraction of the primary network bursts with a spike of the unit
C.ROLE_SCORES = ["waveformScore1", "waveformScore2", "firingPatternScore"];
C.ROLE_COVARIATES = ["noiseSimilarity", "logRateOutside"];  % similarity of the template to the QC noise shape; log10 of the firing rate
                                    % outside the network bursts (+0.01)
C.MIN_BURSTS = 5;                   % organoids with fewer primary network bursts at baseline are left out of the role models
C.MIN_ROLE_UNITS = 20;              % ... and those with fewer units with all the values

% Results: <OUTPUT>/<NAME>_<date and time>.xlsx and .mat, and the axes in <NAME>_<date and time>_axes.mat
C.OUTPUT = [];                      % [] = the folder of this settings file
C.NAME = "comparison";
