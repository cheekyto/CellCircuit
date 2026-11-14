function [name, dose_c, file_c] = label_condition(rec_conditions, cond, base_cond, TREATMENT, RATIO_MATCH)

% label (name), dose and file name of the recording condition "cond", from the file names: "Baseline" for the baseline, otherwise the
% treatments of TREATMENT found in the part of its file name after the name of the baseline ("+CGP+TTX"; "unknown<rest of the name>" when
% none), and the doses in that part (e.g., "_20uM" -> "20uM"; several doses are written one after the other). These are the labels of the
% conditions in the results ("PerCondition(cond).condition" and ".dose"), by which "Compile_MEAproject" groups the recordings
% (revised November 2026: moved, unchanged, from the loop across conditions of "Analyze_spikeTrain" to a local function of it, then to this
% file, so that "plotMEArecording" names the conditions as the compiled tables do)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   rec_conditions:     file name of each recording condition (first column; e.g., "rec_conditions" of the "_QCmanual.mat" file)
%   cond:               index of the condition to label
%   base_cond:          index of the baseline condition (see "baselineCondition")
%   TREATMENT:          labels of the treatments (Analyze_spikeTrain.TREATMENT of the parameter set, e.g., {'CGP', 'Gabazine', 'TTX'})
%   RATIO_MATCH:        minimum matching between a treatment label and the file name (Analyze_spikeTrain.RATIO_MATCH; see "lavelVAR")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   name:               label of the condition ('Baseline', '+TTX', '+CGP+TTX', 'unknown<rest of the name>')
%   dose_c:             dose(s) of the treatment found in the file name ("" if none)
%   file_c:             file name of the condition, without the extension (.brw / .bxr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_treatment = size(TREATMENT,2);
dose_c = "";
aux_cluster = split(rec_conditions{cond,1},".b");                               % assuming the raw data file was a .brw or .bxr file
file_c = aux_cluster(1);                                                        % file name for the condition
pat = caseInsensitivePattern("_" + digitsPattern + ("u"|"m"|"n") + "M");
if cond==base_cond                                                      % the baseline (revised November 2026: identified by its file name; before, the first condition)
    name = 'Baseline';                                       % name of the condition
else
    aux_condition0 = rec_conditions{base_cond,1};                       % file name of the baseline (revised November 2026: before, of the first condition)
    aux0 = split(aux_condition0,".b");
    condition0 = aux0(1);                                               % main name for the experiment
    aux_condition = extractAfter(rec_conditions{cond,1},condition0);
    aux_label = split(aux_condition,".b");                              % specific label for the condition, contained in the file name
    condition_label = [];
    for k=1:n_treatment
        if contains(convertCharsToStrings(aux_label{1}),lavelVAR(TREATMENT{k},RATIO_MATCH),'IgnoreCase',true)
            condition_label = horzcat(condition_label,'+');
            condition_label = horzcat(condition_label,TREATMENT{k});
        end
    end
    if size(condition_label)==0
        name = horzcat('unknown', aux_label{1});
    else
        name = condition_label;
        aux = extractAfter(extract(aux_condition, pat),"_");
        if size(aux,1)>1
            dose_aux = [];
            for jd=1:size(aux,1)                                        % (revised October 2026: "jd" instead of "j", the index of the loop across conditions)
                dose_aux = strcat(dose_aux, string(aux{jd}));
            end
            dose_c = dose_aux; 
        else
            if size(aux,1)>0
                dose_c = string(aux);
            end
        end
    end
    clear condition_label 
end
end
