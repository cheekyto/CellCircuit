function [families, n_families] = primaryFamilies(label1)

% (new, November 2026) groups of primary events used for the per-family results ("Analyze_spikeTrain", "Compile_NetBurstResults",
% "Compile_MEAproject"), from the family of each primary event ("clustering.label(1:nBurst,1)", see "cluster_neuronal_events")
% The groups are the families (labels 1, 2, ...), followed by one group with all the primary events that belong to no family (label 0,
% the events of the groups with fewer than MIN_FAMILY_EVENTS events), if there are any (user's choice: the remaining events are counted
% together, in a single row of the per-family results, marked as "remaining events (no family)", see "Compile_NetBurstResults"). When no
% primary event belongs to a family, the only group is label 0, with all the primary events

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   label1(k):          family of the primary event k (0 = no family)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   families(i):        label of the group i: the families in ascending order, then 0 (the remaining events) if any event has no family;
%                       the primary events of the group i are those with label1 == families(i)
%   n_families:         number of real families of primary events (the group with label 0 is not counted)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

label1 = label1(:);
families = unique(label1(label1 > 0));
n_families = numel(families);
if any(label1 == 0)
    families = [families; 0];                                               % the remaining events (no family), after the families
end
end
