function labelVAR = lavelVAR(label,ratio_match)

% Function to generate all the possible combinations by which the name "label" could have been spelled 

n_match = ceil(size(label,2)*ratio_match);                              % number of characters I am looking for
comb_length = size(label,2) - n_match;
n_combination = size(label,2)-comb_length+2;
label_array = strings(1,n_combination);
idx = 1;
label_array(idx) = label;
for i=1:(size(label,2)-comb_length+1)
    aux = label;
    aux(i:comb_length+i-1)=[];
    idx = idx+1;
    label_array(idx) = aux;
end
labelVAR = label_array;
end