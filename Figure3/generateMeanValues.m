function [mean1, mean2] = generateMeanValues(structure, one, two, fieldname)
    %Use for r_out structures to pull out mean between different traces.
    %structure is structure, one and two are files in structure, and
    %fieldname is a string for the field to extract from
    N = numel(structure);
    mean1 = nan(N,1);
    mean2 = nan(N,1);
    for idx=1:N
        if isempty(two)
            mean1(idx) = mean(structure(idx).file(one).(fieldname));
        else
            mean1(idx) = mean(structure(idx).file(one).(fieldname));
            mean2(idx) = mean(structure(idx).file(two).(fieldname));
        end
    end
end