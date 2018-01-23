function [resample1, resample2] = randomBinaryResample(structure, one, two, fieldname, Nsamps)
    %Use for r_out structures to pull out mean between different traces.
    %structure is structure, one and two are files in structure, and
    %fieldname is a string for the field to extract from (Should be binary
    %string).  Nsamps is the number of samples to pull out per window.
    R = 1000; %Number of Resamples
    TL = numel(structure(1).file(one).(fieldname)); %Total Length of Signal
    EndSamp = TL - Nsamps;
    N = numel(structure);
    resample1 = nan(R,N);
    resample2 = nan(R,N);
    for idx=1:N
        if isempty(two)
            for iter=1:R
                start = randsample(EndSamp,2,true);
                resample1(iter,idx) = sum(structure(idx).file(one).(fieldname)(start(1):start(1)+Nsamps))/Nsamps;
            end
            resample2 = [];
        else
            for iter=1:R
                start = randsample(EndSamp,2,true);
                resample1(iter,idx) = sum(structure(idx).file(one).(fieldname)(start(1):start(1)+Nsamps))/Nsamps;
                resample2(iter,idx) = sum(structure(idx).file(two).(fieldname)(start(2):start(2)+Nsamps))/Nsamps;
            end
        end
    end
end