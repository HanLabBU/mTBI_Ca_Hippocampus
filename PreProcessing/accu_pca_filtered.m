function [result] = accu_pca_filtered(r_in, variance_threshold)
    
    whole_tic = tic;
    
    [bf, af] = butter(6, 2/(20/2)); %Butterworth at 0.1 Hz cutoff for 20Hz Sampling
    
    for file_idx=1:numel(r_in(1).file)
        
        temp_traces = [];
        for roi_idx=1:numel(r_in)
%             temp_traces = cat(1,temp_traces,filtfilt(bf,af,r_in(roi_idx).file(file_idx).trace));
            temp_trace = r_in(roi_idx).file(file_idx).trace;
            temp_trace = (temp_trace-mean(temp_trace))/mean(temp_trace);
            temp_traces = cat(1,temp_traces,temp_trace);
        end
        [coeff,score,latent,tsquared,explained,mu] = pca(temp_traces);
        accu_explained = cumsum(explained);
        for threshold_idx=1:numel(variance_threshold)
            result(threshold_idx,file_idx) = find(accu_explained>=variance_threshold(threshold_idx),1);
        end
        
    end


    fprintf(['Total loading time: ',num2str(toc(whole_tic)),' seconds.\n']);

end