%% Code to find .mat files for TBI data, loop through them and add something to the structure

%Find files
[status, list] = system( 'dir /B /S R_Date.mat' ); 
temp_result = textscan( list, '%s', 'delimiter', '\n' );
filename_list = temp_result{1};

for filename_idx=1:numel(filename_list)
    filename = filename_list{filename_idx};
    try
        clear r_out;
    end
    fprintf(['Loading file :',num2str(filename_idx),'/',num2str(numel(filename_list)),'\n']);
    load(filename);
    
    [pathstr,name,ext] = fileparts(filename);
    save_roi_filename = [pathstr,'\R_final_w_bg_20170522',ext];
    save_pca_filename = [pathstr,'\pca_20170522',ext];
    
     [r_all, r_out] = remove_overlap(r_out, 0.05);
    
    traces = cat(2,r_out.trace);
    stdvals = std(traces);
    [sortstds, sortinds] = sort(stdvals); %Use for plotting/finding set
    keepvals = stdvals > 3; %Remove Artifacts from shifts & low stds.
    keepvals(end) = 1;
    
    othertraces = traces(:,keepvals);
 
    for iter = 1:numel(r_out)
        r_out(iter).stdKeep = keepvals(iter);
        if keepvals(iter) == 1
            r_out(iter).stdKeepTrace = traces(:,iter);
        elseif keepvals(iter) == 0
            r_out(iter).stdKeepTrace = nan(size(traces,1),1);
        end
    end
    
     save(save_roi_filename,'r_out','-v6');
    
%     r_out = r_out([r_out.stdKeep]==1);
    [result] = accu_pca_filtered(r_out(1:end-1), [90 95 99]);

    
    save(save_pca_filename,'result','-v6');
end

