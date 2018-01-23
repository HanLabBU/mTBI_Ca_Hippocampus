function [r_all, r_out] = remove_overlap(r_in, overlap_threshold)

    whole_tic = tic;
    
    frame_per_video = 2047;
    total_number_r = numel(r_in);
    
    
    for r_idx = 1:total_number_r
        
        if r_idx==1
            
            try
                other_pixel = unique(cat(1,r_in(2:total_number_r).pixel_idx));
            catch
                other_pixel = unique(cat(1,r_in(2:total_number_r).PixelIdxList));
            end
            
        elseif r_idx==total_number_r
            
            try
                other_pixel = unique(cat(1,r_in(1:total_number_r-1).pixel_idx));
            catch
                other_pixel = unique(cat(1,r_in(1:total_number_r-1).PixelIdxList));
            end
            
        else
            
            try
                other_pixel = unique(cat(1,r_in([1:r_idx-1 r_idx+1:total_number_r]).pixel_idx));
            catch
                other_pixel = unique(cat(1,r_in([1:r_idx-1 r_idx+1:total_number_r]).PixelIdxList));
            end
            
        end
        
        try
            overlapped_pixel_count = sum(ismember(r_in(r_idx).pixel_idx,other_pixel));
        catch
            overlapped_pixel_count = sum(ismember(r_in(r_idx).PixelIdxList,other_pixel));
        end

        r_all(r_idx).area = numel(r_in(r_idx).pixel_idx);
        
        if r_idx<total_number_r
            r_all(r_idx).overlapped_ratio = overlapped_pixel_count/r_all(r_idx).area;
        else
            r_all(r_idx).overlapped_ratio = 0;
        end

        try
            r_all(r_idx).pixel_idx = r_in(r_idx).pixel_idx;
        catch
            r_all(r_idx).PixelIdxList = r_in(r_idx).PixelIdxList;
        end

        r_all(r_idx).color = r_in(r_idx).color;
        r_all(r_idx).trace = r_in(r_idx).trace;
        r_all(r_idx).file = r_in(r_idx).file;
%         r_all(r_idx).Trace_df = (r_all(r_idx).Trace-mean(r_all(r_idx).Trace))/mean(r_all(r_idx).Trace);
%         r_all(r_idx).Trace_df_100 = 100*(r_all(r_idx).Trace_df)/(max(r_all(r_idx).Trace_df));
%         
%         for video_idx=1:floor(numel(r_all(r_idx).Trace)/frame_per_video)
%             trace_idx = [(frame_per_video*(video_idx-1)+1):(frame_per_video*video_idx)];
%             r_all(r_idx).Trace_detrend(trace_idx) = detrend(r_all(r_idx).Trace(trace_idx));
%             r_all(r_idx).Trace_detrend_dfi(trace_idx) = (r_all(r_idx).Trace_detrend(trace_idx)-mean(r_all(r_idx).Trace_detrend(trace_idx)))/mean(r_all(r_idx).Trace_detrend(trace_idx));
%             r_all(r_idx).Trace_dfi(trace_idx) = (r_all(r_idx).Trace(trace_idx)-mean(r_all(r_idx).Trace(trace_idx)))/mean(r_all(r_idx).Trace(trace_idx));
%             r_all(r_idx).Trace_dfi_detrend(trace_idx) = detrend(r_all(r_idx).Trace_df(trace_idx));
%         end
        
        
    end
    
    r_out = r_all([r_all.overlapped_ratio]<=overlap_threshold);
    
    image_size = [1024 1024];
    
    all_roi_image(1).data = zeros(image_size);
    all_roi_image(2).data = zeros(image_size);
    all_roi_image(3).data = zeros(image_size);
    for i=1:numel(r_all)
        try
            all_roi_image(1).data(r_all(i).pixel_idx) = r_all(i).color(1);
            all_roi_image(2).data(r_all(i).pixel_idx) = r_all(i).color(2);
            all_roi_image(3).data(r_all(i).pixel_idx) = r_all(i).color(3);
        catch
            all_roi_image(1).data(r_all(i).PixelIdxList) = r_all(i).color(1);
            all_roi_image(2).data(r_all(i).PixelIdxList) = r_all(i).color(2);
            all_roi_image(3).data(r_all(i).PixelIdxList) = r_all(i).color(3);
        end
    end
    
%     figure
%     imagesc(cat(3,all_roi_image.data));
%     title('All ROIs');
%     axis image;
%     axis off;
%     
%     all_roi_image(1).data = zeros(image_size);
%     all_roi_image(2).data = zeros(image_size);
%     all_roi_image(3).data = zeros(image_size);
%     for i=1:numel(r_out)
%         all_roi_image(1).data(r_out(i).PixelIdxList) = r_out(i).Color(1);
%         all_roi_image(2).data(r_out(i).PixelIdxList) = r_out(i).Color(2);
%         all_roi_image(3).data(r_out(i).PixelIdxList) = r_out(i).Color(3);
%     end
%     figure
%     imagesc(cat(3,all_roi_image.data));
%     title('Less Overlapped ROIs');
%     axis image;
%     axis off;

    fprintf(['Total loading time: ',num2str(toc(whole_tic)),' seconds.\n']);

end