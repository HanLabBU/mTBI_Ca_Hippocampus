function r_out=extract_trace(r_in, GUIpick, rmfilt, rmgreen, isLabel)
    %GUI pick is a logical (0 or 1).  1 if you'd like to select files via
    %GUI, or 0 if it should autoselect *.tif files with m_ as a prefix in 
    %the pwd as the tifs to use.
    %rmgreen is a logical (0 or 1).  Select 1 if you want to remove all
    %.tif files that have the name "green" in them.
    
    if GUIpick
        [selected_files,selected_dir] = uigetfile('*.tif','MultiSelect','on');
    else
        selected_dir = pwd;
        if rmfilt
            filt_struct = dir('*f_*.tif'); %Remove Filtered Videos from extraction
            all_struct = dir('m_*.tif');
            selected_struct = all_struct;
            for idx = 1:numel(filt_struct)
                match = strcmp({selected_struct.name}, filt_struct(idx).name)==1;
                selected_struct(match) = [];
            end
        else
            selected_struct = dir('m_*.tif');
        end
        if rmgreen
            green_struct = dir('*green*.tif'); %Remove green channels from extraction
            for idx = 1:numel(green_struct)
                match = strcmp({selected_struct.name}, green_struct(idx).name)==1;
                selected_struct(match) = [];
            end
        end
        temp_cell = struct2cell(selected_struct);
        selected_files = temp_cell(1,:);
    end
    
    whole_tic = tic;
    
    if class(selected_files)=='char'
        file_list(1).name = fullfile(selected_dir,selected_files);
    else
        file_list = cell2struct(fullfile(selected_dir,selected_files),'name',1);
    end
    
    for file_idx=1:length(file_list)
        
        
        filename = file_list(file_idx).name;
        fprintf(['Processing ',filename,'....\n']);
        
        InfoImage = imfinfo(filename);
        NumberImages=length(InfoImage);

        f_matrix = zeros(InfoImage(1).Height,InfoImage(1).Width,NumberImages,'uint16');
        for i=1:NumberImages
            f_matrix(:,:,i) = imread(filename,'Index',i,'Info',InfoImage);
        end
        
        f_matrix = double(reshape(f_matrix,InfoImage(1).Height*InfoImage(1).Width,NumberImages));
        
        for roi_idx=1:numel(r_in)
            current_mask = zeros(1,InfoImage(1).Height*InfoImage(1).Width);
            try
                current_mask(r_in(roi_idx).pixel_idx) = 1;
                r_out(roi_idx).pixel_idx = r_in(roi_idx).pixel_idx;
            catch
                current_mask(r_in(roi_idx).PixelIdxList) = 1;
                r_out(roi_idx).pixel_idx = r_in(roi_idx).PixelIdxList;
            end
            image_mask = reshape(current_mask,InfoImage(1).Height,InfoImage(1).Width);
            %Find Centroid code from
            %(https://www.mathworks.com/matlabcentral/answers/322369-find-centroid-of-binary-image)
            [y, x] = ndgrid(1:InfoImage(1).Height, 1:InfoImage(1).Width);
            centroid = mean([x(logical(image_mask)), y(logical(image_mask))]);
            xcent = centroid(1);
            ycent = centroid(2);
            BG_radius = 50; %Number of Pixels for Local Background Radius.  Assume 1 neuron about 10 pixels in diameter
            BG_mask = (y - ycent).^2 + (x-xcent).^2 <= BG_radius.^2;
            onlyBG_mask = BG_mask - image_mask;
            BG_vect = reshape(onlyBG_mask, 1,InfoImage(1).Height*InfoImage(1).Width);
            r_out(roi_idx).BG_idx = find(BG_vect == 1);
            %Take Out Trace Values & Update r_out
            current_trace = (current_mask*f_matrix)/sum(current_mask);
            current_BG = (BG_vect*f_matrix)/sum(BG_vect);
            r_out(roi_idx).file(file_idx).filename = filename;
            r_out(roi_idx).file(file_idx).trace = current_trace;
            r_out(roi_idx).file(file_idx).BGtrace = current_BG;
          
            
            if file_idx==1
                r_out(roi_idx).trace = current_trace;
                r_out(roi_idx).BGtrace = current_BG;
            else
                r_out(roi_idx).trace = cat(2,r_out(roi_idx).trace,current_trace);
                r_out(roi_idx).BGtrace = cat(2,r_out(roi_idx).BGtrace,current_BG);
            end
            
            %Include Label
            if isLabel
                r_out(roi_idx).isLabel = r_in(roi_idx).isLabel;
            end
            
        end
        
    end
    
    for roi_idx=1:numel(r_in)
        r_out(roi_idx).color = rand(1,3);
    end
        
    fprintf(['Total loading time: ',num2str(round(toc(whole_tic),2)),' seconds.\n']);
    
end