% add the x and y of background ROI

sz = [1024,1024];
frameNum = 0;
basedir = 'H:\Hua-an\Data\TBI\Ali'; %Base
cd(basedir);


[~, list] = system( 'dir /B /S Processed_ROIs*_2017_07_09*.mat'); %Find After Files
temp_result = textscan( list, '%s', 'delimiter', '\n' );
result_1 = temp_result{1};
[~, list] = system( 'dir /B /S Processed_ROIs*_2017_07_08*.mat'); %Find 10 Before Files
temp_result = textscan( list, '%s', 'delimiter', '\n' );
result_2 = temp_result{1};


filename_list = [result_1; result_2];

for filename_list_idx=1:size(filename_list,1)
    
    
%     cx=filename_list{filename_list_idx,2};cy=filename_list{filename_list_idx,3};ix=1024;iy=1024;r=7;
%     [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
%     c_mask=((x.^2+y.^2)<=r^2);
%     [b,~,~,~] = bwboundaries(c_mask);
    %circ = plot(b{1}(:,2),b{1}(:,1),'r');
    %%%%%%
%             button = questdlg('Would you like to keep this ROI?',...
%                 'Keep ROI','Yes','No','Yes');
%             
%             if strcmp(button,'Yes')
        %delete(circ)
        %plot(b{1}(:,2),b{1}(:,1),'g');
        c_mask = zeros(sz);
        c_mask(1) = 1;
        roiMaskRP = regionprops(c_mask,...
          'Centroid', 'BoundingBox','Area',...
          'Eccentricity', 'PixelIdxList','Perimeter');
        currRoi = RegionOfInterest(roiMaskRP);
        set(currRoi,'FrameSize',sz);
        %roiList{end+1} = currRoi;
%             else
%                 delete(circ)
%             end


    [pathstr, name, ext] = fileparts(filename_list{filename_list_idx,1});
    fullpathstr = fullfile(basedir,pathstr(3:end)); %Indexing removes X:\ from /S
    cd(pathstr);
    load(strcat(name,ext))
    R = cat(1,R,currRoi);
    
    save('R_w_bg_20170710','R');
    
    
end