%Script to run batch trace extraction for a lot of files

%File inputs
%Single File Set
[roifile, basedir ] = uigetfile('.mat', 'Select ROI Cell List');
% --------------------------------------------------------------
%Multiple File Sets
%Use this section instead of the above section to find a single ROI file
%name recursively within one higher directory and extract traces from all
%those subsequent files/sessions.
%basedir = 'Dir';
%roifile = 'circleROIs.mat';
% --------------------------------------------------------------
%Get File Lists
file_list = findNestedFiles(basedir, roifile);

%Output File Name
savename = ['circleTraces_', datestr(now,'mm_dd_yyyy'), '.mat'];


%Loop through files
for idx = 1:numel(file_list) %Connection shut down.
    %Navigate to folder
    filename = file_list{idx};
    [pathstr, name, ext] = fileparts(filename);
    cd(pathstr);
    
    fprintf(['Loading ',filename,'\n']);
    load(filename); %load ROI (Assumes CellList is name of ROI struct)
    
    r_out = extract_trace(CellList, 1, 1, 1, 0);
    
    save(savename, 'r_out');
end