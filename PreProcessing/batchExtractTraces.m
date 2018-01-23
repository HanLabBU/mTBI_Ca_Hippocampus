%Script to run batch trace extraction for a lot of files

%File inputs
basedir = 'Dir';
roifile = 'circleROIs.mat';
savename = 'circleTraces_Date.mat';
file_list = findNestedFiles(basedir, roifile);

%Loop through files
for idx = 1:numel(file_list) %Connection shut down.
    %Navigate to folder
    filename = file_list{idx};
    [pathstr, name, ext] = fileparts(filename);
    cd(pathstr);
    
    fprintf(['Loading ',filename,'\n']);
    load(filename); %load ROI (Assumes CellList is name of ROI struct)
    
    r_out = extract_trace(CellList, 0, 1, 1, 0);
    
    save(savename, 'r_out');
end