%Make Organized Data Structure
%Mouse is top level
%Each data matrix is [Days, Sessions, Time, Neurons]
%File info is saved within the structure for each mouse as well

%Find files for Before Blast
befbasedir = '/Volumes/eng_research_handata/Hua-an/Data/TBI/Ali/TBI_BEFORE BLASTING DATA ANALYSIS';
[~, beflist] = system(sprintf('find "%s" -type f -name "R_bins_circle_20170804.mat"', befbasedir));

befbreaks = find(beflist == char(10));
befcell_list = cell(1,numel(befbreaks));
for idx = 1:numel(befbreaks)
    if idx == 1
        befcell_list{idx} = beflist(1:(befbreaks(idx)-1));
    else
        befcell_list{idx} = beflist((befbreaks(idx-1)+1):(befbreaks(idx)-1));
    end
end

%Identify Types from cell lists for Before
%Day 1
CTLorder_befD1 = [9, 1, 3, 13, 15];  %8, 10, 11, Moona1, Moona2
TBIorder_befD1 = [7, 11, 5, 17]; %6, 9, 12, Seth5
%Day 2
CTLorder_befD2 = [10, 2, 4, 14, 16];  %8, 10, 11, Moona1, Moona2
TBIorder_befD2 = [8, 12, 6, 18]; %6, 9, 12, Seth5
% %All Cells
% CTLorder_bef = [9, 10, 1, 2, 3, 4, 13, 14, 15, 16];  %8, 10, 11, Moona1, Moona2
% TBIorder_bef = [7, 8, 11, 12, 5, 6, 17, 18]; %6, 9, 12, Seth5

%Find files for After Blast
basedir = '/Volumes/eng_research_handata/Hua-an/Data/TBI/Ali/TBI_AFTER BLASTING DATA ANALYSIS';
[~, list] = system(sprintf('find "%s" -type f -name "R_bins_circle_20170804.mat"', basedir));

breaks = find(list == char(10));
cell_list = cell(1,numel(breaks));
for idx = 1:numel(breaks)
    if idx == 1
        cell_list{idx} = list(1:breaks(idx)-1);
    else
        cell_list{idx} = list(breaks(idx-1)+1:breaks(idx)-1);
    end
end

%Identify Types from cell lists
%Day 1
CTLorder_aftD1 = [9, 1, 3, 13, 15];  %8, 10, 11, Moona1, Moona2
TBIorder_aftD1 = [7, 11, 5, 17]; %6, 9, 12, Seth5
%Day2
CTLorder_aftD2 = [10, 2, 4, 14, 16];  %8, 10, 11, Moona1, Moona2
TBIorder_aftD2 = [8, 12, 6, 18]; %6, 9, 12, Seth5
% %All Cells
% CTLorder_aft = [9, 10, 1, 2, 3, 4, 13, 14, 15, 16];  %8, 10, 11, Moona1, Moona2
% TBIorder_aft = [7, 8, 11, 12, 5, 6, 17, 18]; %6, 9, 12, Seth5


%Load Data & Put into structure of interest
Fs = 20; %20 Hz sampling
ctlMouse = struct();
expMouse = struct();
allfields = {'trace', 'dftrace', 'fnormtrace', 'bintrace2', 'bintrace3'};

for fld = 1:numel(allfields)
    outfield = allfields{fld};
    expcnt = 1;
    ctlcnt = 1;
    %Before Cells
    for idx = 1:numel(befcell_list)
        fn = befcell_list{idx};
        load(fn);
        TL = numel(r_out(1).file(1).(outfield)); %Total Length of Signal
        if ismember(idx, TBIorder_befD1) %Experimental Mice
            traceMatrix = zeros(TL,numel(r_out));
            for cel = 1:numel(r_out)
                traceMatrix(:,cel) = r_out(cel).file(1).(outfield);
            end
            expMouse(expcnt).day(1).session(1).(outfield) = traceMatrix;
            expMouse(expcnt).fnDay1Bef = befcell_list{idx};
        elseif ismember(idx, TBIorder_befD2)
            traceMatrix = zeros(TL,numel(r_out));
            for cel = 1:numel(r_out)
                traceMatrix(:,cel) = r_out(cel).file(1).(outfield);
            end
            expMouse(expcnt).day(2).session(1).(outfield) = traceMatrix;
            expMouse(expcnt).fnDay2Bef = befcell_list{idx};
            expcnt = expcnt+1; %Only iterate after Day 2 done
        elseif ismember(idx, CTLorder_befD1) %Experimental Mice
            traceMatrix = zeros(TL,numel(r_out));
            for cel = 1:numel(r_out)
                traceMatrix(:,cel) = r_out(cel).file(1).(outfield);
            end
            ctlMouse(ctlcnt).day(1).session(1).(outfield) = traceMatrix;
            ctlMouse(ctlcnt).fnDay1Bef = befcell_list{idx};
        elseif ismember(idx, CTLorder_befD2)
            traceMatrix = zeros(TL,numel(r_out));
            for cel = 1:numel(r_out)
                traceMatrix(:,cel) = r_out(cel).file(1).(outfield);
            end
            ctlMouse(ctlcnt).day(2).session(1).(outfield) = traceMatrix;
            ctlMouse(ctlcnt).fnDay2Bef = befcell_list{idx};
            ctlcnt = ctlcnt+1; %Only iterate after Day 2 done
        end
    end

    expcnt = 1;
    ctlcnt = 1;
    %After Cells
    for idx = 1:numel(cell_list)
        fn = cell_list{idx};
        load(fn);
        TL = numel(r_out(1).file(1).(outfield)); %Total Length of Signal
        if ismember(idx, TBIorder_aftD1) %Experimental Mice
            for sess = 1:numel(r_out(1).file)
                traceMatrix = zeros(TL,numel(r_out));
                for cel = 1:numel(r_out)
                    traceMatrix(:,cel) = r_out(cel).file(sess).(outfield);
                end
                expMouse(expcnt).day(1).session(1+sess).(outfield) = traceMatrix;
            end
            expMouse(expcnt).fnDay1Aft = cell_list{idx};
        elseif ismember(idx, TBIorder_aftD2)
            for sess = 1:numel(r_out(1).file)
                traceMatrix = zeros(TL,numel(r_out));
                for cel = 1:numel(r_out)
                    traceMatrix(:,cel) = r_out(cel).file(sess).(outfield);
                end
                expMouse(expcnt).day(2).session(1+sess).(outfield) = traceMatrix;
            end
            expMouse(expcnt).fnDay2Aft = cell_list{idx};
            expcnt = expcnt+1; %Only iterate after Day 2 done
        elseif ismember(idx, CTLorder_aftD1) %Experimental Mice
            for sess = 1:numel(r_out(1).file)
                traceMatrix = zeros(TL,numel(r_out));
                for cel = 1:numel(r_out)
                    traceMatrix(:,cel) = r_out(cel).file(sess).(outfield);
                end
                ctlMouse(ctlcnt).day(1).session(1+sess).(outfield) = traceMatrix;
            end
            ctlMouse(ctlcnt).fnDay1Aft = cell_list{idx};
        elseif ismember(idx, CTLorder_aftD2)
            for sess = 1:numel(r_out(1).file)
                traceMatrix = zeros(TL,numel(r_out));
                for cel = 1:numel(r_out)
                    traceMatrix(:,cel) = r_out(cel).file(sess).(outfield);
                end
                ctlMouse(ctlcnt).day(2).session(1+sess).(outfield) = traceMatrix;
            end
            ctlMouse(ctlcnt).fnDay2Aft = cell_list{idx};
            ctlcnt = ctlcnt+1; %Only iterate after Day 2 done
        end
    end
end
