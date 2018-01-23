function fullRUpdate_All_CircleROIs()
    %% Code to find .mat files for TBI data, loop through them and add something to the structure

    %Find files
    basedir = '/mnt/eng_research_handata/Hua-an/Data/TBI/Ali';
    inname = 'circleTraces_08032017.mat';
    outname = 'R_bins_circle_20170804.mat';
    cell_list = findNestedFiles(basedir, inname);

    % %Moving Average Filter
    % MAwindow=10;
    % bf = ones(1,MAwindow)/MAwindow;
    % af = 1;

    %Butterworth Filter
    [bf, af] = butter(6, 2/(20/2)); %Butterworth at 2 Hz cutoff for 20Hz Sampling

    for idx = 1:numel(cell_list)
        clear r_out
        [dir,~,~] = fileparts(cell_list{idx});
        load(cell_list{idx});
        for iter = 1:numel(r_out)
            %For Individual Files
            nfiles = numel(r_out(iter).file);
            f0 = mean(r_out(iter).file(end).trace);

            for step = 1:nfiles %Trace Normalization
                tr = r_out(iter).file(step).trace;
                r_out(iter).file(step).dftrace = dFfromOtherVid(tr, f0); %dF from Last Video
                r_out(iter).file(step).fnormtrace = FnormOtherVid(tr, f0); %Normalize F by last video
            end
            totstd = concatenatedStd(r_out(iter), 'dftrace'); %Total Standard Deviation across all videos

            for step = 1:nfiles %Binarization
                dftrace = r_out(iter).file(step).dftrace; %Df Trace
                strace = filtfilt(bf, af, dftrace); %Smoothed Trace
                threshTrace2 = binthresholdOtherStd(strace, 2, totstd);
                threshTrace3 = binthresholdOtherStd(strace, 3, totstd);
                [pkstarts2, pkends2] = findPeakStartsEnds(strace, threshTrace2);
                [pkstarts3, pkends3] = findPeakStartsEnds(strace, threshTrace3);
                r_out(iter).file(step).bintrace2 = generateBinoFromStartAndEnd(dftrace, pkstarts2, pkends2);
                r_out(iter).file(step).bintrace3 = generateBinoFromStartAndEnd(dftrace, pkstarts3, pkends3);
            end
            if strfind(cell_list{idx},'AFTER') %Permutation Statistics for Post Blasting
                tr1 = r_out(iter).file(1).fnormtrace;
                tr5 = r_out(iter).file(5).fnormtrace;
                r_out(iter).meandiff5_1 = mean(tr5) - mean(tr1);
            end
        end
        fprintf('Saving %s\n',cell_list{idx});
        save(fullfile(dir,outname),'r_out');
    end
end


%% Added External Functions

function dftrace = dFfromOtherVid(trace, f0)
    %Calculate dF/F using F0 from another trace
    dftrace = (detrend(trace)-mean(detrend(trace)))/f0;
end

function fnorm_trace = FnormOtherVid(trace, f0)
    %Normalize Fluorescence by F0 from another video
    fnorm_trace = trace/f0;
end

function combstd = concatenatedStd(structure, fieldname)
    %Calculate standard deviation of several traces contained in a single
    %structure.  fieldname is a string of the field to pull traces from
    numfiles = numel(structure.file);
    trlength = numel(structure.file(1).trace);
    combtrace = zeros(numfiles*trlength,1);
    start=1;
    for idx = 1:numfiles
        combtrace(start:start+trlength-1,1) = structure.file(idx).(fieldname);
        start=start+trlength;
    end
    combstd = std(combtrace);
end

function bintrace = binthresholdOtherStd(trace, stdThreshold, otherStd)
    %Make everything above stdThreshold*std(trace) into 1s
    bintrace = trace > stdThreshold*otherStd;
end

function [pkstarts, pkends] = findPeakStartsEnds(trace, bintrace)
    %Get peak start indices and peak maximums from a derivative trace
    difftrace = diff(trace); %Derivatives
    pktrace = zeros(1,numel(trace));
    pktrace(bintrace) = trace(bintrace);
    [~, pklocs] = findpeaks(pktrace);
    if isempty(pklocs) %If no peaks
        pkstarts=NaN;
        pkends=NaN;
    else
        pkstarts = nan(numel(pklocs),1);
        pkends = nan(numel(pklocs),1);
        for idx = 1:numel(pklocs) %Loop through found peaks and find where derivative goes negative for both
            start = pklocs(idx)-1; %Minus 1 for shifting indicies by 1 with diff
            before = find(difftrace(1:start)<0);
            after = find(difftrace(start:end)<0);
            if isempty(before) %Handle Endpoints
                before = NaN;
                after = NaN;
            elseif isempty(after)
                after=numel(difftrace);
            end
            pkstarts(idx) = before(end); 
            pkends(idx) = start+after(1);
        end
    end
end

function binoArray = generateBinoFromStartAndEnd(origTrace, highStartPeakIndices, highEndPeakIndices)
    binoArray = zeros(size(origTrace));
    for id = 1:length(highStartPeakIndices)
        try
            binoArray(highStartPeakIndices(id):highEndPeakIndices(id)) = 1;
        catch
            disp('caught');
        end
    end
    if length(binoArray) > length(origTrace)
        binoArray = binoArray(1:length(origTrace));
    end
end

