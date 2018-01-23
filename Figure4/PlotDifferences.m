%% Loading
%load('allData.mat')
%% First Pass Plots
%Initialize Things
ctlBef = zeros(5,2);
ctlAft = zeros(5,5,2);
expBef = zeros(4,2);
expAft = zeros(4,5,2);

%Loops to get points
for idx = 1:numel(ctlMouse)
    for sess = 1:(numel(ctlMouse(1).day(1).session)-1)
        ctlBef(idx,1) = mean(mean(ctlMouse(idx).day(1).session(1).bintrace3));
        ctlAft(idx,sess,1) = mean(mean(ctlMouse(idx).day(1).session(sess+1).bintrace3));
    end
end
for idx = 1:numel(ctlMouse)
    for sess = 1:(numel(ctlMouse(1).day(2).session)-1)
        ctlBef(idx,2) = mean(mean(ctlMouse(idx).day(2).session(1).bintrace3));
        ctlAft(idx,sess,2) = mean(mean(ctlMouse(idx).day(2).session(sess+1).bintrace3));
    end
end
for idx = 1:numel(expMouse)
    for sess = 1:(numel(expMouse(1).day(1).session)-1)
        expBef(idx,1) = mean(mean(expMouse(idx).day(1).session(1).bintrace3));
        expAft(idx,sess,1) = mean(mean(expMouse(idx).day(1).session(sess+1).bintrace3));
    end
end
for idx = 1:numel(expMouse)
    for sess = 1:(numel(expMouse(1).day(2).session)-1)
        expBef(idx,2) = mean(mean(expMouse(idx).day(2).session(1).bintrace3));
        expAft(idx,sess,2) = mean(mean(expMouse(idx).day(2).session(sess+1).bintrace3));
    end
end

%Normalization Values
ctlBefNormD1 = [ctlBef(:,1),ctlBef(:,1),ctlBef(:,1),ctlBef(:,1),ctlBef(:,1)];
ctlBefNormD2 = [ctlBef(:,2),ctlBef(:,2),ctlBef(:,2),ctlBef(:,2),ctlBef(:,2)];
expBefNormD1 = [expBef(:,1),expBef(:,1),expBef(:,1),expBef(:,1),expBef(:,1)];
expBefNormD2 = [expBef(:,2),expBef(:,2),expBef(:,2),expBef(:,2),expBef(:,2)];

%Normalize
ctlMiceD1 = ctlAft(:,:,1) ./ ctlBefNormD1;
ctlMiceD2 = ctlAft(:,:,2) ./ ctlBefNormD2;
expMiceD1 = expAft(:,:,1) ./ expBefNormD1;
expMiceD2 = expAft(:,:,2) ./ expBefNormD2;

%All Points
ctlAllD1 = [ctlBef(:,1), ctlAft(:,:,1)];
ctlAllD2 = [ctlBef(:,2), ctlAft(:,:,2)];
expAllD1 = [expBef(:,1), expAft(:,:,1)];
expAllD2 = [expBef(:,2), expAft(:,:,2)];

% %Normalized Mean Firing Rate Plots
% figure(); plot([1:5],mean(ctlMiceD1),'b',[1:5], mean(expMiceD1),'r')
% hold on;
% plot([1:5],mean(ctlMiceD1)+std(ctlMiceD1),'--b',[1:5], mean(expMiceD1)+std(expMiceD1),'--r')
% plot([1:5],mean(ctlMiceD1)-std(ctlMiceD1),'--b',[1:5], mean(expMiceD1)-std(expMiceD1),'--r')
% title('Mean Firing Rate Normalized by Pre-Blast Day 1')
% 
% figure(); plot([1:5],mean(ctlMiceD2),'b',[1:5], mean(expMiceD2),'r')
% hold on;
% plot([1:5],mean(ctlMiceD2)+std(ctlMiceD2),'--b',[1:5], mean(expMiceD2)+std(expMiceD2),'--r')
% plot([1:5],mean(ctlMiceD2)-std(ctlMiceD2),'--b',[1:5], mean(expMiceD2)-std(expMiceD2),'--r')
% title('Mean Firing Rate Normalized by Pre-Blast Day 2')
% 
% %Total Mean Firing Rate Plots
% figure(); plot([1:6],mean(ctlAllD1),'b',[1:6], mean(expAllD1),'r')
% hold on;
% plot([1:6],mean(ctlAllD1)+std(ctlAllD1),'--b',[1:6], mean(expAllD1)+std(expAllD1),'--r')
% plot([1:6],mean(ctlAllD1)-std(ctlAllD1),'--b',[1:6], mean(expAllD1)-std(expAllD1),'--r')
% title('Mean Firing Rate Day 1')
% 
% figure(); plot([1:6],mean(ctlAllD2),'b',[1:6], mean(expAllD2),'r')
% hold on;
% plot([1:6],mean(ctlAllD2)+std(ctlAllD2),'--b',[1:6], mean(expAllD2)+std(expAllD2),'--r')
% plot([1:6],mean(ctlAllD2)-std(ctlAllD2),'--b',[1:6], mean(expAllD2)-std(expAllD2),'--r')
% title('Mean Firing Rate Day 2')

%Standard Deviation Calculations
%Initialize Things
ctlBefSD = zeros(5,2);
ctlAftSD = zeros(5,5,2);
expBefSD = zeros(4,2);
expAftSD = zeros(4,5,2);

%Loops to get points
for idx = 1:numel(ctlMouse)
    for sess = 1:(numel(ctlMouse(1).day(1).session)-1)
        ctlBefSD(idx,1) = mean(std(ctlMouse(idx).day(1).session(1).dftrace));
        ctlAftSD(idx,sess,1) = mean(std(ctlMouse(idx).day(1).session(sess+1).dftrace));
    end
end
for idx = 1:numel(ctlMouse)
    for sess = 1:(numel(ctlMouse(1).day(2).session)-1)
        ctlBefSD(idx,2) = mean(std(ctlMouse(idx).day(2).session(1).dftrace));
        ctlAftSD(idx,sess,2) = mean(std(ctlMouse(idx).day(2).session(sess+1).dftrace));
    end
end
for idx = 1:numel(expMouse)
    for sess = 1:(numel(expMouse(1).day(1).session)-1)
        expBefSD(idx,1) = mean(std(expMouse(idx).day(1).session(1).dftrace));
        expAftSD(idx,sess,1) = mean(std(expMouse(idx).day(1).session(sess+1).dftrace));
    end
end
for idx = 1:numel(expMouse)
    for sess = 1:(numel(expMouse(1).day(2).session)-1)
        expBefSD(idx,2) = mean(std(expMouse(idx).day(2).session(1).dftrace));
        expAftSD(idx,sess,2) = mean(std(expMouse(idx).day(2).session(sess+1).dftrace));
    end
end

%All Points
ctlAllD1SD = [ctlBefSD(:,1), ctlAftSD(:,:,1)];
ctlAllD2SD = [ctlBefSD(:,2), ctlAftSD(:,:,2)];
expAllD1SD = [expBefSD(:,1), expAftSD(:,:,1)];
expAllD2SD = [expBefSD(:,2), expAftSD(:,:,2)];

% %Total Mean Firing Rate Plots
% figure(); plot([1:6],mean(ctlAllD1SD),'b',[1:6], mean(expAllD1SD),'r')
% hold on;
% plot([1:6],mean(ctlAllD1SD)+std(ctlAllD1SD),'--b',[1:6], mean(expAllD1SD)+std(expAllD1SD),'--r')
% plot([1:6],mean(ctlAllD1SD)-std(ctlAllD1SD),'--b',[1:6], mean(expAllD1SD)-std(expAllD1SD),'--r')
% title('Standard Deviation Day 1')
% 
% figure(); plot([1:6],mean(ctlAllD2SD),'b',[1:6], mean(expAllD2SD),'r')
% hold on;
% plot([1:6],mean(ctlAllD2SD)+std(ctlAllD2SD),'--b',[1:6], mean(expAllD2SD)+std(expAllD2SD),'--r')
% plot([1:6],mean(ctlAllD2SD)-std(ctlAllD2SD),'--b',[1:6], mean(expAllD2SD)-std(expAllD2SD),'--r')
% title('Standard Deviation Day 2')

% Cell Number Comparisons
ctlCelNumsD1 = zeros(5,2);
ctlCelNumsD2 = zeros(5,2);
expCelNumsD1 = zeros(4,2);
expCelNumsD2 = zeros(4,2);
for idx=1:numel(ctlMouse)
    ctlCelNumsD1(idx,1) = size(ctlMouse(idx).day(1).session(1).dftrace,2);
    ctlCelNumsD1(idx,2) = size(ctlMouse(idx).day(1).session(2).dftrace,2);
    ctlCelNumsD2(idx,1) = size(ctlMouse(idx).day(2).session(1).dftrace,2);
    ctlCelNumsD2(idx,2) = size(ctlMouse(idx).day(2).session(2).dftrace,2);
end
for idx=1:numel(expMouse)
    expCelNumsD1(idx,1) = size(expMouse(idx).day(1).session(1).dftrace,2);
    expCelNumsD1(idx,2) = size(expMouse(idx).day(1).session(2).dftrace,2);
    expCelNumsD2(idx,1) = size(expMouse(idx).day(2).session(1).dftrace,2);
    expCelNumsD2(idx,2) = size(expMouse(idx).day(2).session(2).dftrace,2);
end

% % Subtype Response Comparisons (Need to have compareEleSupIncDec.m run with variables in workspace for one day)
% dnum = 1; %Day Number Selection
% expresults = cell(numel(expMouse),5,9);
% for ms = 1:numel(expMouse)
%     for sess = 2:numel(expMouse(ms).day(dnum).session)
%         expresults{ms,sess-1,1} = mean(expMouse(ms).day(dnum).session(sess).bintrace3(:,logical(experimental(ms).upinc)));
%         expresults{ms,sess-1,2} = mean(expMouse(ms).day(dnum).session(sess).bintrace3(:,logical(experimental(ms).updec)));
%         expresults{ms,sess-1,3} = mean(expMouse(ms).day(dnum).session(sess).bintrace3(:,logical(experimental(ms).upnone)));
%         expresults{ms,sess-1,4} = mean(expMouse(ms).day(dnum).session(sess).bintrace3(:,logical(experimental(ms).downinc)));
%         expresults{ms,sess-1,5} = mean(expMouse(ms).day(dnum).session(sess).bintrace3(:,logical(experimental(ms).downdec)));
%         expresults{ms,sess-1,6} = mean(expMouse(ms).day(dnum).session(sess).bintrace3(:,logical(experimental(ms).downnone)));
%         expresults{ms,sess-1,7} = mean(expMouse(ms).day(dnum).session(sess).bintrace3(:,logical(experimental(ms).coninc)));
%         expresults{ms,sess-1,8} = mean(expMouse(ms).day(dnum).session(sess).bintrace3(:,logical(experimental(ms).condec)));
%         expresults{ms,sess-1,9} = mean(expMouse(ms).day(dnum).session(sess).bintrace3(:,logical(experimental(ms).connone)));
%     end
% end
% exp_mresults = cellfun(@mean, expresults);
% exp_mresults(isnan(exp_mresults)) = 0;
% 
% ctlresults = cell(numel(ctlMouse),5,9);
% for ms = 1:numel(ctlMouse)
%     for sess = 2:numel(ctlMouse(ms).day(dnum).session)
%         ctlresults{ms,sess-1,1} = mean(ctlMouse(ms).day(dnum).session(sess).bintrace3(:,logical(control(ms).upinc)));
%         ctlresults{ms,sess-1,2} = mean(ctlMouse(ms).day(dnum).session(sess).bintrace3(:,logical(control(ms).updec)));
%         ctlresults{ms,sess-1,3} = mean(ctlMouse(ms).day(dnum).session(sess).bintrace3(:,logical(control(ms).upnone)));
%         ctlresults{ms,sess-1,4} = mean(ctlMouse(ms).day(dnum).session(sess).bintrace3(:,logical(control(ms).downinc)));
%         ctlresults{ms,sess-1,5} = mean(ctlMouse(ms).day(dnum).session(sess).bintrace3(:,logical(control(ms).downdec)));
%         ctlresults{ms,sess-1,6} = mean(ctlMouse(ms).day(dnum).session(sess).bintrace3(:,logical(control(ms).downnone)));
%         ctlresults{ms,sess-1,7} = mean(ctlMouse(ms).day(dnum).session(sess).bintrace3(:,logical(control(ms).coninc)));
%         ctlresults{ms,sess-1,8} = mean(ctlMouse(ms).day(dnum).session(sess).bintrace3(:,logical(control(ms).condec)));
%         ctlresults{ms,sess-1,9} = mean(ctlMouse(ms).day(dnum).session(sess).bintrace3(:,logical(control(ms).connone)));
%     end
% end
% ctl_mresults = cellfun(@mean, ctlresults);
% ctl_mresults(isnan(ctl_mresults)) = 0;

%% Event Calculation
expevents = cell(2, numel(expMouse),6); %day, mouse, session
for dy = 1:2
    for ms = 1:numel(expMouse)
        for sess = 1:numel(expMouse(ms).day(dy).session)
            for cel = 1:size(expMouse(ms).day(dy).session(sess).bintrace3,2)
                trace = expMouse(ms).day(dy).session(sess).bintrace3(:,cel);
                expevents{dy,ms,sess}(cel) = sum(findPulses(trace)>0);
            end
        end
    end
end
exp_mevents = cellfun(@mean, expevents);
expediff = cellfun(@minus,expevents(:,:,6),expevents(:,:,2),'UniformOutput',false);
expD1 = squeeze(exp_mevents(1,:,:));
expD2 = squeeze(exp_mevents(2,:,:));

ctlevents = cell(2, numel(ctlMouse),6); %day, mouse, session
for dy = 1:2
    for ms = 1:numel(ctlMouse)
        for sess = 1:numel(ctlMouse(ms).day(dy).session)
            for cel = 1:size(ctlMouse(ms).day(dy).session(sess).bintrace3,2)
                trace = ctlMouse(ms).day(dy).session(sess).bintrace3(:,cel);
                ctlevents{dy,ms,sess}(cel) = sum(findPulses(trace)>0);
            end
        end
    end
end
ctl_mevents = cellfun(@mean, ctlevents);
ctlediff = cellfun(@minus,ctlevents(:,:,6),ctlevents(:,:,2),'UniformOutput',false);
ctlD1 = squeeze(ctl_mevents(1,:,:));
ctlD2 = squeeze(ctl_mevents(2,:,:));

% %Total Mean Event Rate Plots
% figure(); plot([1:6],mean(ctlD1),'b',[1:6], mean(expD1),'r')
% hold on;
% plot([1:6],mean(ctlD1)+std(ctlD1),'--b',[1:6], mean(expD1)+std(expD1),'--r')
% plot([1:6],mean(ctlD1)-std(ctlD1),'--b',[1:6], mean(expD1)-std(expD1),'--r')
% title('Mean Event Rate Day 1')
% 
% figure(); plot([1:6],mean(ctlD2),'b',[1:6], mean(expD2),'r')
% hold on;
% plot([1:6],mean(ctlD2)+std(ctlD2),'--b',[1:6], mean(expD2)+std(expD2),'--r')
% plot([1:6],mean(ctlD2)-std(ctlD2),'--b',[1:6], mean(expD2)-std(expD2),'--r')
% title('Mean Event Rate Day 2')

%% Plots for Paper
%Default Colors
barColors = [0,0,1; 1,0,0]; %Blue Control, Red Blast
barNames = {'Pre-Blast', 'Post-Blast 1', 'Post-Blast 2', 'Post-Blast 3', 'Post-Blast 4', 'Post-Blast 5'};

%% Plot Event Rates over Time
eventRate_D1 = [mean(ctlD1); mean(expD1)];
lowerErr_D1 = [abs(mean(ctlD1)-quantile(ctlD1,.25)); abs(mean(expD1)-quantile(expD1,.25))];
upperErr_D1 = [abs(mean(ctlD1)-quantile(ctlD1,.75)); abs(mean(expD1)-quantile(expD1,.75))];
eventRate_D2 = [mean(ctlD2); mean(expD2)];
lowerErr_D2 = [abs(mean(ctlD2)-quantile(ctlD2,.25)); abs(mean(expD2)-quantile(expD2,.25))];
upperErr_D2 = [abs(mean(ctlD2)-quantile(ctlD2,.75)); abs(mean(expD2)-quantile(expD2,.75))];

[bt1, hb1, he1] = errorbar_groups(eventRate_D1, lowerErr_D1, upperErr_D1, ...
    'bar_colors', barColors, 'optional_errorbar_arguments',{'LineWidth', 2.5,'LineStyle','none'});
title('Day 1 Event Rate')
ylim([0 3])
[bt2, hb2, he2] = errorbar_groups(eventRate_D2, lowerErr_D2, upperErr_D2, ...
    'bar_colors', barColors, 'optional_errorbar_arguments',{'LineWidth', 2.5,'LineStyle','none'});
title('Day 2 Event Rate')
ylim([0 3])
%% Plot decreased activity cell event rates
% Event Calculation
activity1 = '\\engnas.ad.bu.edu\Research\eng_research_handata\Kyle\Code\TBI_Paper_Scripts\Figure3\activityResults_20170918_Day1.mat'; %'/Volumes/eng_research_handata/Kyle/Code/TBI_Paper_Scripts/Figure3/activityResults_20170918_Day1.mat';
load(activity1,'exppcts','ctlpcts')
expact{1} = exppcts;
ctlact{1} = ctlpcts;
activity2 = '\\engnas.ad.bu.edu\Research\eng_research_handata\Kyle\Code\TBI_Paper_Scripts\Figure3\activityResults_20170918_Day2.mat'; %'/Volumes/eng_research_handata/Kyle/Code/TBI_Paper_Scripts/Figure3/activityResults_20170918_Day2.mat'
load(activity2,'exppcts','ctlpcts')
expact{2} = exppcts;
ctlact{2} = ctlpcts;
expdecevents = cell(2, numel(expMouse),5); %day, mouse, session
for dy = 1:2
    for ms = 1:numel(expMouse)
        for sess = 2:numel(expMouse(ms).day(dy).session)
            decCells = find(expact{dy}(ms).dec==1);
            for cel = 1:numel(decCells)
                trace = expMouse(ms).day(dy).session(sess).bintrace3(:,decCells(cel));
                expdecevents{dy,ms,sess-1}(cel) = sum(findPulses(trace)>0);
            end
        end
    end
end
expdec_mevents = cellfun(@mean, expdecevents);
expdecD1 = squeeze(expdec_mevents(1,:,:));
expdecD2 = squeeze(expdec_mevents(2,:,:));

ctldecevents = cell(2, numel(ctlMouse),5); %day, mouse, session
for dy = 1:2
    for ms = 1:numel(ctlMouse)
        for sess = 2:numel(ctlMouse(ms).day(dy).session)
            decCells = find(ctlact{dy}(ms).dec==1);
            for cel = 1:numel(decCells)
                trace = ctlMouse(ms).day(dy).session(sess).bintrace3(:,decCells(cel));
                ctldecevents{dy,ms,sess-1}(cel) = sum(findPulses(trace)>0);
            end
        end
    end
end
ctldec_mevents = cellfun(@mean, ctldecevents);
ctldecD1 = squeeze(ctldec_mevents(1,:,:));
ctldecD2 = squeeze(ctldec_mevents(2,:,:));

eventRate_decD1 = [mean(ctldecD1); mean(expdecD1)];
lowerErr_decD1 = [abs(mean(ctldecD1)-quantile(ctldecD1,.25)); abs(mean(expdecD1)-quantile(expdecD1,.25))];
upperErr_decD1 = [abs(mean(ctldecD1)-quantile(ctldecD1,.75)); abs(mean(expdecD1)-quantile(expdecD1,.75))];
eventRate_decD2 = [mean(ctldecD2); mean(expdecD2)];
lowerErr_decD2 = [abs(mean(ctldecD2)-quantile(ctldecD2,.25)); abs(mean(expdecD2)-quantile(expdecD2,.25))];
upperErr_decD2 = [abs(mean(ctldecD2)-quantile(ctldecD2,.75)); abs(mean(expdecD2)-quantile(expdecD2,.75))];

[bt1, hb1, he1] = errorbar_groups(eventRate_decD1, lowerErr_decD1, upperErr_decD1, ...
    'bar_colors', barColors, 'optional_errorbar_arguments',{'LineWidth', 2.5,'LineStyle','none'});
title('Day 1 Decreased Event Rate')
[bt2, hb2, he2] = errorbar_groups(eventRate_decD2, lowerErr_decD2, upperErr_decD2, ...
    'bar_colors', barColors, 'optional_errorbar_arguments',{'LineWidth', 2.5,'LineStyle','none'});
title('Day 2 Decreased Event Rate')

%% Plot Normalized Baseline Changes
expMeans = cell(2, numel(expMouse),5); %day, mouse, session
for dy = 1:2
    for ms = 1:numel(expMouse)
        for sess = 2:numel(expMouse(ms).day(dy).session)
            expMeans{dy,ms,sess-1} = mean(expMouse(ms).day(dy).session(sess).fnormtrace);
        end
    end
end
exp_mBL = cellfun(@mean, expMeans);
expBLD1 = squeeze(exp_mBL(1,:,:));
expBLD2 = squeeze(exp_mBL(2,:,:));

ctlMeans = cell(2, numel(ctlMouse),5); %day, mouse, session
for dy = 1:2
    for ms = 1:numel(ctlMouse)
        for sess = 2:numel(ctlMouse(ms).day(dy).session)
            ctlMeans{dy,ms,sess-1} = mean(ctlMouse(ms).day(dy).session(sess).fnormtrace);
        end
    end
end
ctl_mBL = cellfun(@mean, ctlMeans);
ctlBLD1 = squeeze(ctl_mBL(1,:,:));
ctlBLD2 = squeeze(ctl_mBL(2,:,:));

meanBL_D1 = [mean(ctlBLD1); mean(expBLD1)];
lowerErrBL_D1 = [abs(mean(ctlBLD1)-quantile(ctlBLD1,.25)); abs(mean(expBLD1)-quantile(expBLD1,.25))];
upperErrBL_D1 = [abs(mean(ctlBLD1)-quantile(ctlBLD1,.75)); abs(mean(expBLD1)-quantile(expBLD1,.75))];
meanBL_D2 = [mean(ctlBLD2); mean(expBLD2)];
lowerErrBL_D2 = [abs(mean(ctlBLD2)-quantile(ctlBLD2,.25)); abs(mean(expBLD2)-quantile(expBLD2,.25))];
upperErrBL_D2 = [abs(mean(ctlBLD2)-quantile(ctlBLD2,.75)); abs(mean(expBLD2)-quantile(expBLD2,.75))];

[bt1, hb1, he1] = errorbar_groups(meanBL_D1, lowerErrBL_D1, upperErrBL_D1, ...
    'bar_colors', barColors, 'optional_errorbar_arguments',{'LineWidth', 2.5,'LineStyle','none'});
title('Day 1 Mean Baseline')
[bt2, hb2, he2] = errorbar_groups(meanBL_D2, lowerErrBL_D2, upperErrBL_D2, ...
    'bar_colors', barColors, 'optional_errorbar_arguments',{'LineWidth', 2.5,'LineStyle','none'});
title('Day 2 Mean Baseline')

%% Decreased Cell Baseline Changes
% Event Calculation
baseline1 = '\\engnas.ad.bu.edu\Research\eng_research_handata\Kyle\Code\TBI_Paper_Scripts\Figure3\sigpopulation_20170918_Day1.mat'; %'/Volumes/eng_research_handata/Kyle/Code/TBI_Paper_Scripts/Figure3/sigpopulation_20170918_Day1.mat';
load(baseline1,'expmeans','ctlmeans')
expBL{1} = expmeans;
ctlBL{1} = ctlmeans;
baseline2 = '\\engnas.ad.bu.edu\Research\eng_research_handata\Kyle\Code\TBI_Paper_Scripts\Figure3\sigpopulation_20170918_Day2.mat'; %'/Volumes/eng_research_handata/Kyle/Code/TBI_Paper_Scripts/Figure3/sigpopulation_20170918_Day2.mat';
load(baseline2,'expmeans','ctlmeans')
expBL{2} = expmeans;
ctlBL{2} = ctlmeans;
expdownmeans = cell(2, numel(expMouse),5); %day, mouse, session
for dy = 1:2
    for ms = 1:numel(expMouse)
        for sess = 2:numel(expMouse(ms).day(dy).session)
            downCells = find(expBL{dy}(ms).down==1);
            for cel = 1:numel(downCells)
                expdownmeans{dy,ms,sess-1}(cel) = mean(expMouse(ms).day(dy).session(sess).fnormtrace(:,downCells(cel)));
            end
        end
    end
end
expdown_mBL = cellfun(@mean, expdownmeans);
expdownBLD1 = squeeze(expdown_mBL(1,:,:));
expdownBLD2 = squeeze(expdown_mBL(2,:,:));

ctldownmeans = cell(2, numel(ctlMouse),5); %day, mouse, session
for dy = 1:2
    for ms = 1:numel(ctlMouse)
        for sess = 2:numel(ctlMouse(ms).day(dy).session)
            downCells = find(ctlBL{dy}(ms).down==1);
            for cel = 1:numel(downCells)
                ctldownmeans{dy,ms,sess-1}(cel) = mean(ctlMouse(ms).day(dy).session(sess).fnormtrace(:,downCells(cel)));
            end
        end
    end
end
ctldown_mBL = cellfun(@mean, ctldownmeans);
ctldownBLD1 = squeeze(ctldown_mBL(1,:,:));
ctldownBLD2 = squeeze(ctldown_mBL(2,:,:));

meandownBL_D1 = [nanmean(ctldownBLD1); mean(expdownBLD1)];
lowerErrdownBL_D1 = [abs(nanmean(ctldownBLD1)-quantile(ctldownBLD1,.25)); abs(mean(expdownBLD1)-quantile(expdownBLD1,.25))];
upperErrdownBL_D1 = [abs(nanmean(ctldownBLD1)-quantile(ctldownBLD1,.75)); abs(mean(expdownBLD1)-quantile(expdownBLD1,.75))];
meandownBL_D2 = [nanmean(ctldownBLD2); mean(expdownBLD2)];
lowerErrdownBL_D2 = [abs(nanmean(ctldownBLD2)-quantile(ctldownBLD2,.25)); abs(mean(expdownBLD2)-quantile(expdownBLD2,.25))];
upperErrdownBL_D2 = [abs(nanmean(ctldownBLD2)-quantile(ctldownBLD2,.75)); abs(mean(expdownBLD2)-quantile(expdownBLD2,.75))];

[bt1, hb1, he1] = errorbar_groups(meandownBL_D1, lowerErrdownBL_D1, upperErrdownBL_D1, ...
    'bar_colors', barColors, 'optional_errorbar_arguments',{'LineWidth', 2.5,'LineStyle','none'});
title('Day 1 Mean Baseline')
[bt2, hb2, he2] = errorbar_groups(meandownBL_D2, lowerErrdownBL_D2, upperErrdownBL_D2, ...
    'bar_colors', barColors, 'optional_errorbar_arguments',{'LineWidth', 2.5,'LineStyle','none'});
title('Day 2 Mean Baseline')

%% Difference in Baseline Values
expmdiff = cell(2, numel(expMouse)); %day, mouse, session
for dy = 1:2
    for ms = 1:numel(expMouse)
        for cel = 1:size(expMouse(ms).day(dy).session(2).trace,2)
            trace1 = expMouse(ms).day(dy).session(2).trace(:,cel);
            trace5 = expMouse(ms).day(dy).session(6).trace(:,cel);
            expmdiff{dy,ms}(cel) = mean(trace5)-mean(trace1);
        end
    end
end
expmdiffmean = cellfun(@mean,expmdiff);
expediffmean = cellfun(@mean,expediff);

ctlmdiff = cell(2, numel(ctlMouse)); %day, mouse, session
for dy = 1:2
    for ms = 1:numel(ctlMouse)
        for cel = 1:size(ctlMouse(ms).day(dy).session(2).trace,2)
            trace1 = ctlMouse(ms).day(dy).session(2).trace(:,cel);
            trace5 = ctlMouse(ms).day(dy).session(6).trace(:,cel);
            ctlmdiff{dy,ms}(cel) = mean(trace5)-mean(trace1);
        end
    end
end
ctlmdiffmean = cellfun(@mean,ctlmdiff);
ctlediffmean = cellfun(@mean,ctlediff);

%Plots
%Day 1
figure(); hold on;
for idx=1:numel(expmdiff(1,:))
    scatter(expmdiff{1,idx},expediff{1,idx},5,'or','filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
    scatter(ctlmdiff{1,idx},ctlediff{1,idx},5,'ob','filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
end
if (numel(ctlmdiff(1,:))-numel(expmdiff(1,:))) > 0
    scatter(ctlmdiff{1,idx+1},ctlediff{1,idx+1},5,'ob','filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
end
emdm1 = mean(expmdiffmean(1,:));
eedm1 = mean(expediffmean(1,:));
eeneg1 = abs(eedm1-quantile(expediffmean(1,:),.025));
eepos1 = abs(eedm1-quantile(expediffmean(1,:),.975));
emneg1 = abs(emdm1-quantile(expmdiffmean(1,:),.025));
empos1 = abs(emdm1-quantile(expmdiffmean(1,:),.975));
errorbar(emdm1, eedm1, eeneg1, eepos1, emneg1, empos1, 'ok', 'MarkerSize', 8, 'MarkerFaceColor','r', 'LineWidth', 2)
cmdm1 = mean(ctlmdiffmean(1,:));
cedm1 = mean(ctlediffmean(1,:));
ceneg1 = abs(cedm1-quantile(ctlediffmean(1,:),.025));
cepos1 = abs(cedm1-quantile(ctlediffmean(1,:),.975));
cmneg1 = abs(cmdm1-quantile(ctlmdiffmean(1,:),.025));
cmpos1 = abs(cmdm1-quantile(ctlmdiffmean(1,:),.975));
errorbar(cmdm1, cedm1, ceneg1, cepos1, cmneg1, cmpos1, 'ok', 'MarkerSize', 8, 'MarkerFaceColor','b', 'LineWidth', 2)
legend('Blast','Sham','Location','NorthWest')
xlabel('Difference in Baseline Mean')
ylabel('Difference in Mean Number of Firing Events')
title('Day 1')

%Fit linear model
%Experimental
x = [expmdiff{1,:}]; y = [expediff{1,:}];
expPolyFitD1 = polyfit(x,y,1);
expPolyValD1 = polyval(expPolyFitD1, x);
expYResidD1 = y - expPolyValD1;
expSSResidD1 = sum(expYResidD1.^2);
expSSTotalD1 = (length(y) - 1) * var(y);
expRSqD1 = 1 - expSSResidD1/expSSTotalD1;
%Control
x = [ctlmdiff{1,:}]; y = [ctlediff{1,:}];
ctlPolyFitD1 = polyfit(x,y,1);
ctlPolyValD1 = polyval(ctlPolyFitD1, x);
ctlYResidD1 = y - ctlPolyValD1;
ctlSSResidD1 = sum(ctlYResidD1.^2);
ctlSSTotalD1 = (length(y) - 1) * var(y);
ctlRSqD1 = 1 - ctlSSResidD1/ctlSSTotalD1;

%Day 2
figure(); hold on;
for idx=1:numel(expmdiff(2,:))
    scatter(expmdiff{2,idx},expediff{2,idx},5,'or','filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
    scatter(ctlmdiff{2,idx},ctlediff{2,idx},5,'ob','filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
end
if (numel(ctlmdiff(2,:))-numel(expmdiff(2,:))) > 0
    scatter(ctlmdiff{2,idx+1},ctlediff{2,idx+1},5,'ob','filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
end
emdm2 = mean(expmdiffmean(2,:));
eedm2 = mean(expediffmean(2,:));
eeneg2 = abs(eedm2-quantile(expediffmean(2,:),.025));
eepos2 = abs(eedm2-quantile(expediffmean(2,:),.975));
emneg2 = abs(emdm2-quantile(expmdiffmean(2,:),.025));
empos2 = abs(emdm2-quantile(expmdiffmean(2,:),.975));
errorbar(emdm2, eedm2, eeneg2, eepos2, emneg2, empos2, 'ok', 'MarkerSize', 8, 'MarkerFaceColor','r', 'LineWidth', 2)
cmdm2 = mean(ctlmdiffmean(2,:));
cedm2 = mean(ctlediffmean(2,:));
ceneg2 = abs(cedm2-quantile(ctlediffmean(2,:),.025));
cepos2 = abs(cedm2-quantile(ctlediffmean(2,:),.975));
cmneg2 = abs(cmdm2-quantile(ctlmdiffmean(2,:),.025));
cmpos2 = abs(cmdm2-quantile(ctlmdiffmean(2,:),.975));
errorbar(cmdm2, cedm2, ceneg2, cepos2, cmneg2, cmpos2, 'ok', 'MarkerSize', 8, 'MarkerFaceColor','b', 'LineWidth', 2)
legend('Blast','Sham','Location','NorthWest')
xlabel('Difference in Baseline Mean')
ylabel('Difference in Mean Number of Firing Events')
title('Day 2')

%Fit linear model
%Experimental
x = [expmdiff{2,:}]; y = [expediff{2,:}];
expPolyFitD2 = polyfit(x,y,1);
expPolyValD2 = polyval(expPolyFitD2, x);
expYResidD2 = y - expPolyValD2;
expSSResidD2 = sum(expYResidD2.^2);
expSSTotalD2 = (length(y) - 1) * var(y);
expRSqD2 = 1 - expSSResidD2/expSSTotalD2;
%Control
x = [ctlmdiff{2,:}]; y = [ctlediff{2,:}];
ctlPolyFitD2 = polyfit(x,y,1);
ctlPolyValD2 = polyval(ctlPolyFitD2, x);
ctlYResidD2 = y - ctlPolyValD2;
ctlSSResidD2 = sum(ctlYResidD2.^2);
ctlSSTotalD2 = (length(y) - 1) * var(y);
ctlRSqD2 = 1 - ctlSSResidD2/ctlSSTotalD2;

