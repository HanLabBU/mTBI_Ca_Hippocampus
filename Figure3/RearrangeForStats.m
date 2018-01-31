%Get data into plotting/stats realm.

%Default Colors
barColors = [0,0,1; 1,0,0]; %Blue Control, Red Blast
posmat = [1.7306    0.4850    1.3904    1.1616] * 1e3;
%% Up/Down/Consistent Cells
fn{1} = 'sigpopulation_MMDDYYYY_Day1.mat';
fn{2} = 'sigpopulation_MMDDYYYY_Day2.mat';
expStats = struct();
ctlStats = struct();

for f = 1:numel(fn)
    load(fn{f}, 'expmeans', 'ctlmeans')
    for eidx=1:numel(expmeans)
        expStats(f).pctdown(eidx,1) = expmeans(eidx).pctdown;
        expStats(f).pctup(eidx,1) = expmeans(eidx).pctup;
        expStats(f).pctcon(eidx,1) = (expmeans(eidx).tot-expmeans(eidx).totup-expmeans(eidx).totdown)/expmeans(eidx).tot;
    end
    for cidx=1:numel(ctlmeans)
        ctlStats(f).pctdown(cidx,1) = ctlmeans(cidx).pctdown;
        ctlStats(f).pctup(cidx,1) = ctlmeans(cidx).pctup;
        ctlStats(f).pctcon(cidx,1) = (ctlmeans(cidx).tot-ctlmeans(cidx).totup-ctlmeans(cidx).totdown)/ctlmeans(cidx).tot;
    end
end

BLMeansD1 = [mean(ctlStats(1).pctdown), mean(ctlStats(1).pctup), mean(ctlStats(1).pctcon);...
    mean(expStats(1).pctdown), mean(expStats(1).pctup), mean(expStats(1).pctcon)];
lowerBLErrD1 = [abs(mean(ctlStats(1).pctdown)-quantile(ctlStats(1).pctdown,.25)),...
    abs(mean(ctlStats(1).pctup)-quantile(ctlStats(1).pctup,.25)),...
    abs(mean(ctlStats(1).pctcon)-quantile(ctlStats(1).pctcon,.25));...
    abs(mean(expStats(1).pctdown)-quantile(expStats(1).pctdown,.25)),...
    abs(mean(expStats(1).pctup)-quantile(expStats(1).pctup,.25)),...
    abs(mean(expStats(1).pctcon)-quantile(expStats(1).pctcon,.25)),];
upperBLErrD1 = [abs(mean(ctlStats(1).pctdown)-quantile(ctlStats(1).pctdown,.75)),...
    abs(mean(ctlStats(1).pctup)-quantile(ctlStats(1).pctup,.75)),...
    abs(mean(ctlStats(1).pctcon)-quantile(ctlStats(1).pctcon,.75));...
    abs(mean(expStats(1).pctdown)-quantile(expStats(1).pctdown,.75)),...
    abs(mean(expStats(1).pctup)-quantile(expStats(1).pctup,.75)),...
    abs(mean(expStats(1).pctcon)-quantile(expStats(1).pctcon,.75)),];
BLMeansD2 = [mean(ctlStats(2).pctdown), mean(ctlStats(2).pctup), mean(ctlStats(2).pctcon);...
    mean(expStats(2).pctdown), mean(expStats(2).pctup), mean(expStats(2).pctcon)];
lowerBLErrD2 = [abs(mean(ctlStats(2).pctdown)-quantile(ctlStats(2).pctdown,.25)),...
    abs(mean(ctlStats(2).pctup)-quantile(ctlStats(2).pctup,.25)),...
    abs(mean(ctlStats(2).pctcon)-quantile(ctlStats(2).pctcon,.25));...
    abs(mean(expStats(2).pctdown)-quantile(expStats(2).pctdown,.25)),...
    abs(mean(expStats(2).pctup)-quantile(expStats(2).pctup,.25)),...
    abs(mean(expStats(2).pctcon)-quantile(expStats(2).pctcon,.25)),];
upperBLErrD2 = [abs(mean(ctlStats(2).pctdown)-quantile(ctlStats(2).pctdown,.75)),...
    abs(mean(ctlStats(2).pctup)-quantile(ctlStats(2).pctup,.75)),...
    abs(mean(ctlStats(2).pctcon)-quantile(ctlStats(2).pctcon,.75));...
    abs(mean(expStats(2).pctdown)-quantile(expStats(2).pctdown,.75)),...
    abs(mean(expStats(2).pctup)-quantile(expStats(2).pctup,.75)),...
    abs(mean(expStats(2).pctcon)-quantile(expStats(2).pctcon,.75)),];

[bt1, hb1, he1] = errorbar_groups(BLMeansD1, lowerBLErrD1, upperBLErrD1, ...
    'bar_colors', barColors, 'optional_errorbar_arguments',{'LineWidth', 2.5,'LineStyle','none'});
title('Day 1 Mean Baseline')
ylabel('Percentage of Cells per Mouse')
legend('Sham','Blast','Location','NorthWest')
[bt2, hb2, he2] = errorbar_groups(BLMeansD2, lowerBLErrD2, upperBLErrD2, ...
    'bar_colors', barColors, 'optional_errorbar_arguments',{'LineWidth', 2.5,'LineStyle','none'});
title('Day 2 Mean Baseline')
ylabel('Percentage of Cells per Mouse')
legend('Sham','Blast','Location','NorthWest')

%% Increase/Decrease/Nonmodulated Cells
afn{1} = 'activityResults_MMDDYYYY_Day1.mat';
afn{2} = 'activityResults_MMDDYYYY_Day2.mat';

for f = 1:numel(afn)
    load(afn{f}, 'exppcts', 'ctlpcts')
    for cidx=1:numel(exppcts)
        Ncells = numel(exppcts(cidx).inc);
        expStats(f).pctdec(cidx,1) = sum(exppcts(cidx).dec)/Ncells;
        expStats(f).pctinc(cidx,1) = sum(exppcts(cidx).inc)/Ncells;
        expStats(f).pctnm(cidx,1) = (Ncells-sum(exppcts(cidx).dec)-sum(exppcts(cidx).inc))/Ncells;
    end
    for cidx=1:numel(ctlpcts)
        Ncells = numel(ctlpcts(cidx).inc);
        ctlStats(f).pctdec(cidx,1) = sum(ctlpcts(cidx).dec)/Ncells;
        ctlStats(f).pctinc(cidx,1) = sum(ctlpcts(cidx).inc)/Ncells;
        ctlStats(f).pctnm(cidx,1) = (Ncells-sum(ctlpcts(cidx).dec)-sum(ctlpcts(cidx).inc))/Ncells;
    end
end

FRMeansD1 = [mean(ctlStats(1).pctdec), mean(ctlStats(1).pctinc), mean(ctlStats(1).pctnm);...
    mean(expStats(1).pctdec), mean(expStats(1).pctinc), mean(expStats(1).pctnm)];
lowerFRErrD1 = [abs(mean(ctlStats(1).pctdec)-quantile(ctlStats(1).pctdec,.25)),...
    abs(mean(ctlStats(1).pctinc)-quantile(ctlStats(1).pctinc,.25)),...
    abs(mean(ctlStats(1).pctnm)-quantile(ctlStats(1).pctnm,.25));...
    abs(mean(expStats(1).pctdec)-quantile(expStats(1).pctdec,.25)),...
    abs(mean(expStats(1).pctinc)-quantile(expStats(1).pctinc,.25)),...
    abs(mean(expStats(1).pctnm)-quantile(expStats(1).pctnm,.25)),];
upperFRErrD1 = [abs(mean(ctlStats(1).pctdec)-quantile(ctlStats(1).pctdec,.75)),...
    abs(mean(ctlStats(1).pctinc)-quantile(ctlStats(1).pctinc,.75)),...
    abs(mean(ctlStats(1).pctnm)-quantile(ctlStats(1).pctnm,.75));...
    abs(mean(expStats(1).pctdec)-quantile(expStats(1).pctdec,.75)),...
    abs(mean(expStats(1).pctinc)-quantile(expStats(1).pctinc,.75)),...
    abs(mean(expStats(1).pctnm)-quantile(expStats(1).pctnm,.75)),];
FRMeansD2 = [mean(ctlStats(2).pctdec), mean(ctlStats(2).pctinc), mean(ctlStats(2).pctnm);...
    mean(expStats(2).pctdec), mean(expStats(2).pctinc), mean(expStats(2).pctnm)];
lowerFRErrD2 = [abs(mean(ctlStats(2).pctdec)-quantile(ctlStats(2).pctdec,.25)),...
    abs(mean(ctlStats(2).pctinc)-quantile(ctlStats(2).pctinc,.25)),...
    abs(mean(ctlStats(2).pctnm)-quantile(ctlStats(2).pctnm,.25));...
    abs(mean(expStats(2).pctdec)-quantile(expStats(2).pctdec,.25)),...
    abs(mean(expStats(2).pctinc)-quantile(expStats(2).pctinc,.25)),...
    abs(mean(expStats(2).pctnm)-quantile(expStats(2).pctnm,.25)),];
upperFRErrD2 = [abs(mean(ctlStats(2).pctdec)-quantile(ctlStats(2).pctdec,.75)),...
    abs(mean(ctlStats(2).pctinc)-quantile(ctlStats(2).pctinc,.75)),...
    abs(mean(ctlStats(2).pctnm)-quantile(ctlStats(2).pctnm,.75));...
    abs(mean(expStats(2).pctdec)-quantile(expStats(2).pctdec,.75)),...
    abs(mean(expStats(2).pctinc)-quantile(expStats(2).pctinc,.75)),...
    abs(mean(expStats(2).pctnm)-quantile(expStats(2).pctnm,.75)),];

[bt1, hb1, he1] = errorbar_groups(FRMeansD1, lowerFRErrD1, upperFRErrD1, ...
    'bar_colors', barColors, 'optional_errorbar_arguments',{'LineWidth', 2.5,'LineStyle','none'});
title('Day 1 Activity Changes')
ylabel('Percentage of Cells per Mouse')
ylim([0 0.8])
legend('Sham','Blast','Location','NorthEast')
[bt2, hb2, he2] = errorbar_groups(FRMeansD2, lowerFRErrD2, upperFRErrD2, ...
    'bar_colors', barColors, 'optional_errorbar_arguments',{'LineWidth', 2.5,'LineStyle','none'});
title('Day 2 Activity Changes')
ylabel('Percentage of Cells per Mouse')
ylim([0 0.8])
legend('Sham','Blast','Location','NorthEast')
