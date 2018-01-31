%Script to make scatter plots of mean baseline values compared between
%different traces

%Find files
basedir = '/ParentDirectory/ForAll/processedTraces/AfterBlasting';
[~, list] = system(sprintf('find "%s" -type f -name "processedTraces_*.mat"', basedir)); %Currently written for Unix, but can be done on Windows/Linux

breaks = find(list == char(10));
cell_list = cell(1,numel(breaks));
start = 1;
for idx = 1:numel(breaks)
    if idx == 1
        cell_list{idx} = list(1:breaks(idx)-1);
    else
        cell_list{idx} = list(breaks(idx-1)+1:breaks(idx)-1);
    end
end

%Identify Types from cell lists
%These numbers correspond to the corresponding position of the data found
%in cell_list above from the recursive search
%Day 1
CTLorder_aft = [9, 1, 3, 13, 15];  %8, 10, 11, M1, M2
TBIorder_aft = [7, 11, 5, 17]; %6, 9, 12, S5
% %Day2
% CTLorder_aft = [10, 2, 4, 14, 16];  %8, 10, 11, M1, M2
% TBIorder_aft = [8, 12, 6, 18]; %6, 9, 12, S5
% %All Cells
% CTLorder_aft = [9, 10, 1, 2, 3, 4, 13, 14, 15, 16];  %8, 10, 11, M1, M2
% TBIorder_aft = [7, 8, 11, 12, 5, 6, 17, 18]; %6, 9, 12, S5

%Load Data
expmeans = struct();
ctlmeans = struct();
expcnt = 1;
ctlcnt = 1;
for idx = 1:numel(cell_list)
    fn = cell_list{idx};
    load(fn);
    if ismember(idx, TBIorder_aft) %Experimental Mice
        [expmeans(expcnt).one, expmeans(expcnt).five] = generateMeanValues(r_out, 1, 5, 'trace');
        expmeans(expcnt).filename = fn;
        expcnt = expcnt + 1;
    elseif ismember(idx, CTLorder_aft) %Control Mice
        [ctlmeans(ctlcnt).one, ctlmeans(ctlcnt).five] = generateMeanValues(r_out, 1, 5, 'trace');
        ctlmeans(ctlcnt).filename = fn;
        ctlcnt = ctlcnt + 1;
    end
end

%Plot Population of Cells
Nexp = 0;
Nctl = 0;
figure();
set(gca, 'fontsize',30)
hold on;
for idx = 1:numel(expmeans)
    hs(idx) = scatter(expmeans(idx).one, expmeans(idx).five, 'or');
    Nexp = Nexp + numel(expmeans(idx).one);
    %pause()
end
for idx = 1:numel(ctlmeans)
    hs(numel(expmeans)+idx) = scatter(ctlmeans(idx).one, ctlmeans(idx).five, 'sb');
    Nctl = Nctl + numel(ctlmeans(idx).one);
    %pause()
end
plot([0,2^16],[0,2^16],'--k', 'linewidth', 1.5);
xlabel('Time Period 1 (A.U.)');
ylabel('Time Period 5 (A.U.)');
xlim([0,2^16])
ylim([0,2^16])
legend(hs([1,numel(expmeans)+1]),'Blast','Sham','Location','Northwest')
%title('Mean of Trace for Time Periods across All Cells')

pos = [1016, 527, 814, 658];
set(gcf, 'Position', pos);
pause(1);
set(gca, 'XTick', get(gca, 'YTick'))

%Combine Population of All Cells
allexp = zeros(Nexp,2); %1st Column is 1, 2nd Column is 5
allctl = zeros(Nctl,2);
expstart = 1;
ctlstart = 1;
s1 = 0;
s2 = 0;
for idx = 1:numel(expmeans)
    nexp = numel(expmeans(idx).one);
    allexp(expstart:(expstart+nexp-1), 1) = expmeans(idx).one;
    allexp(expstart:(expstart+nexp-1), 2) = expmeans(idx).five;
    expstart = expstart+nexp;
end
for idx = 1:numel(ctlmeans)
    nctl = numel(ctlmeans(idx).one);
    allctl(ctlstart:(ctlstart+nctl-1), 1) = ctlmeans(idx).one;
    allctl(ctlstart:(ctlstart+nctl-1), 2) = ctlmeans(idx).five;
    ctlstart = ctlstart+nctl;
end

%Confidence Interval
ctldiff = allctl(:,1)-allctl(:,2); %Difference value
CI99 = quantile(ctldiff-mean(ctldiff), [.005, .995]);
CI95 = quantile(ctldiff-mean(ctldiff), [.025, .975]);

%Add 95% CIs to Above Plot
plot([0,2^16],[0+abs(CI95(1)),2^16+abs(CI95(1))],':k','linewidth',1.5);
plot([0+abs(CI95(2)),2^16+abs(CI95(2))],[0,2^16],':k','linewidth',1.5);

%Significant Cells
expdiff = allexp(:,1)-allexp(:,2); %Difference of Significant Cells
expdown = expdiff < CI95(1);
expup = expdiff > CI95(2);
ctldown = ctldiff < CI95(1);
ctlup = ctldiff > CI95(2);

%Significant Cells for Example Mice
%Experimental Example
expdiffM6D1.diffvals = expmeans(3).one - expmeans(3).five; %Idx 3 is Mouse 6 Day 1 (6 in cell_list)
expdiffM6D1.sigdown = expdiffM6D1.diffvals < CI95(1);
expdiffM6D1.sigup =  expdiffM6D1.diffvals > CI95(2);
%Control Example
expdiffM11D1.diffvals = ctlmeans(3).one - ctlmeans(3).five; %Idx 3 is Mouse 11 Day 1 (12 in cell_list)
expdiffM11D1.sigdown = expdiffM11D1.diffvals < CI95(1);
expdiffM11D1.sigup =  expdiffM11D1.diffvals > CI95(2);

%save('sigexamples_MMDDYYYY.mat','expdiffM6D1','expdiffM11D1')

%Each experiment percentage up and down
for idx = 1:numel(expmeans)
    expmeans(idx).diff = expmeans(idx).one-expmeans(idx).five;
    expmeans(idx).up = expmeans(idx).diff >CI95(2);
    expmeans(idx).down  = expmeans(idx).diff < CI95(1);
    expmeans(idx).totup = sum(expmeans(idx).up);
    expmeans(idx).totdown = sum(expmeans(idx).down);
    expmeans(idx).tot = numel(expmeans(idx).up);
    expmeans(idx).pctup = expmeans(idx).totup / expmeans(idx).tot;
    expmeans(idx).pctdown = expmeans(idx).totdown / expmeans(idx).tot;
end
for idx = 1:numel(ctlmeans)
    ctlmeans(idx).diff = ctlmeans(idx).one-ctlmeans(idx).five;
    ctlmeans(idx).up = ctlmeans(idx).diff > CI95(2);
    ctlmeans(idx).down = ctlmeans(idx).diff < CI95(1);
    ctlmeans(idx).totup = sum(ctlmeans(idx).up);
    ctlmeans(idx).totdown = sum(ctlmeans(idx).down);
    ctlmeans(idx).tot = numel(ctlmeans(idx).up);
    ctlmeans(idx).pctup = ctlmeans(idx).totup / ctlmeans(idx).tot;
    ctlmeans(idx).pctdown = ctlmeans(idx).totdown / ctlmeans(idx).tot;
end

%Average and Standard Deviation
exppctup = nan(6,1);
exppctdown = nan(6,1);
exppctnm = nan(6,1);
ctlpctup = nan(6,1);
ctlpctdown = nan(6,1);
ctlpctnm = nan(6,1);
for idx = 1:numel(expmeans)
    exppctup(idx) = expmeans(idx).pctup;
    exppctdown(idx) = expmeans(idx).pctdown;
    exppctnm(idx) = (expmeans(idx).tot - expmeans(idx).totup - expmeans(idx).totdown) /expmeans(idx).tot;
end
for idx = 1:numel(ctlmeans)
    ctlpctup(idx) = ctlmeans(idx).pctup;
    ctlpctdown(idx) = ctlmeans(idx).pctdown;
    ctlpctnm(idx) = (ctlmeans(idx).tot - ctlmeans(idx).totup - ctlmeans(idx).totdown) /ctlmeans(idx).tot;
end

[pup, hup, statsup] = ranksum(ctlpctup, exppctup);
[pdown, hdown, statsdown] = ranksum(ctlpctdown, exppctdown);
[pnm, hnm, statsnm] = ranksum(ctlpctnm, exppctnm);

%Comparisons for Before Data
%Before Blasting Data & Comparisons
%Find files
befbasedir = '/ParentDirectory/ForAll/ProcessedTraces/BeforeBlasting';
[~, beflist] = system(sprintf('find "%s" -type f -name "processedTraces_*.mat"', befbasedir));  %Currently written for Unix, but can be done on Windows/Linux

befbreaks = find(beflist == char(10));
befcell_list = cell(1,numel(befbreaks));
start = 1;
for idx = 1:numel(befbreaks)
    if idx == 1
        befcell_list{idx} = beflist(1:(befbreaks(idx)-1));
    else
        befcell_list{idx} = beflist((befbreaks(idx-1)+1):(befbreaks(idx)-1));
    end
end

%Identify Types from cell lists for Before
%These numbers correspond to the corresponding position of the data found
%in cell_list above from the recursive search
%Day 1
CTLorder_bef = [9, 1, 3, 13, 15];  %8, 10, 11, M1, M2
TBIorder_bef = [7, 11, 5, 17]; %6, 9, 12, S5
% %Day 2
% CTLorder_bef = [10, 2, 4, 14, 16];  %8, 10, 11, M1, M2
% TBIorder_bef = [8, 12, 6, 18]; %6, 9, 12, S5
% %All Cells
% CTLorder_bef = [9, 10, 1, 2, 3, 4, 13, 14, 15, 16];  %8, 10, 11, M1, M2
% TBIorder_bef = [7, 8, 11, 12, 5, 6, 17, 18]; %6, 9, 12, S5

%Load Data
Fs = 20; %20 Hz sampling
samptime = 10; %Seconds to sample for resampling (10 seconds)
Nsamps = samptime * Fs;
befexpmeans = struct();
befctlmeans = struct();
expcnt = 1;
ctlcnt = 1;
for idx = 1:numel(befcell_list)
    fn = befcell_list{idx};
    load(fn);
    %%%%%%% Experimental Mice
    if ismember(idx, TBIorder_bef)
        [befexpmeans(expcnt).one, exptest] = generateMeanValues(r_out, 1, [], 'trace');
        expcnt = expcnt + 1;
    %%%%%%% Control Mice
    elseif ismember(idx, CTLorder_bef)
        [befctlmeans(ctlcnt).one, ctltest] = generateMeanValues(r_out, 1, [], 'trace');
        ctlcnt = ctlcnt + 1;
    end
end

nmice = numel(befexpmeans);
expplots = struct();
ctlplots = struct();
for idx = 1:numel(befexpmeans)
    %Experimental
    [expplots(idx).f(:,1), expplots(idx).xi(:,1)] = ksdensity(befexpmeans(idx).one); %Pre Blast
    [expplots(idx).f(:,2), expplots(idx).xi(:,2)] = ksdensity(expmeans(idx).one); %Post Blast Sessions
    [expplots(idx).f(:,3), expplots(idx).xi(:,3)] = ksdensity(expmeans(idx).five);
end
for idx = 1:numel(befctlmeans)
    %Control
    [ctlplots(idx).f(:,1), ctlplots(idx).xi(:,1)] = ksdensity(befctlmeans(idx).one); %Pre Blast
    [ctlplots(idx).f(:,2), ctlplots(idx).xi(:,2)] = ksdensity(ctlmeans(idx).one); %Post Blast Sessions
    [ctlplots(idx).f(:,3), ctlplots(idx).xi(:,3)] = ksdensity(ctlmeans(idx).five);
end

for choice = 1:numel(befexpmeans)
    col = {[1,0,0], [0,0,1], [.5,.5,.5]};
    alphas = [.3, .6, 1];
    figure(); hold on;
    for idx = 1:3
        patchline(expplots(choice).xi(:,idx), expplots(choice).f(:,idx), 'edgecolor', col{idx}, 'linewidth', 2);
    end
    title('Exp')
    xlabel('Intensity (A.U.)')
    ylabel('PDE')
end

for choice = 1:numel(befctlmeans)
    col = {[1,0,0], [0,0,1], [.5,.5,.5]};
    alphas = [.3, .6, 1];
    figure(); hold on;
    for idx = 1:3
        patchline(ctlplots(choice).xi(:,idx), ctlplots(choice).f(:,idx), 'edgecolor', col{idx}, 'linewidth', 2);
    end
    title('Ctl')
    xlabel('Intensity (A.U.)')
    ylabel('PDE')
end

%save('sigpopulation_MMDDYYYY_DayX.mat','expmeans','ctlmeans','expdown','expup','ctldown','ctlup','befexpmeans','befctlmeans','CTLorder_bef','CTLorder_aft','TBIorder_bef','TBIorder_aft')