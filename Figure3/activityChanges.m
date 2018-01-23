%Script to make scatter plots of mean baseline values compared between
%different traces

%Find files
basedir = '/Volumes/eng_research_handata/Hua-an/Data/TBI/Ali/TBI_AFTER BLASTING DATA ANALYSIS';
[~, list] = system(sprintf('find "%s" -type f -name "R_bins_circle_20170804.mat"', basedir));

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
%Day 1
CTLorder_aft = [9, 1, 3, 13, 15];  %8, 10, 11, Moona1, Moona2
TBIorder_aft = [7, 11, 5, 17]; %6, 9, 12, Seth5
% %Day2
% CTLorder_aft = [10, 2, 4, 14, 16];  %8, 10, 11, Moona1, Moona2
% TBIorder_aft = [8, 12, 6, 18]; %6, 9, 12, Seth5
% %All Cells
% CTLorder_aft = [9, 10, 1, 2, 3, 4, 13, 14, 15, 16];  %8, 10, 11, Moona1, Moona2
% TBIorder_aft = [7, 8, 11, 12, 5, 6, 17, 18]; %6, 9, 12, Seth5

%Load Data
Fs = 20; %20 Hz sampling
samptime = 10; %Seconds to sample for resampling (10 seconds)
Nsamps = samptime * Fs;
exppcts = struct();
ctlpcts = struct();
expcnt = 1;
ctlcnt = 1;
for idx = 1:numel(cell_list)
    fn = cell_list{idx};
    load(fn);
    if ismember(idx, TBIorder_aft) %Experimental Mice
        [exppcts(expcnt).one, exppcts(expcnt).five] = randomBinaryResample(r_out, 1, 5, 'bintrace3', Nsamps);
        expcnt = expcnt + 1;
    elseif ismember(idx, CTLorder_aft) %Control Mice
        [ctlpcts(ctlcnt).one, ctlpcts(ctlcnt).five] = randomBinaryResample(r_out, 1, 5, 'bintrace3', Nsamps);
        ctlcnt = ctlcnt + 1;
    end
end

%Perform Statistical Comparisons
expcnt = 1;
ctlcnt = 1;
for idx = 1:numel(cell_list)
    if ismember(idx, TBIorder_aft)
        ncells = size(exppcts(expcnt).one,2);
        p = nan(ncells,1);
        h = nan(ncells,1);
        stats = struct('zval',[],'ranksum',[]);
        diff = nan(ncells,1);
        inc = zeros(ncells,1);
        dec = zeros(ncells,1);
        odd = zeros(ncells,1);
        for cel = 1:ncells
            [p(cel), h(cel), stats(cel)] = ranksum(exppcts(expcnt).one(:,cel), exppcts(expcnt).five(:,cel));
            diff(cel) = mean(exppcts(expcnt).one(:,cel)-exppcts(expcnt).five(:,cel));
            if h(cel) == 1
                if diff(cel) > 0
                    inc(cel) = 1;
                elseif diff(cel) < 0
                    dec(cel) = 1;
                else
                    odd(cel) = 1;
                end
            end
        end
        exppcts(expcnt).stats.p = p;
        exppcts(expcnt).stats.h = h;
        exppcts(expcnt).stats.stats = stats;
        exppcts(expcnt).diff = diff;
        exppcts(expcnt).inc = inc;
        exppcts(expcnt).dec = dec;
        exppcts(expcnt).odd = odd;
        expcnt = expcnt + 1;
    elseif ismember(idx, CTLorder_aft)
        ncells = size(ctlpcts(ctlcnt).one,2);
        p = nan(ncells,1);
        h = nan(ncells,1);
        stats = struct('zval',[],'ranksum',[]);
        diff = nan(ncells,1);
        inc = zeros(ncells,1);
        dec = zeros(ncells,1);
        odd = zeros(ncells,1);
        for cel = 1:ncells
            [p(cel), h(cel), stats(cel)] = ranksum(ctlpcts(ctlcnt).one(:,cel), ctlpcts(ctlcnt).five(:,cel));
            diff(cel) = mean(ctlpcts(ctlcnt).one(:,cel)-ctlpcts(ctlcnt).five(:,cel));
            if h(cel) == 1
                if diff(cel) > 0
                    inc(cel) = 1;
                elseif diff(cel) < 0
                    dec(cel) = 1;
                else
                    odd(cel) = 1;
                end
            end
        end
        ctlpcts(ctlcnt).stats.p = p;
        ctlpcts(ctlcnt).stats.h = h;
        ctlpcts(ctlcnt).stats.stats = stats;
        ctlpcts(ctlcnt).diff = diff;
        ctlpcts(ctlcnt).inc = inc;
        ctlpcts(ctlcnt).dec = dec;
        ctlpcts(ctlcnt).odd = odd;
        ctlcnt = ctlcnt + 1;
    end
end

%Calculate Percentages %%Run this part after loading output to regenerate
%Table 2 in Paper
nmice_exp = numel(exppcts);
nmice_ctl = numel(ctlpcts);
decexp = nan(nmice_exp,1);
incexp = nan(nmice_exp,1);
decctl = nan(nmice_ctl,1);
incctl = nan(nmice_ctl,1);
totexp = nan(nmice_exp,1);
totctl = nan(nmice_ctl,1);
for idx = 1:nmice_exp
    decexp(idx) = sum(exppcts(idx).dec);
    incexp(idx) = sum(exppcts(idx).inc);
    totexp(idx) = size(exppcts(idx).one,2);
end
for idx = 1:nmice_ctl
    decctl(idx) = sum(ctlpcts(idx).dec);
    incctl(idx) = sum(ctlpcts(idx).inc);
    totctl(idx) = size(ctlpcts(idx).one,2);
end
sigexp = decexp+incexp;
sigctl = decctl+incctl;
nmexp = totexp-sigexp;
nmctl = totctl-sigctl;
pctdecexp = decexp./totexp;
pctincexp = incexp./totexp;
pctnmexp = nmexp./totexp;
pctdecctl = decctl./totctl;
pctincctl = incctl./totctl;
pctnmctl = nmctl./totctl;

[psigdec, hsigdec, statssigdec] = ranksum(pctdecctl, pctdecexp);
[psiginc, hsiginc, statssiginc] = ranksum(pctincctl, pctincexp);
[psignm, hsignm, statssignm] = ranksum(pctnmctl, pctnmexp);

%Before Blasting Data & Comparisons
%Find files
befbasedir = '/Volumes/eng_research_handata/Hua-an/Data/TBI/Ali/TBI_BEFORE BLASTING DATA ANALYSIS';
[~, beflist] = system(sprintf('find "%s" -type f -name "R_bins_circle_20170804.mat"', befbasedir));

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
%Day 1
CTLorder_bef = [9, 1, 3, 13, 15];  %8, 10, 11, Moona1, Moona2
TBIorder_bef = [7, 11, 5, 17]; %6, 9, 12, Seth5
% %Day 2
% CTLorder_bef = [10, 2, 4, 14, 16];  %8, 10, 11, Moona1, Moona2
% TBIorder_bef = [8, 12, 6, 18]; %6, 9, 12, Seth5
% %All Cells
% CTLorder_bef = [9, 10, 1, 2, 3, 4, 13, 14, 15, 16];  %8, 10, 11, Moona1, Moona2
% TBIorder_bef = [7, 8, 11, 12, 5, 6, 17, 18]; %6, 9, 12, Seth5

%Load Data
Fs = 20; %20 Hz sampling
samptime = 10; %Seconds to sample for resampling (10 seconds)
Nsamps = samptime * Fs;
befexppcts = struct();
befctlpcts = struct();
expcnt = 1;
ctlcnt = 1;
for idx = 1:numel(befcell_list)
    fn = befcell_list{idx};
    load(fn);
    %%%%%%% Experimental Mice
    if ismember(idx, TBIorder_bef)
        [befexppcts(expcnt).one, exptest] = randomBinaryResample(r_out, 1, [], 'bintrace3', Nsamps);
        expcnt = expcnt + 1;
    %%%%%%% Control Mice
    elseif ismember(idx, CTLorder_bef)
        [befctlpcts(ctlcnt).one, ctltest] = randomBinaryResample(r_out, 1, [], 'bintrace3', Nsamps);
        ctlcnt = ctlcnt + 1;
    end
end

nmice_exp = numel(befexppcts);
nmice_ctl = numel(befctlpcts);
expplots = struct();
ctlplots = struct();
for idx = 1:nmice_exp
    %Experimental
    befexppcts(idx).means = mean(befexppcts(idx).one); %Pre Blast
    exppcts(idx).onemeans = mean(exppcts(idx).one); %Blasting Periods
    exppcts(idx).fivemeans = mean(exppcts(idx).five);
    [expplots(idx).f(:,1), expplots(idx).xi(:,1)] = ksdensity(befexppcts(idx).means);
    [expplots(idx).f(:,2), expplots(idx).xi(:,2)] = ksdensity(exppcts(idx).onemeans);
    [expplots(idx).f(:,3), expplots(idx).xi(:,3)] = ksdensity(exppcts(idx).fivemeans);
end
for idx = 1:nmice_ctl
    %Control
    befctlpcts(idx).means = mean(befctlpcts(idx).one);
    ctlpcts(idx).onemeans = mean(ctlpcts(idx).one);
    ctlpcts(idx).fivemeans = mean(ctlpcts(idx).five);
    [ctlplots(idx).f(:,1), ctlplots(idx).xi(:,1)] = ksdensity(befctlpcts(idx).means);
    [ctlplots(idx).f(:,2), ctlplots(idx).xi(:,2)] = ksdensity(ctlpcts(idx).onemeans);
    [ctlplots(idx).f(:,3), ctlplots(idx).xi(:,3)] = ksdensity(ctlpcts(idx).fivemeans);
end

for choice = 1:nmice_exp
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

for choice = 1:nmice_ctl
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

%save('activityResults_20170918_Day2.mat','exppcts','ctlpcts','cell_list','befexppcts','befctlpcts','befcell_list','expplots','ctlplots','CTLorder_bef','CTLorder_aft','TBIorder_bef','TBIorder_aft')