%Script to make distribution plots of selected cells
%different traces

%Load Data
%Experimental Session
[expFile, basedir] = uigetfile('.mat', 'Select processed traces file for Experimental Session');
load(fullfile(basedir, expFile));
r_outexp = r_out;
clear r_out
%Control Session
[ctlFile, basedir] = uigetfile('.mat', 'Select processed traces file for Control Session');
load(fullfile(basedir, ctlFile));
r_outctl = r_out;
clear r_out

%Timing Traces
fs=20; %Sampling Frequency (Hz)
Npts = numel(r_outexp(1).file(1).trace);
t = 1:Npts;
time = struct();
increase = 0;
spacer = 0;
space = 200;
for idx=1:5
    time(idx).trace = (t+increase+spacer-1)/fs;
    increase = increase+Npts;
    spacer = spacer+space;
end
    
%%
%%%%% Experimental (Mouse 6 Day 1)
%Selected Pixels (By Hand from Maps)
selexpup = (1024*420) + 676; %(ImHeight*X) + Y
selexpdown = (1024*365) + 526;
selexpmid = (1024*326) + 588;
selexpupidx = 0;
selexpdownidx = 0;
selexpmididx = 0;
for idx = 1:numel(r_outexp)
    if sum(r_outexp(idx).pixel_idx == selexpup)
        selexpupidx = idx;
    elseif sum(r_outexp(idx).pixel_idx == selexpdown)
        selexpdownidx = idx;
    elseif sum(r_outexp(idx).pixel_idx == selexpmid)
        selexpmididx = idx;
    end
end

%Since ROIs changed for test
selexpdownidx = 1;

%Density Plots
%Put density values into structure for plotting
[dplots(1).f(:,1), dplots(1).xi(:,1)] = ksdensity(r_outexp(selexpupidx).file(1).trace);
[dplots(1).f(:,2), dplots(1).xi(:,2)] = ksdensity(r_outexp(selexpupidx).file(5).trace);
dplots(1).m(1) = mean(r_outexp(selexpupidx).file(1).trace);
dplots(1).m(2) = mean(r_outexp(selexpupidx).file(5).trace);
[dplots(2).f(:,1), dplots(2).xi(:,1)] = ksdensity(r_outexp(selexpdownidx).file(1).trace);
[dplots(2).f(:,2), dplots(2).xi(:,2)] = ksdensity(r_outexp(selexpdownidx).file(5).trace);
dplots(2).m(1) = mean(r_outexp(selexpdownidx).file(1).trace);
dplots(2).m(2) = mean(r_outexp(selexpdownidx).file(5).trace);
[dplots(3).f(:,1), dplots(3).xi(:,1)] = ksdensity(r_outexp(selexpmididx).file(1).trace);
[dplots(3).f(:,2), dplots(3).xi(:,2)] = ksdensity(r_outexp(selexpmididx).file(5).trace);
dplots(3).m(1) = mean(r_outexp(selexpmididx).file(1).trace);
dplots(3).m(2) = mean(r_outexp(selexpmididx).file(5).trace);
col = {[1,0,0], [0,0,1], [.5,.5,.5]};

figure();
for idx = 1:3
    subplot(3,1,idx)
    set(gca, 'fontsize',20) %30 is default
    hold on;
    patchline(dplots(idx).xi(:,1), dplots(idx).f(:,1), 'edgecolor', col{idx}, 'linewidth', 2, 'edgealpha', 0.6);
    patchline([dplots(idx).m(1), dplots(idx).m(1)], [2.1e-3, 2.7e-3], 'edgecolor', col{idx}, 'linewidth', 2, 'edgealpha', 0.6);
    patchline(dplots(idx).xi(:,2), dplots(idx).f(:,2), 'edgecolor', col{idx}, 'linewidth', 2);
    patchline([dplots(idx).m(2), dplots(idx).m(2)], [2.1e-3, 2.7e-3], 'edgecolor', col{idx}, 'linewidth', 2);
    xlim([0.6e4, 2.3e4])
    ylim([0, 3e-3])
    xlabel('Intensity (A.U.)')
    ylabel('PDE')
end

pos = [500, 300, 560, 810]; %[1000, 530, 560, 810];
set(gcf, 'Position', pos);

%Trace Plots
seltrace = [selexpupidx, selexpdownidx, selexpmididx];
figure();
for step = 1:3
    subplot(3,1,step)
    set(gca, 'fontsize',30)
    hold on;
    for idx = 1:5
        plot(time(idx).trace, r_outexp(seltrace(step)).file(idx).trace, 'color', col{step}, 'linewidth', 2);
    end
    xlim([0,560])
    ylim([0.6e4, 2.3e4])
    xlabel('Time (sec)')
    ylabel('Intensity (A.U.)')
end

pos = [61, 492, 2382, 842];
set(gcf, 'Position', pos);

%%
%%%%% Control (Mouse 11 Day 1)
%Selected Pixels (By Hand from Maps)
%None of these are actually up/down/mid as all are mid
selctlup = (1024*559) + 509; %(ImHeight*X) + Y
selctldown = (1024*433) + 188;
selctlmid = (1024*637) + 278;
selctlupidx = 0;
selctldownidx = 0;
selctlmididx = 0;
for idx = 1:numel(r_outctl)
    if sum(r_outctl(idx).pixel_idx == selctlup)
        selctlupidx = idx;
    elseif sum(r_outctl(idx).pixel_idx == selctldown)
        selctldownidx = idx;
    elseif sum(r_outctl(idx).pixel_idx == selctlmid)
        selctlmididx = idx;
    end
end

%Density Plots
%Put density values into structure for plotting
[dplots(1).f(:,1), dplots(1).xi(:,1)] = ksdensity(r_outctl(selctlupidx).file(1).trace);
[dplots(1).f(:,2), dplots(1).xi(:,2)] = ksdensity(r_outctl(selctlupidx).file(5).trace);
dplots(1).m(1) = mean(r_outctl(selctlupidx).file(1).trace);
dplots(1).m(2) = mean(r_outctl(selctlupidx).file(5).trace);
[dplots(2).f(:,1), dplots(2).xi(:,1)] = ksdensity(r_outctl(selctldownidx).file(1).trace);
[dplots(2).f(:,2), dplots(2).xi(:,2)] = ksdensity(r_outctl(selctldownidx).file(5).trace);
dplots(2).m(1) = mean(r_outctl(selctldownidx).file(1).trace);
dplots(2).m(2) = mean(r_outctl(selctldownidx).file(5).trace);
[dplots(3).f(:,1), dplots(3).xi(:,1)] = ksdensity(r_outctl(selctlmididx).file(1).trace);
[dplots(3).f(:,2), dplots(3).xi(:,2)] = ksdensity(r_outctl(selctlmididx).file(5).trace);
dplots(3).m(1) = mean(r_outctl(selctlmididx).file(1).trace);
dplots(3).m(2) = mean(r_outctl(selctlmididx).file(5).trace);
col = {[.5,.5,.5], [.5,.5,.5], [.5,.5,.5]};

figure();
for idx = 1:3
    subplot(3,1,idx)
    set(gca, 'fontsize',20)
    hold on;
    patchline(dplots(idx).xi(:,1), dplots(idx).f(:,1), 'edgecolor', col{idx}, 'linewidth', 2, 'edgealpha', 0.6);
    patchline([dplots(idx).m(1), dplots(idx).m(1)], [0.012, 0.014], 'edgecolor', col{idx}, 'linewidth', 2, 'edgealpha', 0.6);
    patchline(dplots(idx).xi(:,2), dplots(idx).f(:,2), 'edgecolor', col{idx}, 'linewidth', 2);
    patchline([dplots(idx).m(2), dplots(idx).m(2)], [0.012, 0.014], 'edgecolor', col{idx}, 'linewidth', 2);
    xlim([3000, 9000])
    ylim([0, .015])
    set(gca, 'YTick', [0, .005, .010, .015], 'YTickLabel',  [0, .005, .010, .015])
    xlabel('Intensity (A.U.)')
    ylabel('PDE')
end

pos = [500, 300, 560, 810];
set(gcf, 'Position', pos);

%Trace Plots
seltrace = [selctlupidx, selctldownidx, selctlmididx];
figure();
for step = 1:3
    subplot(3,1,step)
    set(gca, 'fontsize',30)
    hold on;
    for idx = 1:5
        plot(time(idx).trace, r_outctl(seltrace(step)).file(idx).trace, 'color', col{step}, 'linewidth', 2);
    end
    xlim([0,560])
    ylim([3000, 9000])
    xlabel('Time (sec)')
    ylabel('Intensity (A.U.)')
end

pos = [61, 492, 2382, 842];
set(gcf, 'Position', pos);

