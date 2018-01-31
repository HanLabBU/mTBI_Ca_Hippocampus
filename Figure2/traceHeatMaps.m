%Script to plot heat maps of traces

%Load File
%Select processedTraces_*.mat for before blast to load and plot
[befTracesFile, basedir] = uigetfile('.mat', 'Select processed traces file for Before Blast');
load(fullfile(basedir, befTracesFile));
befTraces = r_out;
clear r_out

%Select processedTraces_*.mat for after blast to load and plot
[aftTracesFile, basedir] = uigetfile('.mat', 'Select processed traces file for After Blast');
load(fullfile(basedir, aftTracesFile));
aftTraces = r_out;
clear r_out

%Reshape Structures for Plots
%Before Blast Reshaping
allBefTraces = zeros(numel(befTraces), numel(befTraces(1).file(1).fnormtrace));
for idx = 1:numel(befTraces)
    allBefTraces(idx,:) = 100 * befTraces(idx).file(1).fnormtrace - 100; %Multiply by 100 to get percent
end
%After Blast Reshaping
allAftTraces = zeros(numel(aftTraces(1).file), numel(aftTraces), numel(aftTraces(1).file(1).fnormtrace));
for fiter = 1:numel(aftTraces(1).file) %File Iteration
    for idx = 1:numel(aftTraces) %Cell Iteration
        allAftTraces(fiter,idx,:) = 100 * aftTraces(idx).file(fiter).fnormtrace - 100; %Multiply by 100 to get percent
    end
end

%Do Sorting for plotting
befMeans = mean(allBefTraces,2);
[~, befInds] = sort(befMeans, 'descend');
aftMeans = mean(squeeze(allAftTraces(1,:,:)),2);
[~, aftInds] = sort(aftMeans, 'descend');

%Before
figure(1);
subplot(1, 6, 1);
imagesc(allBefTraces(befInds,:))
caxis([-20, 20])
colormap(redblue)
ylabel('Cell Number')
xlabel('Time (Frame #)')

%After
for idx = 2:numel(aftTraces(1).file)+1
    subplot(1, 6, idx);
    imagesc(squeeze(allAftTraces(idx-1, aftInds, :)))
    caxis([-20, 20])
    colormap(redblue)
    xlabel('Time (Frame #)')
    if idx == numel(aftTraces(1).file)+1
        colorbar
    end
end