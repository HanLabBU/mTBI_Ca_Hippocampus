%Script to make scatter plots of mean baseline values compared between
%different traces

%Load Data
expfile = '/Volumes/eng_research_handata/Hua-an/Data/TBI/Ali/TBI_AFTER BLASTING DATA ANALYSIS/Ali-TBI-6/Ali-TBI-6-Day1/R_bins_circle_20170804.mat';
load(expfile);
r_outexp = r_out;
ctlfile = '/Volumes/eng_research_handata/Hua-an/Data/TBI/Ali/TBI_AFTER BLASTING DATA ANALYSIS/Ali-TBI-11/Ali-Control-TBI-11-Day1/R_bins_circle_20170804.mat';
load(ctlfile);
r_outctl = r_out;

%Run Code for Plot
[expmean1, expmean5] = generateMeanValues(r_outexp,1,5,'trace');
[ctlmean1, ctlmean5] = generateMeanValues(r_outctl,1,5,'trace');

figure();
set(gca, 'fontsize',30)
scatter(expmean1, expmean5, 'or');
hold on;
scatter(ctlmean1, ctlmean5, 'sb');
plot([0,2.5e4],[0,2.5e4],'--k','linewidth',1.5);
xlabel('Time Period 1 (A.U.)');
ylabel('Time Period 5 (A.U.)');
xlim([0,2.5e4])
legend('Blast','Sham','Location','Northwest')
%title('Mean of Trace for Time Periods')

pos = [1016, 527, 814, 658];
set(gcf, 'Position', pos);
pause(1);
yticks = get(gca, 'YTick');
set(gca, 'XTick', yticks);
