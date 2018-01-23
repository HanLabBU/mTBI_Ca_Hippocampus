%Code to plot color coded up/down spatial maps of cells.
%Uses sigexamples.mat to gain up/down information
%Needs r_out for both the example cells for selection.

%Load Data
%r_out sections
expfile = '/Volumes/eng_research_handata/Hua-an/Data/TBI/Ali/TBI_AFTER BLASTING DATA ANALYSIS/Ali-TBI-6/Ali-TBI-6-Day1/R_bins_circle_20170804.mat';
load(expfile);
r_outexp = r_out;
ctlfile = '/Volumes/eng_research_handata/Hua-an/Data/TBI/Ali/TBI_AFTER BLASTING DATA ANALYSIS/Ali-TBI-11/Ali-Control-TBI-11-Day1/R_bins_circle_20170804.mat';
load(ctlfile);
r_outctl = r_out;
clear r_out;
%Up/Down Data
load('sigexamples_20170804.mat')

%Generate Values
sz = [1024,1024]; %Image Size
expmapR = zeros(sz); expmapG = zeros(sz); expmapB = zeros(sz);
ctlmapR = zeros(sz); ctlmapG = zeros(sz); ctlmapB = zeros(sz);
%Gray Values with Up Values (Red) and Down Values (Blue)
for idx = 1:numel(r_outexp)
    pix = r_outexp(idx).pixel_idx;
    if expdiffM6D1.sigup(idx)
        expmapR(pix) = 1;
    elseif expdiffM6D1.sigdown(idx)
        expmapB(pix) = 1;
    elseif expdiffM6D1.sigup(idx) && expdiffM6D1.sigdown(idx)
        print('Uh Oh') %Shouldn't have cell classified as both
    else
        expmapR(pix) = 0.5;
        expmapG(pix) = 0.5;
        expmapB(pix) = 0.5;
    end
end
expmap = cat(3,expmapR,expmapG,expmapB);
for idx = 1:numel(r_outctl)
    pix = r_outctl(idx).pixel_idx;
    if expdiffM11D1.sigup(idx)
        ctlmapR(pix) = 1;
    elseif expdiffM11D1.sigdown(idx)
        ctlmapB(pix) = 1;
    elseif expdiffM11D1.sigup(idx) && expdiffM11D1.sigdown(idx)
        print('Uh Oh') %Shouldn't have cell classified as both
    else
        ctlmapR(pix) = 0.5;
        ctlmapG(pix) = 0.5;
        ctlmapB(pix) = 0.5;
    end
end
ctlmap = cat(3,ctlmapR,ctlmapG,ctlmapB);

%Make Plots
figure(); imshow(expmap);
figure(); imshow(ctlmap);