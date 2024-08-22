clear;clc;close all;

addpath(genpath([pwd,'/model-fits/']))

% get fitted values
fits = load('Figure3_Fits.mat');

int_coef = [];
magprob_weight = [];
inv_temp = [];
for i=1:numel(fits.COMP.outputFull_B)
    int_coef = [int_coef fits.COMP.outputFull_B{i}.Xfit(1)];
    magprob_weight = [magprob_weight fits.COMP.outputFull_B{i}.Xfit(3)];
    inv_temp = [inv_temp fits.COMP.outputFull_B{i}.Xfit(2)];
end

% plot histograms
figure;
h1 = histogram(int_coef);
h1.FaceColor = [100 143 255]/256;
h1.LineWidth = 3;
set(gca,'fontsize',26);
set(gca,'linewidth',3);
xlim([0,1])
ylim([0,60])
title('Integration Coefficient')
box off
set(gca,'TickDir','out');

figure;
h2 = histogram(magprob_weight);
h2.FaceColor = [255 176 0]/256;
h2.LineWidth = 3;
set(gca,'fontsize',26);
set(gca,'linewidth',3);
xlim([0,1])
ylim([0,60])
title('Magnitude/Probability Weighting')
box off
set(gca,'TickDir','out');

figure;
h3 = histogram(inv_temp);
h3.FaceColor = [220 38 127]/255;
h3.LineWidth = 3;
set(gca,'fontsize',26);
set(gca,'linewidth',3);
xlim([0,50])
ylim([0,60])
title('Inverse Temperature')
box off
set(gca,'TickDir','out');
