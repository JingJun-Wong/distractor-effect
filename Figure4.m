clear;clc;close all;

addpath(genpath([pwd,'/model-fits/']))

% get fitted values and betas
fits = load('Figure3_Fits.mat');
betas = load('Figure2_Betas.mat');

int_coef = [];
beta_dvhv_t = [];
for i=1:numel(fits.COMP.outputFull_B)
    int_coef = [int_coef fits.COMP.outputFull_B{i}.Xfit(1)];
    beta_dvhv_t = [beta_dvhv_t betas.betas{3}(i,7)];
end

% correlation and scatterplot
x = beta_dvhv_t;
y = int_coef;
xtext= 'Beta (DV-HV)T';
ytext= 'Inegration Coefficient';

[rvalue_beta_ic,pvalue_beta_ic]=corr(x',y','Type','Pearson')

coefficients = polyfit(x, y, 1);
xFit = linspace(-0.9, 0.9, 1000);
yFit = polyval(coefficients , xFit);

figure;
sz = 200;
c = [0.5 0.5 0.5];
scatter(x,y,sz,c,'filled')
hold on;
plot(xFit, yFit, 'k-', 'LineWidth', 3);
xlabel(xtext)
ylabel(ytext)
xticks([-1:0.5:1])
ylim([0,1])
yticks([0:0.2:1])
set(gca,'TickDir','out');
title('Correlation between Integration Coefficient and Beta (DV-HV)T')
box off;
set(gca,'fontsize',26);
set(gca,'linewidth',3);
