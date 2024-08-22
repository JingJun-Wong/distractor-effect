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

% correlations and scatterplots
x1 = int_coef;
y1 = magprob_weight;
x1text = 'Integration Coefficient';
y1text = 'Magnitude/Probability Weighting';

x2 = int_coef;
y2 = inv_temp;
x2text = 'Integration Coefficient';
y2text = 'Inverse Temperature';

x3 = magprob_weight;
y3 = inv_temp;
x3text = 'Magnitude/Probability Weighting';
y3text = 'Inverse Temperature';

[rvalue_ic_mpw,pvalue_ic_mpw]=corr(x1',y1','Type','Pearson')
[rvalue_ic_it,pvalue_ic_it]=corr(x2',y2','Type','Pearson')
[rvalue_mpw_it,pvalue_mpw_it]=corr(x3',y3','Type','Pearson')

coefficients1 = polyfit(x1, y1, 1);
x1Fit = linspace(0.05, 0.95, 1000);
y1Fit = polyval(coefficients1 , x1Fit);

figure;
sz = 200;
c = [0.5 0.5 0.5];
scatter(x1,y1,sz,c,'filled')
hold on;
plot(x1Fit, y1Fit, 'k-', 'LineWidth', 3);
xlabel(x1text)
ylabel(y1text)
xlim([0,1])
xticks([0:0.2:1])
ylim([0,1])
yticks([0:0.2:1])
set(gca,'TickDir','out');
box off;
set(gca,'fontsize',26);
set(gca,'linewidth',3);

coefficients2 = polyfit(x2, y2, 1);
x2Fit = linspace(0.05, 0.95, 1000);
y2Fit = polyval(coefficients2 , x2Fit);

figure;
sz = 200;
c = [0.5 0.5 0.5];
scatter(x2,y2,sz,c,'filled')
hold on;
plot(x2Fit, y2Fit, 'k-', 'LineWidth', 3);
xlabel(x2text)
ylabel(y2text)
xlim([0,1])
xticks([0:0.2:1])
ylim([0,50])
yticks([0:10:50])
set(gca,'TickDir','out');
box off;
set(gca,'fontsize',26);
set(gca,'linewidth',3);

coefficients3 = polyfit(x3, y3, 1);
x3Fit = linspace(0.05, 0.95, 1000);
y3Fit = polyval(coefficients3 , x3Fit);

figure;
sz = 200;
c = [0.5 0.5 0.5];
scatter(x3,y3,sz,c,'filled')
hold on;
plot(x3Fit, y3Fit, 'k-', 'LineWidth', 3);
xlabel(x3text)
ylabel(y3text)
xlim([0,1])
xticks([0:0.2:1])
ylim([0,50])
yticks([0:10:50])
set(gca,'TickDir','out');
box off;
set(gca,'fontsize',26);
set(gca,'linewidth',3);
