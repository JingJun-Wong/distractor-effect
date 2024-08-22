% This script modified based on the scripts provided by Cao & Tsetsos, 2022
% https://github.com/YinanCao/multiattribute-distractor

clear;clc;close all;

addpath(genpath([pwd,'/model-fits/']))

fits = load('Figure3_Fits.mat');
ns=144;
model = {fits.AU,fits.EV,fits.EVDN,fits.COMP,fits.PT};
xx = [];
for j = 1:numel(model)
    CV_LL_B=model{j}.testLL_B(1:ns,:,:);
    xx = cat(2,xx,CV_LL_B);
end
xx = mean(xx,3);
xx = xx(1:ns,:);

% model comparison
[posterior1,out1] = VBA_groupBMC(xx');

% plot model attributions
hold on;
figure;
colormap(flipud(bone));
hi = imagesc(posterior1.r');
xlabel('Models')
ylabel('Subjects')
xlim([0.5,5.5])
ylim([0,144])
xticks([1:5])
xticklabels({'AU','EV','EV+DN','COMP','PT'})
yticks([1 50 100 144]);
set(gca,'TickDir','out');
box off;
colorbar('Ticks',[0 0.2 0.4 0.6 0.8 1],'location','NorthOutside');
set(gca,'fontsize',26);
set(gca,'linewidth',3);

% plot estimated model frequencies
[haf,hf,hp] = plotUncertainTimeSeries(out1.Ef,out1.Vf);
xlim([0.5,5.5])
ylim([0,1])
xticks([1:5])
xticklabels({'AU','EV','EV+DN','COMP','PT'})
plot([0.5,5.5],[1,1]/5,'r')
set(gca,'TickDir','out');
box off;
set(gca,'fontsize',26);
set(gca,'linewidth',3);
