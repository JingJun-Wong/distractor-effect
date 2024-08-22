clear;clc;close all;

m=0:.01:1;
p=0:.01:1;
[M,P]=ndgrid(m,p);

ticklabels=0:0.2:1;

tcl=tiledlayout(1,3);

nexttile();
title({'Multiplicative (\eta=1)'});
data=M.*P;
imagesc(data);
hold on;
[~,h1]=contour(data,'LineColor','k','LevelList',0.045:0.065:0.955);
h1.LineWidth=2;
xticks = linspace(1, size(data, 2), numel(ticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', ticklabels);
yticks = linspace(1, size(data, 1), numel(ticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', ticklabels)
xlabel({'Probability'})
ylabel({'Magnitude'})
set(gca,'TickDir','out');
box off;
set(gca,'fontsize',26);
set(gca,'linewidth',3);

nexttile();
title({'Composite (\eta=0.5)'});
data=(M.*P+(M+P)/2)/2;
imagesc(data);
hold on;
[~,h2]=contour(data,'LineColor','k','LevelList',0.045:0.065:0.955);
h2.LineWidth=2;
xticks = linspace(1, size(data, 2), numel(ticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', ticklabels);
yticks = linspace(1, size(data, 1), numel(ticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', ticklabels)
xlabel({'Probability'})
set(gca,'TickDir','out');
box off;
set(gca,'fontsize',26);
set(gca,'linewidth',3);

nexttile();
title({'Additive (\eta=0)'});
data=(M+P)/2;
imagesc(data);
hold on;
[~,h3]=contour(data,'LineColor','k','LevelList',0.045:0.065:0.955);
h3.LineWidth=2;
xticks = linspace(1, size(data, 2), numel(ticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', ticklabels);
yticks = linspace(1, size(data, 1), numel(ticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', ticklabels)
xlabel({'Probability'})
set(gca,'TickDir','out');
box off;
set(gca,'fontsize',26);
set(gca,'linewidth',3);

cb=colorbar;
cb.Layout.Tile = 'east';
