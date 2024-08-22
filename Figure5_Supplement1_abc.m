% This script modified based on the scripts provided by Cao & Tsetsos, 2022
% https://github.com/YinanCao/multiattribute-distractor

clear;clc;close all;

datadir = [pwd,'/multiattribute-distractor-main/datasets/'];
datafile = {
       'behav_fmri.mat'
       'gluth_exp4.mat'
       'gluth_exp3.mat'
       'gluth_exp2_HP.mat' % gluth high pressure
       'gluth_exp1.mat'
       };

% aggregating data
accuracy = [];trial_type = [];
probs = [];rews = [];RT = [];
for whichf = 1:length(datafile)
    D1 = load([datadir,datafile{whichf}]);
    accuracy = cat(2,accuracy,D1.behavior.accuracy);
    trial_type = cat(2,trial_type,D1.behavior.trial_type);
    probs = cat(2,probs,D1.behavior.probs);
    rews = cat(2,rews,D1.behavior.rews);
    RT = cat(2,RT,D1.behavior.RT); % in ms
end

n_subj = length(rews);

addpath(genpath([pwd,'/model-fits/']))

% get fitted values
fits = load('Figure3_Fits.mat');

int_coef = [];
for i=1:numel(fits.COMP.outputFull_B)
    int_coef = [int_coef fits.COMP.outputFull_B{i}.Xfit(1)];
end

% split data based on integration coefficient
m_ic = nanmean(int_coef);
multi_idx = find(int_coef>m_ic);
add_idx = find(int_coef<m_ic);

for s = multi_idx
    
    % trial type (2 = ternary trial, 1 = binary)
    % 1 = binary
    tt = trial_type{s};
    if length(unique(tt))>2 % some datasets have irrelevant conditions
        tt(tt~=0&tt~=10) = nan;
        tt(tt==0) = 2;
        tt(tt==10) = 1;
    end
    P = probs{s}; % reward probability (H, L, D)
    X = rews{s};  % reward magnitude (H, L, D)
    P(P<=0) = nan;
    X(X<=0) = nan;
    % rescale to (0,1)
    Pnorm = bsxfun(@times, P, 1./prctile(P,100));
    Xnorm = bsxfun(@times, X, 1./prctile(X,100));
    
    data_acc = accuracy{s}; % p(H over L) relative accuracy

    % for regression:
    Ttrl = [Pnorm(tt==2,:),Xnorm(tt==2,:),data_acc(tt==2,1)]; % ternary,last col = accuacy
    Btrl = [Pnorm(tt==1,:),Xnorm(tt==1,:),data_acc(tt==1,1)]; % binary

    ntrl = size(Ttrl,1);
    regreMatrix = []; % now, let's prepare stuff for the regression
    for trl = ntrl:-1:1 % each T trial, find matched B trials
        bi_id = find(sum(abs(bsxfun(@minus,Btrl(:,[1,2,4,5]),Ttrl(trl,[1,2,4,5]))),2)<1e-10);
        Bmatch = Btrl(bi_id,end); % B responses, 1 or 0
        if ~isempty(Bmatch(~isnan(Bmatch)))
            Hchoice = nansum(Bmatch); % n of H choice
            Ntot = nansum(~isnan(Bmatch)); % total N observations
        else % only 1 trl, but missed response
            Hchoice = nan;
            Ntot = nan;
        end
        Ntot_ternary = 1;
        regreMatrix(trl,:) = [Ttrl(trl,:),Ntot_ternary,Hchoice,Ntot];
    end

    % choice response:
    % two-column for glmfit: number of successes in the corresponding number of trials in n
    Ty = regreMatrix(:,7:8);  % two-column input, ternary
    By = regreMatrix(:,9:10); % two-column input, binary
    
    % prepare GLM regressors
    % regreMatrix: first 6 cols are: ph, pl, pd, xh, xl, xd
    
    hv = regreMatrix(:,1).*regreMatrix(:,4); % get EV
    lv = regreMatrix(:,2).*regreMatrix(:,5);
    dv = regreMatrix(:,3).*regreMatrix(:,6);

    % regression for T and B individually
    regX = zscore([hv+lv hv-lv dv-hv zscore(hv-lv).*zscore(dv-hv)]); % GLM regressors
    
    linkFun = 'logit';
    
    % regression for T & B combined
    task = [zeros(size(regX,1),1);ones(size(regX,1),1)]; % D present
    Xcomb = [repmat(regX,2,1),task,bsxfun(@times,repmat(regX,2,1),task)];
    Xcomb = Xcomb(:,[1:5 7:9]);
    ycomb = [By;Ty];
    multi_betas(s,:) = glmfit(Xcomb,ycomb,'binomial',linkFun); % T & B combined
end

for s = add_idx
    
    % trial type (2 = ternary trial, 1 = binary)
    % 1 = binary
    tt = trial_type{s};
    if length(unique(tt))>2 % some datasets have irrelevant conditions
        tt(tt~=0&tt~=10) = nan;
        tt(tt==0) = 2;
        tt(tt==10) = 1;
    end
    P = probs{s}; % reward probability (H, L, D)
    X = rews{s};  % reward magnitude (H, L, D)
    P(P<=0) = nan;
    X(X<=0) = nan;
    % rescale to (0,1)
    Pnorm = bsxfun(@times, P, 1./prctile(P,100));
    Xnorm = bsxfun(@times, X, 1./prctile(X,100));
    
    data_acc = accuracy{s}; % p(H over L) relative accuracy

    % for regression:
    Ttrl = [Pnorm(tt==2,:),Xnorm(tt==2,:),data_acc(tt==2,1)]; % ternary,last col = accuacy
    Btrl = [Pnorm(tt==1,:),Xnorm(tt==1,:),data_acc(tt==1,1)]; % binary

    ntrl = size(Ttrl,1);
    regreMatrix = []; % now, let's prepare stuff for the regression
    for trl = ntrl:-1:1 % each T trial, find matched B trials
        bi_id = find(sum(abs(bsxfun(@minus,Btrl(:,[1,2,4,5]),Ttrl(trl,[1,2,4,5]))),2)<1e-10);
        Bmatch = Btrl(bi_id,end); % B responses, 1 or 0
        if ~isempty(Bmatch(~isnan(Bmatch)))
            Hchoice = nansum(Bmatch); % n of H choice
            Ntot = nansum(~isnan(Bmatch)); % total N observations
        else % only 1 trl, but missed response
            Hchoice = nan;
            Ntot = nan;
        end
        Ntot_ternary = 1;
        regreMatrix(trl,:) = [Ttrl(trl,:),Ntot_ternary,Hchoice,Ntot];
    end

    % choice response:
    % two-column for glmfit: number of successes in the corresponding number of trials in n
    Ty = regreMatrix(:,7:8);  % two-column input, ternary
    By = regreMatrix(:,9:10); % two-column input, binary
    
    % prepare GLM regressors
    % regreMatrix: first 6 cols are: ph, pl, pd, xh, xl, xd
    
    hv = regreMatrix(:,1).*regreMatrix(:,4); % get EV
    lv = regreMatrix(:,2).*regreMatrix(:,5);
    dv = regreMatrix(:,3).*regreMatrix(:,6);

    % regression for T and B individually
    regX = zscore([hv+lv hv-lv dv-hv zscore(hv-lv).*zscore(dv-hv)]); % GLM regressors
    
    linkFun = 'logit';
    
    % regression for T & B combined
    task = [zeros(size(regX,1),1);ones(size(regX,1),1)]; % D present
    Xcomb = [repmat(regX,2,1),task,bsxfun(@times,repmat(regX,2,1),task)];
    Xcomb = Xcomb(:,[1:5 7:9]);
    ycomb = [By;Ty];
    add_betas(s,:) = glmfit(Xcomb,ycomb,'binomial',linkFun); % T & B combined
end

% set exclusion criteria to exclude subjects for overfitting
exclude_criteria=10;
excl1=(abs(multi_betas)>exclude_criteria);
ind1=find(sum(excl1,2)>0 | sum(multi_betas,2)==0);
multi_betas(ind1,:)=NaN;
excl2=(abs(add_betas)>exclude_criteria);
ind2=find(sum(excl2,2)>0 | sum(add_betas,2)==0);
add_betas(ind2,:)=NaN;

multi_betas = multi_betas(multi_idx,:);
add_betas = add_betas(add_idx,:);

multi_beta_mean = nanmean(multi_betas);
add_beta_mean = nanmean(add_betas);
multi_beta_sd = nanstd(multi_betas);
add_beta_sd = nanstd(add_betas);
multi_beta_se = multi_beta_sd./sqrt(size(multi_betas,1));
add_beta_se = add_beta_sd./sqrt(size(add_betas,1));

[~,multi_add_ttest_p,~,multi_add_ttest_stats] = ttest2(multi_betas(:,end-1),add_betas(:,end-1))
[~,multi_ttest_p,~,multi_ttest_stats] = ttest(multi_betas(:,2:end))
[~,add_ttest_p,~,add_ttest_stats] = ttest(add_betas(:,2:end))

% plot graphs
figure;
y1=[multi_beta_mean(end-1) add_beta_mean(end-1)];
x1=1:length(y1);
se1=[multi_beta_se(end-1) add_beta_se(end-1)];
b=bar(x1,y1);
b.FaceColor = 'flat';
b.CData(1,:) = [100 143 255]/256;
b.CData(2,:) = [255 176 0]/256;
b.BarWidth = 0.5;
b.LineWidth = 2;
hold on;
errorbar(x1,y1,se1,'k','LineStyle','none','linewidth',2,'CapSize',18);
xlim([0.5,2.5])
ylim([-0.15,0.15])
xticks([1:8])
yticks([-0.15:0.05:0.15])
xticklabels({'Multiplicative','Additive'})
ylabel({'Effect size of (D-HV)T on Accuracy';'(a.u.)'})
set(gca,'TickDir','out');
box off;
set(gca,'fontsize',26);
set(gca,'linewidth',3);

figure;
y2=[multi_beta_mean(2:end)];
x2=1:length(y2);
se2=[multi_beta_se(2:end)];
b=bar(x2,y2);
b.FaceColor = 'flat';
for i=1:length(y2)
    b.CData(i,:)=[0.5 0.5 0.5];
end
b.CData(length(y2)-1,:) = [100 143 255]/256;
b.BarWidth = 0.5;
b.LineWidth = 2;
hold on;
errorbar(x2,y2,se2,'k','LineStyle','none','linewidth',2,'CapSize',18);
xlim([0.5,length(y2)+0.5])
ylim([-0.3,0.8])
xticks([1:length(y2)])
xticklabels({'HV+LV','HV-LV','DV-HV','(HV-LV)(DV-HV)','T','(HV-LV)T','(DV-HV)T','(HV-LV)(DV-HV)T'})
xtickangle(45)
ylabel({'Effect size on Accuracy';'(a.u.)'})
set(gca,'TickDir','out');
title('Multiplicative Group')
box off;
set(gca,'fontsize',26);
set(gca,'linewidth',3);

figure;
y3=[add_beta_mean(2:end)];
x3=1:length(y3);
se3=[add_beta_se(2:end)];
b=bar(x3,y3);
b.FaceColor = 'flat';
for i=1:length(y3)
    b.CData(i,:)=[0.5 0.5 0.5];
end
b.CData(length(y3)-1,:) = [255 176 0]/256;
b.BarWidth = 0.5;
b.LineWidth = 2;
hold on;
errorbar(x3,y3,se3,'k','LineStyle','none','linewidth',2,'CapSize',18);
xlim([0.5,length(y3)+0.5])
ylim([-0.3,0.8])
xticks([1:length(y3)])
xticklabels({'HV+LV','HV-LV','DV-HV','(HV-LV)(DV-HV)','T','(HV-LV)T','(DV-HV)T','(HV-LV)(DV-HV)T'})
xtickangle(45)
ylabel({'Effect size on Accuracy';'(a.u.)'})
set(gca,'TickDir','out');
title('Additive Group')
box off;
set(gca,'fontsize',26);
set(gca,'linewidth',3);
