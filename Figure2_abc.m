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

for s = 1:n_subj
    
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
    regX = zscore([hv-lv dv-hv zscore(hv-lv).*zscore(dv-hv)]); % GLM regressors
    
    linkFun = 'logit';
    betas{1}(s,:) = glmfit(regX,Ty,'binomial',linkFun); % T
    betas{2}(s,:) = glmfit(regX,By,'binomial',linkFun); % B
    
    % regression for T & B combined
    task = [zeros(size(regX,1),1);ones(size(regX,1),1)]; % D present
    Xcomb = [repmat(regX,2,1),task,bsxfun(@times,repmat(regX,2,1),task)];
    ycomb = [By;Ty];
    betas{3}(s,:) = glmfit(Xcomb,ycomb,'binomial',linkFun); % T & B combined
    
    r(:,:,s)=corrcoef(regX);
end


% set exclusion criteria to exclude subjects for overfitting
exclude_criteria=10;

for x=1:3
    excl=(abs(betas{x})>exclude_criteria);
    ind=find(sum(excl,2)>0 | sum(betas{x},2)==0);
    betas{x}(ind,:)=NaN;
end

% statistics and plots
title_labels=[{'Distractor Trials'} {'Two-Option Trials'} {'All Trials'}];

for panel=1:3
    figure;
    curval=betas{panel};

    [~,testp,~,teststats]=ttest(curval);

    beta=nanmean(curval);
    sd=nanstd(curval);
    se=sd./sqrt(n_subj);
    df=teststats.df;
    p_val=testp;
    t_stat=teststats.tstat;

    % statistics table including intercept
    results=table(beta',sd',t_stat',df',p_val')

    beta(1)=[];
    se(1)=[];

    x=1:length(beta);
    y=beta;

    bar(x,y,'FaceColor',[0.5 0.5 0.5],'barwidth',0.5,'linewidth',2);
    hold on;
    errorbar(x,y,se,'k','LineStyle','none','linewidth',2,'CapSize',18);
    
    if panel~=3
        xlim([0.5,3.5])
        xticks([1:3])
        xticklabels({'HV-LV','DV-HV','(HV-LV)(DV-HV)'})

    else
        xlim([0.5,7.5])
        xticks([1:7])
        xticklabels({'HV-LV','DV-HV','(HV-LV)(DV-HV)','T','(HV-LV)T','(DV-HV)T','(HV-LV)(DV-HV)T'})
    end

    xtickangle(45)
    ylim([-0.2,0.8])
    ylabel({'Effect size on Accuracy';'(a.u.)'})
    set(gca,'TickDir','out');
    title(title_labels{panel});
    box off;
    set(gca,'fontsize',26);
    set(gca,'linewidth',3);
end