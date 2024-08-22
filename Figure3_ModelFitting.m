% This script modified based on the scripts provided by Cao & Tsetsos, 2022
% https://github.com/YinanCao/multiattribute-distractor

clear;clc;close all;

addpath(genpath([pwd,'/multiattribute-distractor-main/functions/']))
datadir = [pwd,'/multiattribute-distractor-main/datasets/'];
addpath(genpath([pwd,'/model-functions/'])) % add composite model

datafile = {
       'behav_fmri.mat'
       'gluth_exp4.mat'
       'gluth_exp3.mat'
       'gluth_exp2_HP.mat'
       'gluth_exp1.mat'
       };

% aggregating data
accuracy = [];trial_type = [];probs = [];rews = [];RT = [];
for whichf = 1:length(datafile)
    D1 = load([datadir,datafile{whichf}]);
    trial_type = cat(2,trial_type,D1.behavior.trial_type);
    probs = cat(2,probs,D1.behavior.probs);
    rews = cat(2,rews,D1.behavior.rews);
    accuracy = cat(2,accuracy,D1.behavior.accuracy);
    RT = cat(2,RT,D1.behavior.RT); % in ms
end
n_subj = length(rews);

Nfit = 10;

for s = 1:n_subj
    
    disp(['subj: ',num2str(s)])
    tt = trial_type{s}; % trial type (2 = distractor trial)
    
    if length(unique(tt))>2
        tt(tt~=0&tt~=10) = -99;
        tt(tt==0) = 2;
        tt(tt==10) = 1;
    end
    
    P = probs{s}; % reward probability (HV, LV, D)
    X = rews{s};  % reward magnitude (HV, LV, D)
    P(P<=0) = nan;
    X(X<=0) = nan;
    % normalise reward P and X to (0,1]
    Pnorm = bsxfun(@times, P, 1./prctile(P,100));
    Xnorm = bsxfun(@times, X, 1./prctile(X,100));
    data_acc = accuracy{s}; % p(HV over LV)
    data_rt = RT{s}/1000;
    data_rt(data_rt<0.1) = nan;
    miss =  isnan(data_acc) | isnan(data_rt);
    data_acc(miss) = nan;
    p_data = [data_acc, 1-data_acc]; % data matrix
    attribute = [Pnorm, Xnorm];
    
    distort_flag = 0; % linear (default)
    
    parmEst = [];

    % AU:
    AU.outputFull_B{s} = fitFunc_Binary_static_SI(attribute(tt==1,[1,2,4,5]),p_data(tt==1,:),Nfit,...
        0,distort_flag,parmEst);
    
    % EV:
    DN_flag = 0;
    EV.outputFull_B{s} = fitFunc_Binary_static_EV(attribute(tt==1,[1,2,4,5]),p_data(tt==1,:),Nfit,...
        distort_flag, DN_flag, parmEst);
    
    % EV + DN:
    DN_flag = 1;
    EVDN.outputFull_B{s} = fitFunc_Binary_static_EV(attribute(tt==1,[1,2,4,5]),p_data(tt==1,:),Nfit,...
        distort_flag, DN_flag, parmEst);
    
    % COMP:
    DN_flag = 0;
    COMP.outputFull_B{s} = fitFunc_composite(attribute(tt==1,[1,2,4,5]),p_data(tt==1,:),Nfit,...
        distort_flag, DN_flag, parmEst);
    
    %%%%%%%%%%%%%%%%%%%
    % cross-validation:
    Y = [data_acc,1-data_acc]; % data matrix
    Y = Y(tt==1,:); % b trials
    ntrl = size(Y,1);
    attribute = [Pnorm(tt==1,1:2),Xnorm(tt==1,1:2)]; % attribute matrix
    nfold = 5;
    CV_trl = reshape(1:ntrl,ntrl/nfold,nfold);
    
    % cross-validation:
    for fold = 1:nfold % loop over 5 folds
    trltest = CV_trl(:,fold); trltrain = setdiff(1:ntrl,trltest);

    SI_flag = 0; % AU
    output = fitFunc_Binary_static_SI(attribute(trltrain,:),Y(trltrain,1:2),Nfit,SI_flag,distort_flag,[]);
    parmEst = output.Xfit; % get best parameters
    % cross-fit test data:
    AU.testLL_B(s,1,fold) = fitFunc_Binary_static_SI(attribute(trltest,:),Y(trltest,1:2),Nfit,SI_flag,distort_flag,parmEst);
    
    DN_flag = 0; % EV
    output = fitFunc_Binary_static_EV(attribute(trltrain,:),Y(trltrain,1:2),Nfit,distort_flag,DN_flag,[]);
    parmEst = output.Xfit; % get best parameters
    EV.testLL_B(s,1,fold) = fitFunc_Binary_static_EV(attribute(trltest,:),Y(trltest,1:2),Nfit,...
        distort_flag,DN_flag,parmEst);
    
    DN_flag = 1; % EV + DN
    output = fitFunc_Binary_static_EV(attribute(trltrain,:),Y(trltrain,1:2),Nfit,distort_flag,DN_flag,[]);
    parmEst = output.Xfit; % get best parameters
    EVDN.testLL_B(s,1,fold) = fitFunc_Binary_static_EV(attribute(trltest,:),Y(trltest,1:2),Nfit,...
        distort_flag,DN_flag,parmEst);
    
    DN_flag = 0; % COMP
    output = fitFunc_composite(attribute(trltrain,:),Y(trltrain,1:2),Nfit,distort_flag,DN_flag,[]);
    parmEst = output.Xfit; % get best parameters
    COMP.testLL_B(s,1,fold) = fitFunc_composite(attribute(trltest,:),Y(trltest,1:2),Nfit,...
        distort_flag,DN_flag,parmEst);

    
    end
    %%%%%%%%%%%%%%%%%%%
    
end

time_str = strrep(mat2str(fix(clock)),' ','_');
save(['Figure3_Fits_',time_str,'.mat'],'AU','EV','EVDN','COMP')
disp('done')
