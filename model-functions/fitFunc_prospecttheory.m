function output = fitFunc_prospecttheory(attribute, p_data, Nfit, distort_flag, DN_flag, parmEst)

obFunc = @(x) LL_func(attribute,p_data,DN_flag,x(1),x(2),x(3));

if isempty(parmEst)
    
    T=[1e-5, 5];
    alpha=[0.5 1.5];
    gamma=[0.5 1.5];
    
    T_lim=[0, 1e3];
    alpha_lim=[0, 1e3];
    gamma_lim=[0, 1e3];
    
    B = [T; alpha; gamma];
    B_lim = [T_lim; alpha_lim; gamma_lim];

    LB = B_lim(:,1); UB = B_lim(:,2);

    % grid search starting points:
    Nall = 100;
    X0 = zeros(Nall,size(B,1));
    for i = 1:size(B,1)
        a = B(i,1); b = B(i,2);
        X0(:,i) = a + (b-a).*rand(Nall,1);
    end
    Np = sum(std(X0)~=0); % number of free parameters
    feval = [5000, 5000]; % max number of function evaluations and iterations
    options = optimset('MaxFunEvals',feval(1),'MaxIter',feval(2),'TolFun',1e-20,'TolX',1e-20,'Display','none');

    X0_valid = [];
    for iter = 1:Nall
        init_fval = obFunc(X0(iter,:));
        if isreal(init_fval) && ~isnan(init_fval) && ~isinf(init_fval)
            X0_valid = [X0_valid; X0(iter,:)];
        end
    end

    tic
    parfor iter = 1:Nfit
        [Xfit_grid(iter,:), NegLL_grid(iter)] = fmincon(obFunc,X0_valid(iter,:),[],[],[],[],LB,UB,[],options);
    end
    toc

    [~,best] = min(NegLL_grid);
    Xfit = Xfit_grid(best,:);
    NegLL = NegLL_grid(best);
    LL = -NegLL;

    [~,p_pred,relacc] = obFunc(Xfit);

    n = sum(~isnan(relacc));
    BIC = Np*log(n) + 2*NegLL;
    AIC = Np*2 + 2*NegLL;
    AICc = 2*NegLL + 2*Np + 2*Np*(Np+1)/(n-Np-1);

    output.Xfit = Xfit;
    output.Xfit_grid = Xfit_grid;
    output.NegLL_grid = NegLL_grid;
    output.LL = LL;
    output.BIC = BIC;
    output.AIC = AIC;
    output.AICc = AICc;
    output.pout = p_pred;
    output.relacc = relacc;
else
    output = -obFunc(parmEst);
end

end

%%
% log-likelihood function:
function [negLL, p_pred, rel_acc] = LL_func(attribute,p_data,DN_flag,T,alpha,gamma)

miss = isnan(p_data(:,1));

att_P = attribute(:,[1,2]);
att_X = attribute(:,[3,4]);

sub_P = exp((-(-log(att_P)).^alpha));
sub_X = att_X.^gamma;

% Utility:
compU = sub_X.*sub_P;

v = compU*T;
v = bsxfun(@minus, v, prctile(v,100,2));
p_pred = exp(v) ./ nansum(exp(v),2);

L = p_pred.^p_data(:,1:2);
L_valid = L;
L_valid(miss,:) = [];
L_valid(L_valid==0) = eps;
negLL = -sum(nansum(log(L_valid),2));
if isnan(negLL)
    negLL = 1e10;
end

rel_acc = p_pred(:,1)./sum(p_pred(:,1:2),2);
rel_acc(isnan(p_data(:,1))) = NaN;

end


