%%%% Multiple Testing Problem %%%%
clc;
clear all;
close all;
addpath(genpath('E:\research\toolboxes\general'));

%% Test Scenario Setting
n = 1e1; % Number of hypotheses
n1 = 0; % Number of false null
K = 1e4; % Number of experiments
alpha = 0.05; 
eff = 2;
mu = [eff*ones(1,n1),zeros(1,n-n1)]; % Drawing from Normal distribution
sigma = 1; 

%% No Correction
sen = zeros(K,1);
spec = zeros(K,1);
nFalse = zeros(K,1);
for k = 1:K 
    x = normrnd(mu,sigma);
    p = 2*min(normcdf(x),1-normcdf(x)); % Estimate p-values for two-sided test
    nFalse(k) = sum(p<=alpha & mu==0)/max(1,sum(p<=alpha));
    sen(k) = sum(p(1:n1)<=alpha)/n1; % Proportion of false null with value exceeding threshold
    spec(k) = sum(p(n1+1:end)>alpha)/(n-n1); % Proportion of false null with value exceeding threshold
end
disp('No Correction');
sensitivity = mean(sen)
specificity = mean(spec)
FWER = sum(nFalse>0)/K

%% Bonferroni Correction
sen = zeros(K,1);
spec = zeros(K,1);
nFalse = zeros(K,1);
for k = 1:K 
    x = normrnd(mu,sigma);
    p = 2*min(normcdf(x),1-normcdf(x)); % Estimate p-values for two-sided test
    nFalse(k) = sum(p<=alpha/n & mu==0)/max(1,sum(p<=alpha/n));
    sen(k) = sum(p(1:n1)<=alpha/n)/n1; % Proportion of false null with value exceeding threshold
    spec(k) = sum(p(n1+1:end)>alpha/n)/(n-n1); % Proportion of false null with value exceeding threshold
end
disp('Bonferroni');
sensitivity = mean(sen)
specificity = mean(spec)
FWER = sum(nFalse>0)/K

%% Step-up Procedure
i = (n:-1:1);
sen = zeros(K,1);
spec = zeros(K,1);
nFalse = zeros(K,1);
for k = 1:K 
    x = normrnd(mu,sigma);
    p = 2*min(normcdf(x),1-normcdf(x)); % Estimate p-values for two-sided test
    [pS,indS] = sort(p,2,'descend'); % Sort p-values in descending order
    ind = find(pS<=alpha./(n-i+1),1,'first'); % Index of largest p-value <= alpha/(n-i+1)
    sig = zeros(1,n);
    sig(indS(ind:end)) = 1;
    nFalse(k) = sum(sig & mu==0)/max(1,sum(sig));
    sen(k) = sum(sig & mu~=0)/n1;
    spec(k) = sum(sig==0 & mu==0)/(n-n1);
end
disp('Step Up Procedure');
sensitivity = mean(sen)
specificity = mean(spec)
FWER = sum(nFalse>0)/K

%% FDR Correction
i = (1:n);
sen = zeros(K,1);
spec = zeros(K,1);
nFalse = zeros(K,1);
for k = 1:K 
    x = normrnd(mu,sigma);
    p = 2*min(normcdf(x),1-normcdf(x)); % Estimate p-values for two-sided test
    [pS,indS] = sort(p);
    ind = find(pS<=alpha*i/n,1,'last');
    sig = zeros(1,n);
    sig(indS(1:ind)) = 1;
    nFalse(k) = sum(sig & mu==0)/max(1,sum(sig));
    sen(k) = sum(sig & mu~=0)/n1;
    spec(k) = sum(sig==0 & mu==0)/(n-n1);
end
disp('FDR');
sensitivity = mean(sen)
specificity = mean(spec)
FWER = sum(nFalse>0)/K


%% Max-z Permutation Test
n = 1e1; % Number of hypotheses
n1 = 0; % Number of false null
m = 200; % Number of samples
K = 1e2; % Number of experiments
nPerm = 1e3; % Number of Permutations
alpha = 0.05; 
eff = 2;
mu = [eff*ones(m,n1),zeros(m,n-n1)];
sigma = eye(n); % Independent hypotheses
sen = zeros(K,1);
spec = zeros(K,1);
nFalse = zeros(K,1);
tic;
for k = 1:K 
    disp(['Experiment',int2str(k)]);
    % Generate Data
    X = mvnrnd(mu,sigma); % x drawn from standard normal distribution, i.e. mean x = 0
    [~,~,~,stats] = ttest(X);
    tval = stats.tstat;
    % Estimate Null Distribution
    maxTval = zeros(nPerm,1);
    minTval = zeros(nPerm,1);
    parfor perm = 1:nPerm
        sgn = sign(randn(m,1));
        [~,~,~,stats] = ttest(X.*(sgn*ones(1,n)));
        maxTval(perm) = max(stats.tstat);
        minTval(perm) = min(stats.tstat);
    end
    param = evfit(-maxTval);
    threshPos = -evinv(alpha/2,param(1),param(2));
    param = evfit(minTval);
    threshNeg = evinv(alpha/2,param(1),param(2));
%     threshPos = prctile(maxTval,(1-alpha/2)*100);
%     threshNeg = prctile(minTval,alpha/2*100);
    sig = (tval>=threshPos) | (tval<=threshNeg);
    nFalse(k) = sum(sig & mu(1,:)==0)/max(1,sum(sig));
    sen(k) = sum(sig & mu(1,:)~=0)/n1;
    spec(k) = sum(sig==0 & mu(1,:)==0)/(n-n1);
end
toc;
disp('Max-t Permutation Test');
sensitivity = mean(sen)
specificity = mean(spec)
FWER = sum(nFalse>0)/K % Number of experiments for which one or more hypotheses were declared significant
