%Murphyk HMM
clear;

load murphykHMMtrainData.mat;
load murphykHMMinit.mat;

[LL, A0_post, A_post, O_post] = dhmm_em(train, A0, A', O', 'max_iter', 50);


load murphykHMMtestData.mat;

% use model to compute log likelihood
loglik = zeros(size(test,1),1);

for i=1:size(test,1)
    loglik(i) = dhmm_logprob(test{i}, A0_post, A_post, O_post)/length(test{i});
end
