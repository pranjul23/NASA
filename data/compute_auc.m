%script to e get AUC scores for classifier
close all

%ID = [40:59 0:39 60:79];
ID = 0:39;

N = length(ID);

loc = '../libdai/examples/data/results/';

HSMMscores = [];
for i=1:N
    file = strcat(loc, 'HSMMlikelihood_test_', num2str(ID(i)), '.txt');
    HSMMscores = [HSMMscores getScores(file)];
end

% 
% HMMscores = [];
% for i=1:N
%     file = strcat(loc, 'HMMlikelihood_test_', num2str(ID(i)), '.txt');
%     HMMscores = [HMMscores getScores(file)];
% end



%plot results
S = 4;
A = 5;
auc = zeros(A,S);
dev = zeros(A,S);

for i=1:A
    for j=1:S
        auc(i,j) = mean(HSMMscores(i, (j-1)*10+1 : j*10));
        dev(i,j) = std(HSMMscores(i, (j-1)*10+1 : j*10));
    end

    figure;
    errorbar(1:S, auc(i,:), dev(i,:));
    xlabel('Number of hidden states');
    ylabel('AUC with st. dev.');
    ylim([0 1]);
end



