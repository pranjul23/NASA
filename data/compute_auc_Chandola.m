%script to e get AUC scores for classifier
close all
%clear;

ID = 2:19;

N = length(ID);

loc = '../libdai/examples/data/resultChandola/';

auc = [];

for i=1:N
    file = strcat(loc, 'fsaScores', num2str(ID(i)));
    auc = [auc getScores(file)];
end





