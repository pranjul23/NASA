%script to e get AUC scores for classifier

ID = 0:39;

N = length(ID);

loc = '../libdai/examples/data/';

HSMMscores = [];
for i=1:N
    file = strcat(loc, 'HSMMlikelihood_test_', num2str(ID(i)), '.txt');
    HSMMscores = [HSMMscores getScores(file)];
end


HMMscores = [];
for i=1:N
    file = strcat(loc, 'HMMlikelihood_test_', num2str(ID(i)), '.txt');
    HMMscores = [HMMscores getScores(file)];
end