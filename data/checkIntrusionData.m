
%directory to read training data from
dir_train = 'lpr-mit-live-normal/';

list_train = dir(dir_train);
numSeqTrain = length(list_train)-3;


%directory to read anomaly data from
dir_test = 'lpr-mit-live-anomaly/';

list_test = dir(dir_test);
numSeqTest = length(list_test)-3;

A = [];

for i=1:numSeqTrain
    obs = load(strcat(dir_train, list_train(i+3).name));
    obs = obs(:,2);
    
    A = union(A,unique(obs));    
end


for i=1:numSeqTest
    obs = load(strcat(dir_test, list_test(i+3).name));
    obs = obs(:,2);
    
    A = union(A,unique(obs));    
end


A = A + 1;
map = (-1)*ones(183,1);

for i=1:length(A)
   map(A(i)) = i-1;
end