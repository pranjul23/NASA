function [] = createTesting_real(ID)

%directory to read anomaly data from
dir_test = 'lpr-mit-live-anomaly/';

list_test = dir(dir_test);
numSeq = length(list_test)-3;

%the result will be written to 
loc = strcat('../libdai/examples/data/HSMMtesting_', num2str(ID), '.txt');
fidhsmm = fopen(loc, 'w');

loc = strcat('../libdai/examples/data/HMMtesting_', num2str(ID), '.txt');
fidhmm = fopen(loc, 'w');

%data needed to visualize system call commands
%names = importdata('lpr-mit-live-callnames.txt');
%out = fopen('commands.txt', 'w');

%numSeq = 20;

fprintf(fidhsmm, '%d\n', numSeq+100);
fprintf(fidhmm, '%d\n', numSeq+100);

%data for Murphyk HMM function
test = cell(numSeq+100, 1);
iterator = 1;


%directory to read training data from
dir_train = 'lpr-mit-live-normal/';

list_train = dir(dir_train);


%symblos which occur in the training and testing data
load map.mat;


%first 10 sequences are normal ones
for i=150:250
    
    obs = load(strcat(dir_train, list_train(i+3).name));
    obs = obs(:,2);
    
    %remap observations to have only observable system calls
    obs = map(obs+1);
    
    %plot(obs, 'r*')
    
    %sequence = names(obs+ones(size(obs)));
    %fprintf(out, '%s\n', sequence{:});
    
    lenObs = length(obs);            
    
    fprintf(fidhsmm, '%d\n', lenObs);
    fprintf(fidhmm, '%d\n', lenObs);
    
    dataHSMM = [4:3:3*lenObs+1; obs'];
    dataHMM = [2:2:2*lenObs; obs'];
    
    format = repmat('%d\t', 1, lenObs-1);
    
    fprintf(fidhsmm, [format,'%d\n'], dataHSMM'); 
    fprintf(fidhmm, [format,'%d\n'], dataHMM');
    
    %save data for Murphyk HMM function    
    test{iterator} = obs' + ones(size(obs'));    
    iterator = iterator + 1;
end


%second 10 sequences are anomalous ones
%take 5 from top
for i=1:numSeq
    
    obs = load(strcat(dir_test, list_test(i+3).name));
    obs = obs(:,2);
    
    %remap observations to have only observable system calls
    obs = map(obs+1);
    
    %plot(obs, 'r*')
    
    %sequence = names(obs+ones(size(obs)));
    %fprintf(out, '%s\n', sequence{:});
    
    lenObs = length(obs);            
    
    fprintf(fidhsmm, '%d\n', lenObs);
    fprintf(fidhmm, '%d\n', lenObs);
    
    dataHSMM = [4:3:3*lenObs+1; obs'];
    dataHMM = [2:2:2*lenObs; obs'];
    
    format = repmat('%d\t', 1, lenObs-1);
    
    fprintf(fidhsmm, [format,'%d\n'], dataHSMM'); 
    fprintf(fidhmm, [format,'%d\n'], dataHMM');
    
    %save data for Murphyk HMM function    
    test{iterator} = obs' + ones(size(obs'));    
    iterator = iterator + 1;    
end


% %take 5 from bottom
% for i = numSeq:-1:numSeq-4
%     
%     obs = load(strcat(dir_test, list_test(i+3).name));
%     obs = obs(:,2);
%     
%     %remap observations to have only observable system calls
%     obs = map(obs+1);
%     
%     %plot(obs, 'r*')
%     
%     %sequence = names(obs+ones(size(obs)));
%     %fprintf(out, '%s\n', sequence{:});
%     
%     lenObs = length(obs);            
%     
%     fprintf(fidhsmm, '%d\n', lenObs);
%     fprintf(fidhmm, '%d\n', lenObs);
%     
%     dataHSMM = [4:3:3*lenObs+1; obs'];
%     dataHMM = [2:2:2*lenObs; obs'];
%     
%     format = repmat('%d\t', 1, lenObs-1);
%     
%     fprintf(fidhsmm, [format,'%d\n'], dataHSMM'); 
%     fprintf(fidhmm, [format,'%d\n'], dataHMM');
%     
%     %save data for Murphyk HMM function    
%     test{iterator} = obs' + ones(size(obs'));    
%     iterator = iterator + 1;  
% end

save('murphykHMMtestData.mat', 'test');

fclose(fidhsmm);
fclose(fidhmm);











