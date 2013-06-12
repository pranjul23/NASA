function [] = createTesting_real(ID, trainInd)

%directory to read anomaly data from
dir_test = 'lpr-mit-live-anomaly/';

list_test = dir(dir_test);
numSeq_anom = length(list_test)-3;

%the result will be written to 
loc = strcat('../libdai/examples/data/HSMMtesting_', num2str(ID), '.txt');
fidhsmm = fopen(loc, 'w');

loc = strcat('../libdai/examples/data/HMMtesting_', num2str(ID), '.txt');
fidhmm = fopen(loc, 'w');

%data needed to visualize system call commands
%names = importdata('lpr-mit-live-callnames.txt');
%out = fopen('commands.txt', 'w');


%% ==================== PREPARE DATA =======================================

%directory to read training data from
dir_norm = 'lpr-mit-live-normal/';

list_norm = dir(dir_norm);
numSeq_norm = length(list_norm)-3;


%for now find indeces of sequences that are short
reduced_data = 1;

if reduced_data   
    ind_norm = [];
    for i=1:numSeq_norm
        obs = load(strcat(dir_norm, list_norm(i+3).name));
        obs = obs(:,2);
        lenObs = length(obs);
        
        %satisfies length requirement and did not appear in training
        if lenObs <= 700 && ~any(trainInd==i)
            ind_norm = [ind_norm i];
        end
    end
    numSeq_norm = length(ind_norm);
else
    ind_norm = 1:numSeq_norm;
end


ind_anom = [1 2 14 15 16 846 847 1001];







fprintf(fidhsmm, '%d\n', numSeq_norm+length(ind_anom));
fprintf(fidhmm, '%d\n', numSeq_norm+length(ind_anom));

%data for Murphyk HMM function
test = cell(length(ind_norm)+length(ind_anom), 1);
iterator = 1;


%symblos which occur in the training and testing data
load map.mat;

%% ========================= NORMAL DATA ==================================

%first sequences are normal ones
for i=1:numSeq_norm
    
    obs = load(strcat(dir_norm, list_norm(ind_norm(i)+3).name));
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

%% ====================== ANOMALY DATA ====================================

%second sequences are anomalous ones
for i=1:length(ind_anom)
    
    obs = load(strcat(dir_test, list_test(ind_anom(i)+3).name));
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


save('murphykHMMtestData.mat', 'test');

fclose(fidhsmm);
fclose(fidhmm);











