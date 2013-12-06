function test = createTesting_real(ID, p)

%directory to read anomaly data from
dir_test = '../../data/lpr-mit-live-anomaly/';

list_test = dir(dir_test);
numSeq_anom = length(list_test)-3;

%the result will be written to 
loc = strcat('../../libdai/examples/data/HSMMtraining_',num2str(ID),'.txt');
fidhsmm = fopen(loc, 'w');

%data needed to visualize system call commands
%names = importdata('lpr-mit-live-callnames.txt');
%out = fopen('commands.txt', 'w');


%% ==================== PREPARE DATA =======================================

%directory to read training data from
dir_norm = '../../data/lpr-mit-live-normal/';

list_norm = dir(dir_norm);
numSeq_norm = length(list_norm)-3;

%for now find indeces of sequences that are short
reduced_data = 0;

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
    num = ceil(numSeq_norm*p);
    ind_norm = numSeq_norm-num+1:numSeq_norm;
end


%test only on a subset of testing data
ind_anom = [1 2 14 15 16 846 847 1001];


fprintf(fidhsmm, '%d\n', length(ind_norm)+length(ind_anom));

test = cell(length(ind_norm)+length(ind_anom), 1);

iterator = 1;

%symblos which occur in the training and testing data
load ../../data/map.mat;

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
    
    dataHSMM = [[2:3:3*(lenObs-1)-1 3*(lenObs-1)+1]; obs'];
    
    format = repmat('%d\t', 1, lenObs-1);
    
    fprintf(fidhsmm, [format,'%d\n'], dataHSMM'); 
    
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
    
    dataHSMM = [[2:3:3*(lenObs-1)-1 3*(lenObs-1)+1]; obs'];
    
    format = repmat('%d\t', 1, lenObs-1);
    
    fprintf(fidhsmm, [format,'%d\n'], dataHSMM'); 
    
    test{iterator} = obs' + ones(size(obs'));    
    iterator = iterator + 1;    
end

fclose(fidhsmm);




