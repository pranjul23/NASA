function train = createTraining_real(ID, p)


%directory to read training data from
dir_train = '../../data/lpr-mit-live-normal/';

list_train = dir(dir_train);
numSeq = length(list_train)-3;

%the result will be written to 
loc = strcat('../../libdai/examples/data/HSMMtraining_',num2str(ID),'.txt');
fidhsmm = fopen(loc, 'w');

%data needed to visualize system call commands
%names = importdata('lpr-mit-live-callnames.txt');
%out = fopen('commands.txt', 'w');


%for now find indeces of sequences that are short
reduced_data = 0;

%symblos which occur in the training and testing data
load ../../data/map.mat;

if reduced_data   
    ind = [];
    for i=1:numSeq
        obs = load(strcat(dir_train, list_train(i+3).name));
        obs = obs(:,2);
        lenObs = length(obs);
        
        if lenObs <= 300
            ind = [ind i];
        end
    end
    numSeq = length(ind);
else
    numSeq = floor(numSeq*p);
    ind = 1:numSeq;
end

fprintf(fidhsmm, '%d\n',numSeq);


train = cell(numSeq,1);
for i=1:numSeq
    
    obs = load(strcat(dir_train, list_train(ind(i)+3).name));
    obs = obs(:,2);
    
    %remap observations to have only observable system calls
    obs = map(obs+1);
        
%     plot(obs, 'r*')
    
    %sequence = names(obs+ones(size(obs)));
    %fprintf(out, '%s\n', sequence{:});
    
    lenObs = length(obs);            
    
    fprintf(fidhsmm, '%d\n', lenObs);
        
    dataHSMM = [[2:3:3*(lenObs-1)-1 3*(lenObs-1)+1]; obs'];
    
    format = repmat('%d\t', 1, lenObs-1);
    
    fprintf(fidhsmm, [format,'%d\n'], dataHSMM'); 
    
    %add 1 to all observations (in Matlab we assume it starts with 1)
    train{i} = obs' + ones(size(obs'));    
end

fclose(fidhsmm);















