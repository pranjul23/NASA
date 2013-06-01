function [] = createTraining_real()


%directory to read training data from
dir_train = 'lpr-mit-live-normal/';

list_train = dir(dir_train);
numSeq = length(list_train)-3;

%the result will be written to 
fidhsmm = fopen('../libdai/examples/HSMMtraining.txt', 'w');
fidhmm = fopen('../libdai/examples/HMMtraining.txt', 'w');

%data needed to visualize system call commands
%names = importdata('lpr-mit-live-callnames.txt');
%out = fopen('commands.txt', 'w');


%for now find indeces of sequences that are short
reduced_data = 1;

if reduced_data
    
    ind = [];
    for i=1:numSeq
        obs = load(strcat(dir_train, list_train(i+3).name));
        obs = obs(:,2);
        lenObs = length(obs);
        
        if lenObs <= 160
            ind = [ind i];
        end
    end
    numSeq = length(ind);
else
    ind = 1:numSeq;
end


fprintf(fidhsmm, '%d\n',numSeq);
fprintf(fidhmm, '%d\n',numSeq);


%save data for Murphyk HMM function 
train = cell(numSeq,1);

for i=1:numSeq
    
    obs = load(strcat(dir_train, list_train(ind(i)+3).name));
    obs = obs(:,2);
    
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
    train{i} = obs' + ones(size(obs'));    
end

save('murphykHMMtrainData.mat', 'train');

fclose(fidhsmm);
fclose(fidhmm);















