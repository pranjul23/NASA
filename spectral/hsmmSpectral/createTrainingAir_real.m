function train = createTrainingAir_real(ID)

load('data.mat');

seq_ind = 1:125;

%the result will be written to 
loc = strcat('../../libdai/examples/data/HSMMtraining_',num2str(ID),'.txt');
fidhsmm = fopen(loc, 'w');

fprintf(fidhsmm, '%d\n', length(seq_ind));

train = cell(length(seq_ind),1);

for k=1:length(seq_ind)
    
    i = seq_ind(k);
    
    obs = actions{i};
    lenObs = length(obs);            
    
    fprintf(fidhsmm, '%d\n', lenObs);
    
    %remove 1 from all observations (in C++ we assume it starts with 0)    
    dataHSMM = [[2:3:3*(lenObs-1)-1 3*(lenObs-1)+1]; obs-1];
    
    format = repmat('%d\t', 1, lenObs-1);
    
    fprintf(fidhsmm, [format,'%d\n'], dataHSMM'); 
    
    train{i} = obs;    
end

fclose(fidhsmm);















