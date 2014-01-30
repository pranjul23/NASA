function test = createTestingAir_real(ID)

load('dataAnom.mat');

seq_anom_ind = 101:300;
seq_norm_ind = 1:100;

%the result will be written to 
loc = strcat('../../libdai/examples/data/HSMMtesting_',num2str(ID),'.txt');
fidhsmm = fopen(loc, 'w');

%% ========================= NORMAL DATA ==================================

fprintf(fidhsmm, '%d\n', length(seq_norm_ind)+length(seq_anom_ind));

test = cell(length(seq_norm_ind)+length(seq_anom_ind), 1);

iterator = 1;

%first sequences are normal ones
for k=1:length(seq_norm_ind)
    
    i = seq_norm_ind(k);
    
    obs = actions{i};
    lenObs = length(obs);            
    
%     fprintf(fidhsmm, '%d\n', lenObs);   !!!!!!!!! MODIFY THIS      !!!!!!!!!!!!!!!!!!!!!!!!
    
    %remove 1 from all observations (in C++ we assume it starts with 0)
    dataHSMM = [[2:3:3*(lenObs-1)-1 3*(lenObs-1)+1]; obs-1];
    
    format = repmat('%d\t', 1, lenObs-1);
    
%     fprintf(fidhsmm, [format,'%d\n'], dataHSMM');      !!!!!!!!! MODIFY THIS      !!!!!!!!!!!!!!!!!!!!!!!!
    
    test{iterator} = obs;    
    iterator = iterator + 1;
end



%% ====================== ANOMALY DATA ====================================

% load('dataAnom.mat');

%second sequences are anomalous ones
for k=1:length(seq_anom_ind)
    
    i = seq_anom_ind(k);
    
    obs = actions{i};
    
    lenObs = length(obs);            
    
%     fprintf(fidhsmm, '%d\n', lenObs);     !!!!!!!!! MODIFY THIS      !!!!!!!!!!!!!!!!!!!!!!!!
    
    %remove 1 from all observations (in C++ we assume it starts with 0)
    dataHSMM = [[2:3:3*(lenObs-1)-1 3*(lenObs-1)+1]; obs-1];
    
    format = repmat('%d\t', 1, lenObs-1);
    
%     fprintf(fidhsmm, [format,'%d\n'], dataHSMM');     !!!!!!!!! MODIFY THIS      !!!!!!!!!!!!!!!!!!!!!!!!
    
    test{iterator} = obs;    
    iterator = iterator + 1;    
end

fclose(fidhsmm);




