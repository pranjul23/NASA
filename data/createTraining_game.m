function [] = createTraining_game(ID)

data = 'GameTraining.txt';
 
fid = fopen(data);

tline = fgets(fid);
len=1;
while ischar(tline)
    tline = fgets(fid);
    len = len+1;
end
len = len-1;

fclose(fid);

loc = strcat('../libdai/examples/data/HSMMtraining_',num2str(ID),'.txt');
fidhsmm = fopen(loc, 'w');

loc = strcat('../libdai/examples/data/HMMtraining_',num2str(ID),'.txt');
fidhmm = fopen(loc, 'w');

fprintf(fidhsmm, '%d\n',len);
fprintf(fidhmm, '%d\n',len);

%save data for Murphyk HMM function 
train = cell(len,1);

fid = fopen(data); 

tline = fgets(fid);
counter=1;

while ischar(tline)
    
    D = str2num(tline);
    
    lenObs = length(D);
    fprintf(fidhsmm, '%d\n', lenObs);
    fprintf(fidhmm, '%d\n', lenObs);
    
    dataHSMM = [4:3:3*lenObs+1; D];
    dataHMM = [2:2:2*lenObs; D];
    
    format = repmat('%d\t', 1, lenObs-1);
    
    fprintf(fidhsmm, [format,'%d\n'], dataHSMM');
    fprintf(fidhmm, [format,'%d\n'], dataHMM');
    
    tline = fgets(fid);    
      
    %save data for Murphyk HMM function    
    train{counter} = D + ones(size(D));  
    
    counter = counter+1;
end

save('murphykHMMtrainData.mat', 'train');

fclose(fid);
fclose(fidhsmm);
fclose(fidhmm);























    