function [] = createTrainingData_game(ind, ID)

%data = 'PranjulGame/GameTraining.txt';
data = 'IgorGame/Normal.txt';

len = length(ind);

loc = strcat('../libdai/examples/data/HSMMtraining_',num2str(ID),'.txt');
fidhsmm = fopen(loc, 'w');

loc = strcat('../libdai/examples/data/HMMtraining_',num2str(ID),'.txt');
fidhmm = fopen(loc, 'w');

fprintf(fidhsmm, '%d\n',len);
fprintf(fidhmm, '%d\n',len);

fid = fopen(data); 

tline = fgets(fid);
counter=1;

while ischar(tline)
    
    %keep only sequences that are in "ind"
    if all(counter ~= ind)
        tline = fgets(fid);          
        counter = counter+1;
        continue;
    end
    
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
      
    counter = counter+1;
end

fclose(fid);
fclose(fidhsmm);
fclose(fidhmm);























    