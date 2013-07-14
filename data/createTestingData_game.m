function [] = createTestingData_game(ind, ID)

%data = 'PranjulGame/GameTesting.txt';
dataAnom = 'IgorGame/Anomalies.txt'; 
dataNorm = 'IgorGame/Normal.txt';

fid = fopen(dataAnom);

tline = fgets(fid);
len=1;
while ischar(tline)
    tline = fgets(fid);
    len = len+1;
end
len = len-1;

fclose(fid);


%total number of testing sequences
%includes anomalies and normal points
len = len + length(ind);

loc = strcat('../libdai/examples/data/HSMMtesting_',num2str(ID),'.txt');
fidhsmm = fopen(loc, 'w');

loc = strcat('../libdai/examples/data/HMMtesting_',num2str(ID),'.txt');
fidhmm = fopen(loc, 'w');


fprintf(fidhsmm, '%d\n',len);
fprintf(fidhmm, '%d\n',len);


%% ====== Anomaly data ========
fidAnom = fopen(dataAnom);
tline = fgets(fidAnom);

iterator = 1;

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
    
    tline = fgets(fidAnom);
    
    iterator = iterator + 1;
end

%% ====== Normal data ========
fidNorm = fopen(dataNorm);
tline = fgets(fidNorm);

iterator = 1;

while ischar(tline)
    
    %keep only sequences that are in "ind"
    if all(iterator ~= ind)            
        tline = fgets(fidNorm);    
        iterator = iterator + 1;
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
    
    tline = fgets(fidNorm);
    
    iterator = iterator + 1;
end


fclose(fidAnom);
fclose(fidNorm);
fclose(fidhsmm);
fclose(fidhmm);

















    