function [] = createTraining_game(ID)

data = 'FlightDataGame.txt';
 
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

fprintf(fidhsmm, '%d\n',len-110);
fprintf(fidhmm, '%d\n',len-110);


fid = fopen(data); 

tline = fgets(fid);
counter=1;

while ischar(tline)
    
    if counter > 110
        D = str2num(tline);
        
        lenObs = length(D);
        fprintf(fidhsmm, '%d\n', lenObs);
        fprintf(fidhmm, '%d\n', lenObs);
        
        dataHSMM = [4:3:3*lenObs+1; D];
        dataHMM = [2:2:2*lenObs; D];
        
        format = repmat('%d\t', 1, lenObs-1);
        
        fprintf(fidhsmm, [format,'%d\n'], dataHSMM');
        fprintf(fidhmm, [format,'%d\n'], dataHMM');        
    end
    
    
    tline = fgets(fid);
    counter = counter+1;    
end

fclose(fid);
fclose(fidhsmm);
fclose(fidhmm);




loc = strcat('../libdai/examples/data/HSMMtesting_',num2str(ID),'.txt');
fidhsmm = fopen(loc, 'w');

loc = strcat('../libdai/examples/data/HMMtesting_',num2str(ID),'.txt');
fidhmm = fopen(loc, 'w');


fprintf(fidhsmm, '%d\n',len);
fprintf(fidhmm, '%d\n',len);

fid = fopen(data);
tline = fgets(fid);

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
    counter = counter+1;    
end

fclose(fid);
fclose(fidhsmm);
fclose(fidhmm);





















    