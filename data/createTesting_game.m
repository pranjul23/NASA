function [] = createTesting_game(ID)

data = 'GameTesting.txt';
 
fid = fopen(data);

tline = fgets(fid);
len=1;
while ischar(tline)
    tline = fgets(fid);
    len = len+1;
end
len = len-1;

fclose(fid);


loc = strcat('../libdai/examples/data/HSMMtesting_',num2str(ID),'.txt');
fidhsmm = fopen(loc, 'w');

loc = strcat('../libdai/examples/data/HMMtesting_',num2str(ID),'.txt');
fidhmm = fopen(loc, 'w');


fprintf(fidhsmm, '%d\n',len);
fprintf(fidhmm, '%d\n',len);

%data for Murphyk HMM function
test = cell(len, 1);


fid = fopen(data);
tline = fgets(fid);

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
    
    tline = fgets(fid);
    
    %save data for Murphyk HMM function    
    test{iterator} = D + ones(size(D));  
    iterator = iterator + 1;
end

save('murphykHMMtestData.mat', 'test');

fclose(fid);
fclose(fidhsmm);
fclose(fidhmm);

















    