data = 'GameTesting.txt';
 
fid1 = fopen(data);


loc = 'GameTestingChan.txt';

fid2 = fopen(loc, 'w');


tline = fgets(fid1);

while ischar(tline)
    
    D = str2num(tline);
    
    D = D + ones(size(D));
    
    lenObs = length(D);
        
    format = repmat('%d ', 1, lenObs-1);
    
    fprintf(fid2, [format,'%d\n'], D');
    
    tline = fgets(fid1);              
end

fclose(fid1);
fclose(fid2);
