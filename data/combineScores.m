loc = '../libdai/examples/data/resultGame4/';

ID = 0:19;

N = length(ID);

for i=1:N
    file1 = strcat(loc, 'HSMMmarginal_test_', num2str(ID(i)), '.txt');
    file2 = strcat(loc, 'HSMMlikelihood_test_', num2str(ID(i)), '.txt');
    
    fid1 = fopen(file1);
    fid2 = fopen(file2,'w');
    
    tline = fgets(fid1);
    
    while ischar(tline)
        
        D = str2num(tline);
        
        %remove small numbers
        %D(D<1e-6)=1e-6;
        
        %remove infinities
        tmp = D;
        tmp(tmp == -inf) = [];
        D(D == -inf) = min(tmp) - 0.5;
        
        windScores = sum(D)/length(D);
        
        tline = fgets(fid1);
        
        fprintf(fid2, '%.5f\n', windScores);
    end
    
    fclose(fid1);
    fclose(fid2);
    
end


