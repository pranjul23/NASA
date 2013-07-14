loc = '../libdai/examples/data/resultChandola/';

ID = 1:29;

N = length(ID);

for i=1:N
    file1 = strcat(loc, 'stideScore', num2str(ID(i)));
    file2 = strcat(loc, 'stideScores', num2str(ID(i)));
    
    fid1 = fopen(file1);
    fid2 = fopen(file2,'w');    
    
    tline = fgets(fid1);
    
    while ischar(tline)
        
        D = str2num(tline);
        
        %remove small numbers
        D(D<1e-6)=1e-6;
        
        windScores = sum(D)/length(D);
        
        tline = fgets(fid1);
                
        fprintf(fid2, '%.5f\n', windScores);
    end
                        
    fclose(fid1);
    fclose(fid2);
            
end


