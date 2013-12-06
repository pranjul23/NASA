function res = computeO2X1(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);


res = zeros(o, x);

for ind_o2 = 1:o
for ind_x1 = 1:x   
   
    
    for ind_d1 = 1:d
    for ind_x2 = 1:x    
    
    res(ind_o2, ind_x1) = res(ind_o2, ind_x1) + ...
        O(ind_o2, ind_x2) * A(ind_x2, ind_x1, ind_d1) * D1(ind_d1, ind_x1);
        
    end
    end    
    
end
end



