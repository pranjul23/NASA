function res = computeO1X2(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);


res = zeros(o, x);

for ind_o1 = 1:o
for ind_x2 = 1:x   
   
    
    for ind_d1 = 1:d
    for ind_x1 = 1:x    
    
    res(ind_o1, ind_x2) = res(ind_o1, ind_x2) + ...
        O(ind_o1, ind_x1) * A(ind_x2, ind_x1, ind_d1) * D1(ind_d1, ind_x1) * A1(ind_x1);
        
    end
    end    
    
end
end



