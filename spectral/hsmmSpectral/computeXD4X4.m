function O2O3 = computeXD4X4(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);


O2O3 = zeros(x, x, x, d);

for ind_d4 = 1:d    
for ind_x4 = 1:x   
   
    
    for ind_d1 = 1:d
    for ind_d2 = 1:d
    for ind_d3 = 1:d    
    for ind_x1 = 1:x
    for ind_x2 = 1:x    
    for ind_x3 = 1:x    

    O2O3(ind_x1, ind_x2, ind_x2, ind_d2) = O2O3(ind_x1, ind_x2, ind_x2, ind_d2) + ...
        A(ind_x4, ind_x3, ind_d3) * D(ind_d4, ind_x4, ind_d3) * A(ind_x3, ind_x2, ind_d2) * D(ind_d3, ind_x3, ind_d2) * D(ind_d2, ind_x2, ind_d1) * A(ind_x2, ind_x1, ind_d1) * D1(ind_d1, ind_x1) * A1(ind_x1);
        
    end
    end
    end
    end
    end
    end    
    
end
end



