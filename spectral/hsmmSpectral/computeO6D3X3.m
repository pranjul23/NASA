function O6 = computeO6D3X3(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);


O6 = zeros(x, x, x, x, d);

for ind_x6 = 1:x
for ind_d3 = 1:d 
for ind_x3 = 1:x
       
    for ind_d4 = 1:d
    for ind_d5 = 1:d            
    for ind_x4 = 1:x
    for ind_x5 = 1:x
    
    O6(ind_x4, ind_x5, ind_x6, ind_x3, ind_d3) = O6(ind_x4, ind_x5, ind_x6, ind_x3, ind_d3) + ...
         A(ind_x6, ind_x5, ind_d5) * A(ind_x5, ind_x4, ind_d4) * D(ind_d5, ind_x5, ind_d4)* A(ind_x4, ind_x3, ind_d3) * D(ind_d4, ind_x4, ind_d3);
    
    end
    end
    end
    end
    
end
end
end