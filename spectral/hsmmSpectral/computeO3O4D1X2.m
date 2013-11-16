function O3O4 = computeO3O4D1X2(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);


O3O4 = zeros(o, o, d, x);

for ind_o3 = 1:o
for ind_o4 = 1:o
for ind_d1 = 1:d    
for ind_x2 = 1:x   
   
    
    for ind_d2 = 1:d
    for ind_d3 = 1:d
    for ind_x3 = 1:x
    for ind_x4 = 1:x

    O3O4(ind_o3, ind_o4, ind_d1, ind_x2) = O3O4(ind_o3, ind_o4, ind_d1, ind_x2) + ...
        O(ind_o3, ind_x3) * O(ind_o4, ind_x4) * A(ind_x4, ind_x3, ind_d3) * A(ind_x3, ind_x2, ind_d2) * D(ind_d3, ind_x3, ind_d2) * D(ind_d2, ind_x2, ind_d1);
        
    end
    end
    end
    end    
    
end
end
end
end



