function O1O3O5O7 = computeO1O3O5O7(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);


O1O3O5O7 = zeros(o, o, o, o);

for ind_o7 = 1:o
for ind_o5 = 1:o
for ind_o3 = 1:o   
for ind_o1 = 1:o   
    
    for ind_d1 = 1:d
    for ind_d2 = 1:d
    for ind_d3 = 1:d    
    for ind_d4 = 1:d
    for ind_d5 = 1:d
    for ind_d6 = 1:d
    for ind_x1 = 1:x    
    for ind_x2 = 1:x    
    for ind_x3 = 1:x        
    for ind_x4 = 1:x
    for ind_x5 = 1:x
    for ind_x6 = 1:x    
    for ind_x7 = 1:x   

    O1O3O5O7(ind_o1, ind_o3, ind_o5, ind_o7) = O1O3O5O7(ind_o1, ind_o3, ind_o5, ind_o7) + ...
        O(ind_o1, ind_x1) * O(ind_o3, ind_x3) * O(ind_o5, ind_x5) * O(ind_o7, ind_x7) *A(ind_x7, ind_x6, ind_d6) *A(ind_x6, ind_x5, ind_d5) * D(ind_d6, ind_x6, ind_d5) * A(ind_x5, ind_x4, ind_d4) * D(ind_d5, ind_x5, ind_d4)* A(ind_x4, ind_x3, ind_d3) * D(ind_d4, ind_x4, ind_d3) * A(ind_x3, ind_x2, ind_d2) * D(ind_d3, ind_x3, ind_d2) * A(ind_x2, ind_x1, ind_d1) * D(ind_d2, ind_x2, ind_d1) *D1(ind_d1, ind_x1) * A1(ind_x1);
    
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end    
    
end
end
end
end