function O6O7 = computeO6O7D4X4(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);


O6O7 = zeros(o, o, d, x);

for ind_o7 = 1:o
for ind_o6 = 1:o
for ind_d4 = 1:d   
for ind_x4 = 1:x   
    
    for ind_d5 = 1:d
    for ind_d6 = 1:d
    for ind_x5 = 1:x
    for ind_x6 = 1:x    
    for ind_x7 = 1:x   

    O6O7(ind_o6, ind_o7, ind_d4, ind_x4) = O6O7(ind_o6, ind_o7, ind_d4, ind_x4) + ...
        O(ind_o6, ind_x6) * O(ind_o7, ind_x7) *A(ind_x7, ind_x6, ind_d6) *A(ind_x6, ind_x5, ind_d5) * D(ind_d6, ind_x6, ind_d5) * A(ind_x5, ind_x4, ind_d4) * D(ind_d5, ind_x5, ind_d4);
    
    end
    end
    end
    end
    end
    
end
end
end
end