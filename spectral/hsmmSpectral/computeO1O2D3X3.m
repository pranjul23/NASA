function O1O2 = computeO1O2D3X3(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);


O1O2 = zeros(o, o, d, x);

for ind_o1 = 1:o
for ind_o2 = 1:o
for ind_d3 = 1:d    
for ind_x3 = 1:x   
   
    
    for ind_d1 = 1:d
    for ind_d2 = 1:d
    for ind_x1 = 1:x
    for ind_x2 = 1:x

    O1O2(ind_o1, ind_o2, ind_d3, ind_x3) = O1O2(ind_o1, ind_o2, ind_d3, ind_x3) + ...
        O(ind_o1, ind_x1) * O(ind_o2, ind_x2) * A(ind_x3, ind_x2, ind_d2) * D(ind_d3, ind_x3, ind_d2) * D(ind_d2, ind_x2, ind_d1) * A(ind_x2, ind_x1, ind_d1) *D1(ind_d1, ind_x1) * A1(ind_x1);
        
    end
    end
    end
    end    
    
end
end
end
end



