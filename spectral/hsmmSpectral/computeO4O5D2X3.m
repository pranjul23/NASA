function O4O5 = computeO4O5D2X3(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);


O4O5 = zeros(o, o, d, x);

for ind_o5 = 1:o
for ind_o4 = 1:o
for ind_d2 = 1:d    
for ind_x3 = 1:x   
   
    
    for ind_d3 = 1:d
    for ind_d4 = 1:d
    for ind_x4 = 1:x
    for ind_x5 = 1:x

    O4O5(ind_o4, ind_o5, ind_d2, ind_x3) = O4O5(ind_o4, ind_o5, ind_d2, ind_x3) + ...
        O(ind_o4, ind_x4) * O(ind_o5, ind_x5) * A(ind_x5, ind_x4, ind_d4) * A(ind_x4, ind_x3, ind_d3) * D(ind_d4, ind_x4, ind_d3) * D(ind_d3, ind_x3, ind_d2);
    
    end
    end
    end
    end    
    
end
end
end
end



