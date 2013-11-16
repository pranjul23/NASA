function O5O6 = computeO5O6D3X3(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);


O5O6 = zeros(o, o, d, x);

for ind_o6 = 1:o
for ind_o5 = 1:o
for ind_d3 = 1:d 
for ind_x3 = 1:x
       
    for ind_d4 = 1:d
    for ind_d5 = 1:d            
    for ind_x4 = 1:x
    for ind_x5 = 1:x
    for ind_x6 = 1:x    
    
    O5O6(ind_o5, ind_o6, ind_d3, ind_x3) = O5O6(ind_o5, ind_o6, ind_d3, ind_x3) + ...
        O(ind_o5, ind_x5)* O(ind_o6, ind_x6) * A(ind_x6, ind_x5, ind_d5) * A(ind_x5, ind_x4, ind_d4) * D(ind_d5, ind_x5, ind_d4)* A(ind_x4, ind_x3, ind_d3) * D(ind_d4, ind_x4, ind_d3);
    
    end
    end
    end
    end
    end
    
end
end
end
end