function O2O4O7O9 = computeO2O4O7O9(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);


O2O4O7O9 = zeros(o, o, o, o);

% L = kron(kron(kron(O,O),O),O);


for ind_o9 = 1:o
for ind_o7 = 1:o
for ind_o4 = 1:o
for ind_o2 = 1:o   
   
    
    for ind_d1 = 1:d
    for ind_d2 = 1:d
    for ind_d3 = 1:d    
    for ind_d4 = 1:d
    for ind_d5 = 1:d
    for ind_d6 = 1:d
    for ind_d7 = 1:d    
    for ind_d8 = 1:d    
    for ind_x1 = 1:x    
    for ind_x2 = 1:x    
    for ind_x3 = 1:x        
    for ind_x4 = 1:x
    for ind_x5 = 1:x
    for ind_x6 = 1:x    
    for ind_x7 = 1:x        
    for ind_x8 = 1:x
    for ind_x9 = 1:x    

        
    O2O4O7O9(ind_o2, ind_o4, ind_o7, ind_o9) = O2O4O7O9(ind_o2, ind_o4, ind_o7, ind_o9) + ...
        O(ind_o2, ind_x2) * O(ind_o4, ind_x4) * O(ind_o7, ind_x7) * O(ind_o9, ind_x9) *A(ind_x9, ind_x8, ind_d8)*A(ind_x8, ind_x7, ind_d7) * D(ind_d8, ind_x8, ind_d7)*A(ind_x7, ind_x6, ind_d6) * D(ind_d7, ind_x7, ind_d6)*A(ind_x6, ind_x5, ind_d5) * D(ind_d6, ind_x6, ind_d5) * A(ind_x5, ind_x4, ind_d4) * D(ind_d5, ind_x5, ind_d4)* A(ind_x4, ind_x3, ind_d3) * D(ind_d4, ind_x4, ind_d3) * A(ind_x3, ind_x2, ind_d2) * D(ind_d3, ind_x3, ind_d2) * A(ind_x2, ind_x1, ind_d1) * D(ind_d2, ind_x2, ind_d1) *D1(ind_d1, ind_x1) * A1(ind_x1);
%         L(ind_o2+o*(ind_o4-1)+o*o*(ind_o7-1)+o*o*o*(ind_o9-1), ind_x2+x*(ind_x4-1)+x*x*(ind_x7-1)+x*x*x*(ind_x9-1))*A(ind_x9, ind_x8, ind_d8)*A(ind_x8, ind_x7, ind_d7) * D(ind_d8, ind_x8, ind_d7)*A(ind_x7, ind_x6, ind_d6) * D(ind_d7, ind_x7, ind_d6)*A(ind_x6, ind_x5, ind_d5) * D(ind_d6, ind_x6, ind_d5) * A(ind_x5, ind_x4, ind_d4) * D(ind_d5, ind_x5, ind_d4)* A(ind_x4, ind_x3, ind_d3) * D(ind_d4, ind_x4, ind_d3) * A(ind_x3, ind_x2, ind_d2) * D(ind_d3, ind_x3, ind_d2) * A(ind_x2, ind_x1, ind_d1) * D(ind_d2, ind_x2, ind_d1) *D1(ind_d1, ind_x1) * A1(ind_x1);
        
    
    
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
    
    
    
end
end
end
end