%compute P(O1,O3)

function O1O3 = computeO1O3(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);

O3X1 = zeros(o,x);

for ind_o3 = 1:o
    for ind_x1 = 1:x
        
        for ind_d1 = 1:d
        for ind_d2 = 1:d
        for ind_x2 = 1:x
        for ind_x3 = 1:x
            
            O3X1(ind_o3, ind_x1) = O3X1(ind_o3, ind_x1) + ...
                O(ind_o3, ind_x3) * A(ind_x3, ind_x2, ind_d2) * D(ind_d2, ind_x2, ind_d1) * A(ind_x2, ind_x1, ind_d1) * D1(ind_d1, ind_x1) * A1(ind_x1);
        end
        end
        end
        end
    end
end


%normalize
O3X1 = O3X1/sum(sum(O3X1));



O1O3 = zeros(o,o);

for ind_o1 = 1:o
    for ind_o3 = 1:o
        
        for ind_x1 = 1:x
            O1O3(ind_o1, ind_o3) = O1O3(ind_o1, ind_o3) + ...
                O(ind_o1, ind_x1) * O3X1(ind_o3, ind_x1);
        end
    end
end
    

%% ========================================================================

O1X3 = zeros(o,x);

for ind_o1 = 1:o
    for ind_x3 = 1:x
        
        for ind_d1 = 1:d
        for ind_d2 = 1:d
        for ind_x1 = 1:x
        for ind_x2 = 1:x
            
            O1X3(ind_o1, ind_x3) = O1X3(ind_o1, ind_x3) + ...
                O(ind_o1, ind_x1) * A(ind_x3, ind_x2, ind_d2) * D(ind_d2, ind_x2, ind_d1) *  A(ind_x2, ind_x1, ind_d1) * D1(ind_d1, ind_x1) * A1(ind_x1);
        end
        end
        end
        end
    end
end

%normalize
O1X3 = O1X3/sum(sum(O1X3));

V = zeros(o,o);

for ind_o1 = 1:o
    for ind_o3 = 1:o
        
        for ind_x3 = 1:x
            V(ind_o1, ind_o3) = V(ind_o1, ind_o3) + ...
                O(ind_o3, ind_x3) * O1X3(ind_o1, ind_x3);
        end
    end
end
