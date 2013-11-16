%compute p(O3,4O4) for HSMM

function O3O4 = computeO3O4(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);

O4X3 = zeros(o,x);

for ind_o4 = 1:o
    for ind_x3 = 1:x
        
        for ind_x1 = 1:x
        for ind_x2 = 1:x
        for ind_x4 = 1:x
        for ind_d1 = 1:d    
        for ind_d2 = 1:d
        for ind_d3 = 1:d
            
                O4X3(ind_o4, ind_x3) = O4X3(ind_o4, ind_x3) + ...
                    O(ind_o4, ind_x4) * A(ind_x4, ind_x3, ind_d3) * D(ind_d3, ind_x3, ind_d2) * D(ind_d2, ind_x2, ind_d1) * A(ind_x2, ind_x1, ind_d1) * D1(ind_d1, ind_x1) *  A1(ind_x1);
        end
        end
        end
        end        
        end
        end
        
    end
end

%normalize
for ind_x3 = 1:x
    O4X3(:, ind_x3) = O4X3(:, ind_x3)/sum(O4X3(:, ind_x3));
end


X3 = zeros(x,1);

for ind_x3 = 1:x
    for ind_x1 = 1:x
    for ind_x2 = 1:x
    for ind_d1 = 1:d
    for ind_d2 = 1:d
        
        X3(ind_x3) = X3(ind_x3) + ...
            A(ind_x3, ind_x2, ind_d2) * D(ind_d2, ind_x2, ind_d1) * A(ind_x2, ind_x1, ind_d1) * D1(ind_d1, ind_x1) * A1(ind_x1);
    end
    end
    end
    end
end


%normalize
X3 = X3/sum(X3);

O3O4 = zeros(o, o);

for ind_o3 = 1:o
    for ind_o4 = 1:o
        
        for ind_x3 = 1:x
            O3O4(ind_o3, ind_o4) = O3O4(ind_o3, ind_o4) + ...
                O(ind_o3, ind_x3) * O4X3(ind_o4, ind_x3) * X3(ind_x3);
        end
        
    end
end

