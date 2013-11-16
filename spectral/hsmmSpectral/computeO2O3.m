%compute p(O2,O3) for HSMM

function O2O3 = computeO2O3(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);

O3X2 = zeros(o,x);

for ind_o3 = 1:o
    for ind_x2 = 1:x
        
        for ind_x1 = 1:x
        for ind_x3 = 1:x
        for ind_d1 = 1:d    
        for ind_d2 = 1:d
            
                O3X2(ind_o3, ind_x2) = O3X2(ind_o3, ind_x2) + ...
                    O(ind_o3, ind_x3) * A(ind_x3, ind_x2, ind_d2) * D(ind_d2, ind_x2, ind_d1) * D1(ind_d1, ind_x1) *  A1(ind_x1);
        end
        end        
        end
        end
        
    end
end

%normalize
for ind_x2 = 1:x
    O3X2(:, ind_x2) = O3X2(:, ind_x2)/sum(O3X2(:, ind_x2));
end


X2 = zeros(x,1);

for ind_x2 = 1:x
    for ind_x1 = 1:x
    for ind_d1 = 1:d
        
        X2(ind_x2) = X2(ind_x2) + ...
            A(ind_x2, ind_x1, ind_d1) * D1(ind_d1, ind_x1) * A1(ind_x1);
        
    end
    end
end

%normalize
X2 = X2/sum(X2);


O2O3 = zeros(o, o);

for ind_o2 = 1:o
    for ind_o3 = 1:o
        
        for ind_x2 = 1:x
            O2O3(ind_o2, ind_o3) = O2O3(ind_o2, ind_o3) + ...
                O(ind_o2, ind_x2) * O3X2(ind_o3, ind_x2) * X2(ind_x2);
        end
        
    end
end

