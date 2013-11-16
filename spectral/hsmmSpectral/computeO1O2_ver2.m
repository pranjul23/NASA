%compute p(O1,O2) for HSMM

function O1O2 = computeO1O2_ver2(A1, A, D1, O)


[o x] = size(O);
d = size(D1,1);


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



%copmute p(O1|x2)

O1X2 = zeros(o,x);

for ind_o1 = 1:o
    for ind_x2 = 1:x
        
        for ind_d1 = 1:d
            for ind_x1 = 1:x
                
                O1X2(ind_o1, ind_x2) = O1X2(ind_o1, ind_x2) + ...
                    (O(ind_o1, ind_x1) * A(ind_x2, ind_x1, ind_d1) * D1(ind_d1, ind_x1) * A1(ind_x1))/X2(ind_x2);
                
            end
        end
        
    end
end

%normalize
for ind_x2 = 1:x
    O1X2(:, ind_x2) = O1X2(:, ind_x2)/sum(O1X2(:, ind_x2));
end


O1O2 = zeros(o, o);

for ind_o1 = 1:o
    for ind_o2 = 1:o
        
        for ind_x2 = 1:x
            O1O2(ind_o1, ind_o2) = O1O2(ind_o1, ind_o2) + O(ind_o2, ind_x2) * O1X2(ind_o1, ind_x2) * X2(ind_x2);
        end
        
    end
end

