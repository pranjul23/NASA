%compute p(O1,O2) for HSMM

function O1O2 = computeO1O2(A1, A, D1, D, O)


[o x] = size(O);
d = size(D1,1);
%copmute p(O2|x1)

O2X1 = zeros(o,x);

for ind_o2 = 1:o
    for ind_x1 = 1:x
        
        for ind_d1 = 1:d
            for ind_x2 = 1:x
                
                O2X1(ind_o2, ind_x1) = O2X1(ind_o2, ind_x1) + D1(ind_d1, ind_x1) * A(ind_x2, ind_x1, ind_d1) * O(ind_o2, ind_x2);
                
            end
        end
        
    end
end

O1O2 = zeros(o, o);

for ind_o1 = 1:o
    for ind_o2 = 1:o
        
        for ind_x1 = 1:x
            O1O2(ind_o1, ind_o2) = O1O2(ind_o1, ind_o2) + O(ind_o1, ind_x1) * O2X1(ind_o2, ind_x1) * A1(ind_x1);
        end
        
    end
end

