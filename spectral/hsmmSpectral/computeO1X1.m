function res = computeO1X1(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);


res = zeros(o, x);

for ind_o1 = 1:o
for ind_x1 = 1:x   
   
    
    res(ind_o1, ind_x1) = res(ind_o1, ind_x1) + ...
        O(ind_o1, ind_x1) * A1(ind_x1);
        
    
end
end



