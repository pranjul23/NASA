function res = embedA1(A,D1,A1)


x = size(A,1);
d = size(A,3);

%x_t x_t x_t-1 d_t-1
res = zeros(x,x,x,d);

for ind1 = 1:x
for ind2 = 1:x
for ind3 = 1:d
    
    res(ind1, ind1, ind2, ind3) = A(ind1,ind2,ind3)*D1(ind3,ind2)*A1(ind2);
    
end
end
end

