function res = embed(A)


x = size(A,1);
d = size(A,3);

%%%x_t x_t x_t-1 d_t-1 d_t-1

%x_t-1 x_t x_t-1 d_t-1
res = zeros(x,d,x,d);

for ind1 = 1:x
for ind2 = 1:x
for ind3 = 1:d
    
    res(ind1, ind3, ind2, ind3) = A(ind1,ind2,ind3);
    
end
end
end

