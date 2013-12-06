function res = embed2(A)


x = size(A,1);
d = size(A,4);


%x_t x_t-1 x_t-2 x_t-2 d_t-2
res = zeros(x,x,x,x,d);

for ind1 = 1:x
for ind2 = 1:x
for ind3 = 1:x
for ind4 = 1:d    
    
    res(ind3, ind1, ind2, ind3, ind4) = A(ind1,ind2,ind3,ind4);
    
end
end
end
end

