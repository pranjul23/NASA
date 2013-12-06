function res = embedD(D)


d = size(D,1);
x = size(D,2);


%%% x_t d_t x_t d_t-1

%d_t x_t x_t d_t-1
res = zeros(d,x,x,d);

for ind1 = 1:d
for ind2 = 1:x
for ind3 = 1:d
    
    res(ind1,ind2,ind2,ind3) = D(ind1,ind2,ind3);
    
end
end
end

