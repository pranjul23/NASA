%compute conditional probability OiOj..Ok|XcDc
%type 0: |XcDc
%type 1: |XcDc-1

function res = getCondProb(O, D, X, ind, type)


assert(length(ind) >= 3);

[o x] = size(O);
d = size(D,1);


%embed D
% x_t d_t x_t d_t-1
res = zeros(x,d,x,d);

for ind1 = 1:d
for ind2 = 1:x
for ind3 = 1:d
    
    res(ind2,ind1,ind2,ind3) = D(ind1,ind2,ind3);
    
end
end
end

D = tensor(res);
D = tenmat(D,[1 2]);
D = D.data;


%create start tensor
T = tensor(X);
T = tenmat(X, 1);
T = T.data;

%embed X
%x_t-1 x_t x_t-1 d_t-1
res = zeros(x,d,x,d);

for ind1 = 1:x
for ind2 = 1:x
for ind3 = 1:d
    
    res(ind1, ind3, ind2, ind3) = X(ind1,ind2,ind3);
    
end
end
end

X = tensor(res);
X = tenmat(X,[1 2]);
X = X.data;


%initialization
l = ind(end)-ind(end-1)-1;
for j=1:l
    T = T*D*X;
end


%main loop
for i=length(ind)-1:-1:2
    
    T = kr(T, kron(ones(1,d), eye(x)));
    
    l = ind(i)-ind(i-1);    
    for j=1:l
        T = T*D*X;
    end            
end


%extra term if we compute |XcDc-1
if type == 1
    T = T*D;
end


%prepare observation matrix
K = O;
for i=1:length(ind)-2
    K = kron(K, O);
end

Z = K*T;

%transform result into tensor
dim = [ones(1,length(ind)-1)*o x d];

res = tensor(zeros(dim));
res = tenmat(res, 1 : length(ind)-1);
res(:) = Z;
res = tensor(res);



