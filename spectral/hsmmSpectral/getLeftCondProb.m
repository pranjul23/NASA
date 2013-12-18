%compute conditional probability OiOj..Ok|XcDc (to the left of XcDc)
%type 0: |XcDc
%type 1: |XcDc-1

function res = getLeftCondProb(O, D, X, Prior, ind, type)


assert(length(ind) >= 3);

[o x] = size(O);
d = size(D,1);

%initialization
T = tensor(D(:,:,1)*diag(Prior));
T = tenmat(T, [2 1]);
T = T.data;


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


%embed X
%x_t d_t-1 x_t-1 d_t-1
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


%first step
for k=1:ind(1)-1
    T = D*X*T;
end

curr_ind = ind(1);


%need to do it differently if type==1
%run until ind(end)-1 and manually do the last multiplication
if type == 1
    ind(end) = ind(end)-1;
end

%main loop
for i=2:length(ind)
    
    T = kr(T', kron(ones(1,d), eye(x)))';
    
    for k=curr_ind:ind(i)-1
        T = D*X*T;
    end
        
    curr_ind = ind(i);
end


%manually do the last multiplications to get |XcDc-1
if type == 1
    T = X*T;
end


%reverse the order of indeces for variables
%(after the above they will be in reverse order: X3X2X1|X4D4)
dim = [ones(1,length(ind)-1)*x x d];
res = tensor(zeros(dim));
res = tenmat(res, 1 : length(ind)-1);
res(:) = T';
res = tensor(res);
T = tenmat(res, length(ind)-1:-1:1);
T = T.data;

% 
% res=T;
% 
% [U S V] = svd(T);
% 
% min(diag(S))

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




