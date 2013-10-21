function [A0 A Afull D O] = genHSMMparam_init(Nobs, Nhid, Dmin, Dmax)

%initial state distribution
A0 = rand(Nhid,1);
A0 = A0./sum(A0); %normalize


%========= transition distribution ==========
%element i,j,k is transition from j to i when d_{t-1} = k: p(i|j,k)

A1 = rand(Nhid);
Afull = A1;
A1 = tril(A1,-1)+triu(A1,1);

for i=1:Nhid
    A1(:,i) = A1(:,i)/sum(A1(:,i));
    Afull(:,i) = Afull(:,i)/sum(Afull(:,i));
end

A = cat(3, A1, eye(Nhid));
for i=3:Dmax
    A = cat(3, A, eye(Nhid));
end

%========= duration distribution ============
%element i,j is duration of i time units, given we are in state j and d_{t-1} = k

D1 = rand(Dmax, Nhid);
D1(1:(Dmin-1),:) = 0;

for i=1:Nhid
    D1(:,i) = D1(:,i)/sum(D1(:,i));
end
%D1 = [[.1 .1 .5 .3 0 0 0]' [0 0 0 .2 .8 0 0]' [0 0 0 0 .4 .4 .2]' [.1 .3 .3 .3 0 0 0]' [0 0 0 0 .3 .7 0]'];

tmp = ones(1, Nhid);
M = [tmp; zeros(Dmax-1, Nhid)];
D = cat(3, D1, M);

for i=2:Dmax-1
    M = [zeros(i-1, Nhid); tmp; zeros(Dmax-i, Nhid)];
    D = cat(3, D, M);
end

%========= observation distribution =========
%element i,j is observation of symbol i, given we are in state j and d_{t-1} = k

O = rand(Nobs, Nhid);
for i=1:Nhid
    O(:,i) = O(:,i)/sum(O(:,i));
end


