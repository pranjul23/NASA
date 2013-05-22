function [A0 A O] = genHMMparam(Nobs, Nhid)

%initial state distribution
A0 = rand(Nhid,1);
A0 = A0./sum(A0); %normalize


%========= transition distribution ==========
%element i,j is transition from j to i: p(i|j)

A = rand(Nhid);
for i=1:Nhid
    A(:,i) = A(:,i)/sum(A(:,i));
end

%========= observation distribution =========
%element i,j is observation of symbol i, given we are in state j 

O = rand(Nobs, Nhid);
for i=1:Nhid
    O(:,i) = O(:,i)/sum(O(:,i));
end
