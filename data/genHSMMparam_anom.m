function [A0 D0 A D O] = genHSMMparam_anom(Nobs, Nhid, Dmax, Dmin, A0_true, D0_true, A_true, D_true, O_true)

%initial state distribution
A0 = rand(Nhid,1);
A0 = A0./sum(A0); %normalize

A0 = A0_true;


%initial duration distribution
D0 = rand(Dmax,1);
%D0(1:Dmin) = 0;

D0 = D0./sum(D0); %normalize



%========= transition distribution ==========
%element i,j,k is transition from j to i when d_{t-1} = k: p(i|j,k)

A1 = rand(Nhid);
for i=1:Nhid
    A1(:,i) = A1(:,i)/sum(A1(:,i));
end

A1 = A_true(:,:,1);

A = cat(3, A1, eye(Nhid));
for i=3:Dmax
    A = cat(3, A, eye(Nhid));
end



%========= duration distribution ============
%element i,j is duration of i time units, given we are in state j and d_{t-1} = k

D1 = full(sprand(Dmax, Nhid, 0.5));
%D1(1:Dmin,:) = 0;
for i=1:Nhid
    D1(:,i) = D1(:,i)/sum(D1(:,i));
end


tmp = ones(1, Nhid);
M = [tmp; zeros(Dmax-1, Nhid)];
D = cat(3, D1, M);

for i=2:Dmax-1
    M = [zeros(i-1, Nhid); tmp; zeros(Dmax-i, Nhid)];
    D = cat(3, D, M);
end


%========= observation distribution =========
%element i,j is observation of symbol i, given we are in state j and d_{t-1} = k

O = full(sprand(Nobs, Nhid, 0.5));
for i=1:Nhid
    
    if sum(O(:,i)) == 0
        O(randi(Nobs,1),i)=1;
    end
    
    O(:,i) = O(:,i)/sum(O(:,i));
end

%O = O_true;










