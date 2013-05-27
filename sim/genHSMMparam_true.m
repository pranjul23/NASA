function [A0 D0 A D O] = genHSMMparam_true(Nobs, Nhid, Dmax)

%initial state distribution
A0 = full(sprand(Nhid, 1, 0.6));
if sum(A0) == 0
    A0(randi(Nhid,1))=1;
end

A0 = A0./sum(A0); %normalize

%initial duration distribution
D0 = full(sprand(Dmax, 1, 0.5));
if sum(D0) == 0
    D0(randi(Dmax,1))=1;
end

D0 = D0./sum(D0); %normalize


%========= transition distribution ==========
%element i,j,k is transition from j to i when d_{t-1} = k: p(i|j,k)

A1 = full(sprand(Nhid, Nhid, 0.6));
A1 = tril(A1,-1)+triu(A1,1);

for i=1:Nhid
    
    if sum(A1(:,i)) == 0
        ind = 1:Nhid;
        ind(ind==i)=[];
        A1(randi(ind,1),i)=1;
    end

    A1(:,i) = A1(:,i)/sum(A1(:,i));
end

A = cat(3, A1, eye(Nhid));
for i=3:Dmax
    A = cat(3, A, eye(Nhid));
end

%========= duration distribution ============
%element i,j is duration of i time units, given we are in state j and d_{t-1} = k

D1 = full(sprand(Dmax, Nhid, 0.5));
for i=1:Nhid
    
    if sum(D1(:,i)) == 0
        D1(randi(Dmax,1),i)=1;
    end
    
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



%============================================
%write results to file
filename = '../libdai/examples/true_parameters.txt';

fid = fopen(filename, 'w');
fprintf(fid, '%s\n', 'A0 = ');

for i=1:Nhid
    fprintf(fid,'%.4f\n',A0(i));
end
fprintf(fid, '%s\n', '=============================================== ');


fprintf(fid, '%s\n', 'D = ');
for i=1:Dmax
    fprintf(fid,'%.4f\n',D0(i));
end
fprintf(fid, '%s\n', '=============================================== ');


fprintf(fid, '%s\n', 'A = ');
for i=1:Nhid
    format = repmat('%.4f\t', 1, Nhid-1);
    fprintf(fid,[format,'%.4f\n'], A(i,:,1));
end
fprintf(fid, '%s\n', '=============================================== ');



fprintf(fid, '%s\n', 'D = ');
for i=1:Dmax
    format = repmat('%.4f\t', 1, Nhid-1);
    fprintf(fid,[format,'%.4f\n'], D(i,:,1));
end
fprintf(fid, '%s\n', '=============================================== ');



fprintf(fid, '%s\n', 'O = ');
for i=1:Nobs
    format = repmat('%.4f\t', 1, Nhid-1);
    fprintf(fid,[format,'%.4f\n'], O(i,:));
end
fclose(fid);





















