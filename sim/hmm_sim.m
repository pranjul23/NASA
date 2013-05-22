%simulator script to generate structure for HMM

clear;
clc;


load matlabHMMdata.mat;

%number of observation symbols
Nobs = 5;


%% ================ INIT PARAMETERS =======================================

%number of hidden states
Nhid_init = 3;

%initial state distribution
A0 = rand(Nhid_init,1);
A0_init = A0./sum(A0); %normalize


%========= transition distribution ==========
%element i,j is transition from j to i: p(i|j)

A_init = rand(Nhid_init);
for i=1:Nhid_init
    A_init(:,i) = A_init(:,i)/sum(A_init(:,i));
end

%========= observation distribution =========
%element i,j is observation of symbol i, given we are in state j 

O_init = rand(Nobs, Nhid_init);
for i=1:Nhid_init
    O_init(:,i) = O_init(:,i)/sum(O_init(:,i));
end


%==========================================================================
%write to factor graph file

fid = fopen('../libdai/examples/hmm_factor_graph_init.fg', 'w');
fprintf(fid, '%d\n\n',3);



%factor: p(A0)
fprintf(fid, '%d\n', 1); %number of variables
fprintf(fid, '%d\n', 1); %variable ids in the factor
fprintf(fid, '%d\n', Nhid_init); %dim of variables in the factor
fprintf(fid, '%d\n', nnz(A0_init)); %number of nonzero elements

elem = A0_init;
nnzind = find(elem > 0);
nnzind = nnzind - ones(length(nnzind),1);
elem(elem==0)=[];

data = [nnzind elem];
fprintf(fid, '%d %.4g\n', data');

fprintf(fid, '\n');



%factor: p(A1|A0)
fprintf(fid, '%d\n', 2); %number of variables
fprintf(fid, '%d\t%d\n', 1, 0); %variable ids in the factor
fprintf(fid, '%d\t%d\n', Nhid_init, Nhid_init); %dim of variables in the factor
fprintf(fid, '%d\n', nnz(A_init)); %number of nonzero elements

elem = reshape(A_init, numel(A_init),1);
nnzind = find(elem > 0);
nnzind = nnzind - ones(length(nnzind),1);
elem(elem==0)=[];

data = [nnzind elem];
fprintf(fid, '%d %.4g\n', data');

fprintf(fid, '\n');




%factor: p(O1|A1)
fprintf(fid, '%d\n', 2); %number of variables
fprintf(fid, '%d\t%d\n', 2, 1); %variable ids in the factor
fprintf(fid, '%d\t%d\n', Nobs, Nhid_init); %dim of variables in the factor
fprintf(fid, '%d\n', nnz(O_init)); %number of nonzero elements

elem = reshape(O_init, numel(O_init),1);
nnzind = find(elem > 0);
nnzind = nnzind - ones(length(nnzind),1);
elem(elem==0)=[];

data = [nnzind elem];
fprintf(fid, '%d %.4g\n', data');

fprintf(fid, '\n');

fclose(fid);


