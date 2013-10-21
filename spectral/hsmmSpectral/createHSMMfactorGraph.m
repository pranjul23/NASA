function [] = createHSMMfactorGraph(Nobs, Nhid, Dmax, A1, A, D, O, ID, type)


%type of model paramters: init or true
loc = strcat('../../libdai/examples/data/hsmm_factor_graph_', type, '_',  num2str(ID), '.fg');
fid = fopen(loc, 'w');


fprintf(fid, '%d\n\n', 5); %total number of factors in this graph



%factor: p(A1)
fprintf(fid, '%d\n', 1); %number of variables
fprintf(fid, '%d\n', 1); %variable ids in the factor
fprintf(fid, '%d\n', Nhid); %dim of variables in the factor
fprintf(fid, '%d\n', nnz(A1)); %number of nonzero elements
fprintf(fid, '%d\n', 0); %number of nonzero elements that are possible
fprintf(fid, '%d\n', 0); %type of factor: 0 - regular, dense; 1 - sparse

elem = A1;
nnzind = find(elem > 0);
nnzind = nnzind - ones(length(nnzind),1);
elem(elem==0)=[];

data = [nnzind elem];
fprintf(fid, '%d %.4g\n', data');

fprintf(fid, '\n');


%factor: p(D1|A1)
fprintf(fid, '%d\n', 2); %number of variables
fprintf(fid, '%d\t%d\n', 0, 1); %variable ids in the factor
fprintf(fid, '%d\t%d\n', Dmax, Nhid); %dim of variables in the factor
fprintf(fid, '%d\n', nnz(D(:,:,1))); %number of nonzero elements
fprintf(fid, '%d\n', 0); %number of nonzero elements that are possible
fprintf(fid, '%d\n', 0); %type of factor: 0 - regular, dense; 1 - sparse

elem = reshape(D(:,:,1), numel(D(:,:,1)),1);
nnzind = find(elem > 0);
nnzind = nnzind - ones(length(nnzind),1);
elem(elem==0)=[];

data = [nnzind elem];
fprintf(fid, '%d %.4g\n', data');

fprintf(fid, '\n');


%factor: p(D2|A2, D1)
fprintf(fid, '%d\n', 3); %number of variables
fprintf(fid, '%d\t%d\t%d\n', 3, 4, 0); %variable ids in the factor
fprintf(fid, '%d\t%d\t%d\n', Dmax, Nhid, Dmax); %dim of variables in the factor
fprintf(fid, '%d\n', nnz(D)); %number of nonzero elements
fprintf(fid, '%d\n', Dmax*Nhid + Nhid*(Dmax-1)); %number of nonzero elements that are possible
fprintf(fid, '%d\n', 1); %type of factor: 0 - regular, dense; 1 - sparse

elem = reshape(D, numel(D),1);
nnzind = find(elem > 0);
nnzind = nnzind - ones(length(nnzind),1);
elem(elem==0)=[];

data = [nnzind elem];
fprintf(fid, '%d %.4g\n', data');

fprintf(fid, '%d\n', [ 0:Dmax*Nhid-1   (Dmax*Nhid-1) + find(D(:,:,2:end))' ] ); %indeces of nonzero elements

fprintf(fid, '\n');



%factor: p(A2|A1, D1)
fprintf(fid, '%d\n', 3); %number of variables
fprintf(fid, '%d\t%d\t%d\n', 4, 1, 0); %variable ids in the factor
fprintf(fid, '%d\t%d\t%d\n', Nhid, Nhid, Dmax); %dim of variables in the factor
fprintf(fid, '%d\n', nnz(A)); %number of nonzero elements
fprintf(fid, '%d\n', Nhid*Nhid + Nhid*(Dmax-1)); %number of nonzero elements that are possible
fprintf(fid, '%d\n', 1); %type of factor: 0 - regular, dense; 1 - sparse

elem = reshape(A, numel(A),1);
nnzind = find(elem > 0);
nnzind = nnzind - ones(length(nnzind),1);
elem(elem==0)=[];

data = [nnzind elem];
fprintf(fid, '%d %.4g\n', data');

fprintf(fid, '%d\n', [ 0:Nhid*Nhid-1   (Nhid*Nhid-1)+find(A(:,:,2:end))' ] ); %indeces of nonzero elements

fprintf(fid, '\n');



%factor: p(O1|A1)
fprintf(fid, '%d\n', 2); %number of variables
fprintf(fid, '%d\t%d\n', 2, 1); %variable ids in the factor
fprintf(fid, '%d\t%d\n', Nobs, Nhid); %dim of variables in the factor
fprintf(fid, '%d\n', nnz(O)); %number of nonzero elements
fprintf(fid, '%d\n', 0); %number of nonzero elements that are possible
fprintf(fid, '%d\n', 0); %type of factor: 0 - regular, dense; 1 - sparse

elem = reshape(O, numel(O),1);
nnzind = find(elem > 0);
nnzind = nnzind - ones(length(nnzind),1);
elem(elem==0)=[];

data = [nnzind elem];
fprintf(fid, '%d %.4g\n', data');

fprintf(fid, '\n');

fclose(fid);
