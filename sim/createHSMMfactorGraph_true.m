function [] = createHSMMfactorGraph_true(Nobs, Nhid, Dmax, A0, D0, A, D, O)


fid = fopen('../libdai/examples/hsmm_factor_graph_init_true.fg', 'w');
fprintf(fid, '%d\n\n',5);


%factor: p(D0)
fprintf(fid, '%d\n', 1); %number of variables
fprintf(fid, '%d\n', 0); %variable ids in the factor
fprintf(fid, '%d\n', Dmax); %dim of variables in the factor
fprintf(fid, '%d\n', nnz(D0)); %number of nonzero elements

elem = D0;
nnzind = find(elem > 0);
nnzind = nnzind - ones(length(nnzind),1);
elem(elem==0)=[];

data = [nnzind elem];
fprintf(fid, '%d %.4g\n', data');

fprintf(fid, '\n');


%factor: p(A0)
fprintf(fid, '%d\n', 1); %number of variables
fprintf(fid, '%d\n', 1); %variable ids in the factor
fprintf(fid, '%d\n', Nhid); %dim of variables in the factor
fprintf(fid, '%d\n', nnz(A0)); %number of nonzero elements

elem = A0;
nnzind = find(elem > 0);
nnzind = nnzind - ones(length(nnzind),1);
elem(elem==0)=[];

data = [nnzind elem];
fprintf(fid, '%d %.4g\n', data');

fprintf(fid, '\n');



%factor: p(D1|A1, D0)
fprintf(fid, '%d\n', 3); %number of variables
fprintf(fid, '%d\t%d\t%d\n', 2, 3, 0); %variable ids in the factor
fprintf(fid, '%d\t%d\t%d\n', Dmax, Nhid, Dmax); %dim of variables in the factor
fprintf(fid, '%d\n', nnz(D)); %number of nonzero elements

elem = reshape(D, numel(D),1);
nnzind = find(elem > 0);
nnzind = nnzind - ones(length(nnzind),1);
elem(elem==0)=[];

data = [nnzind elem];
fprintf(fid, '%d %.4g\n', data');

fprintf(fid, '\n');




%factor: p(A1|A0, D0)
fprintf(fid, '%d\n', 3); %number of variables
fprintf(fid, '%d\t%d\t%d\n', 3, 1, 0); %variable ids in the factor
fprintf(fid, '%d\t%d\t%d\n', Nhid, Nhid, Dmax); %dim of variables in the factor
fprintf(fid, '%d\n', nnz(A)); %number of nonzero elements

elem = reshape(A, numel(A),1);
nnzind = find(elem > 0);
nnzind = nnzind - ones(length(nnzind),1);
elem(elem==0)=[];

data = [nnzind elem];
fprintf(fid, '%d %.4g\n', data');

fprintf(fid, '\n');




%factor: p(O1|A1)
fprintf(fid, '%d\n', 2); %number of variables
fprintf(fid, '%d\t%d\n', 4, 3); %variable ids in the factor
fprintf(fid, '%d\t%d\n', Nobs, Nhid); %dim of variables in the factor
fprintf(fid, '%d\n', nnz(O)); %number of nonzero elements

elem = reshape(O, numel(O),1);
nnzind = find(elem > 0);
nnzind = nnzind - ones(length(nnzind),1);
elem(elem==0)=[];

data = [nnzind elem];
fprintf(fid, '%d %.4g\n', data');

fprintf(fid, '\n');

fclose(fid);




