%simulator script to generate structure for 
%explicit duration HSMM

clear;
clc;

%simulation time
T = 3;

%number of observation symbols
%including NULL
Nobs = 3;

%number of hidden states
Nhid = Nobs-1;


%number of observation steps
Dmax = 3;

%initial state distribution
A0 = rand(Nhid,1);
A0 = A0./sum(A0); %normalize

% pi = [1 0 0 0 0]';

%initial duration distribution
D0 = rand(Dmax,1);
D0 = D0./sum(D0); %normalize


%========= transition distribution ==========
%element i,j,k is transition from j to i when d_{t-1} = k: p(i|j,k)

A1 = rand(Nhid);
for i=1:Nhid
    A1(:,i) = A1(:,i)/sum(A1(:,i));
end
%A1 = [[.3 .3 .4 0 0]' [0 0 0 1 0]' [0 0 .2 .8 0]' [0 .4 0 0 .6]' [0 0 1 0 0]'];

A = cat(3, A1, eye(Nhid));
for i=3:Dmax
    A = cat(3,A,eye(Nhid));
end

%========= duration distribution ============
%element i,j is duration of i time units, given we are in state j and d_{t-1} = k

D1 = rand(Dmax, Nhid);
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

O1 = eye(Nhid, Nhid);
O1 = [O1; zeros(1, Nhid)];
%O1 = rand(Nobs, Nhid);
%for i=1:Nhid
%    O1(:,i) = O1(:,i)/sum(O1(:,i));
%end

tmp = zeros(Nhid, Nhid);
M = [tmp; ones(1, Nhid)];
O = cat(3, O1, M);

for i=3:Dmax
    O = cat(3, O, M);
end

%==========================================================================
%write to factor graph file

fid = fopen('../libdai/examples/hsmm_factor_graph.fg', 'w');
fprintf(fid, '%d\n\n',3*T + 2);



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




%factor: p(O1|A1, D0)
fprintf(fid, '%d\n', 3); %number of variables
fprintf(fid, '%d\t%d\t%d\n', 4, 3, 0); %variable ids in the factor
fprintf(fid, '%d\t%d\t%d\n', Nobs, Nhid, Dmax); %dim of variables in the factor
fprintf(fid, '%d\n', nnz(O)); %number of nonzero elements

elem = reshape(O, numel(O),1);
nnzind = find(elem > 0);
nnzind = nnzind - ones(length(nnzind),1);
elem(elem==0)=[];

data = [nnzind elem];
fprintf(fid, '%d %.4g\n', data');

fprintf(fid, '\n');




%now write the rest of factors
for i=2:T
    
    %factor: p(dt|at, d_t-1)
    fprintf(fid, '%d\n', 3); %number of variables
    fprintf(fid, '%d\t%d\t%d\n', 3*i-1, 3*i, 3*i-4); %variable ids in the factor
    fprintf(fid, '%d\t%d\t%d\n', Dmax, Nhid, Dmax); %dim of variables in the factor
    fprintf(fid, '%d\n', nnz(D)); %number of nonzero elements
    
    elem = reshape(D, numel(D),1);
    nnzind = find(elem > 0);
    nnzind = nnzind - ones(length(nnzind),1);
    elem(elem==0)=[];
    
    data = [nnzind elem];
    fprintf(fid, '%d %.4g\n', data');
    
    fprintf(fid, '\n');
    
    
    
    
    %factor: p(at|a_t-1, d_t-1)
    fprintf(fid, '%d\n', 3); %number of variables
    fprintf(fid, '%d\t%d\t%d\n', 3*i, 3*i-3, 3*i-4); %variable ids in the factor
    fprintf(fid, '%d\t%d\t%d\n', Nhid, Nhid, Dmax); %dim of variables in the factor
    fprintf(fid, '%d\n', nnz(A)); %number of nonzero elements
    
    elem = reshape(A, numel(A),1);
    nnzind = find(elem > 0);
    nnzind = nnzind - ones(length(nnzind),1);
    elem(elem==0)=[];
    
    data = [nnzind elem];
    fprintf(fid, '%d %.4g\n', data');
    
    fprintf(fid, '\n');
    
    
    
    
    %factor: p(ot|at, d_t-1)
    fprintf(fid, '%d\n', 3); %number of variables
    fprintf(fid, '%d\t%d\t%d\n', 3*i+1, 3*i, 3*i-4); %variable ids in the factor
    fprintf(fid, '%d\t%d\t%d\n', Nobs, Nhid, Dmax); %dim of variables in the factor
    fprintf(fid, '%d\n', nnz(O)); %number of nonzero elements
    
    elem = reshape(O, numel(O),1);
    nnzind = find(elem > 0);
    nnzind = nnzind - ones(length(nnzind),1);
    elem(elem==0)=[];
    
    data = [nnzind elem];
    fprintf(fid, '%d %.4g\n', data');
    
    fprintf(fid, '\n');
    
end



%generate output sequence

prevState = randsample(1:Nhid, 1, true, A0);
prevDur = randsample(1:Dmax, 1, true, D0);

currTime = 0;
observations = [];
states = [];


while(currTime <= T)
    
    %select curr state
    currState = randsample(1:Nhid, 1, true, A(:, prevState, prevDur))';
    
    %select duration (# of timesteps) until next observation
    dur = randsample(1:Dmax, 1, true, D(:,currState, prevDur))';
    
    %perform dur observations in current state
    currObs = randsample(1:Nobs, 1, true, O(:,currState, prevDur))';
    
    %advance time
    currTime = currTime + 1;
    
    prevDur = dur;    
    
    %save data
    states = [states, currState];
    observations = [observations; currObs];
end

%observations - ones(size(observations))



































