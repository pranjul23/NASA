%simulator script to generate structure for
%explicit duration HSMM

clear;
clc;

%number of observation symbols
%including NULL
Nobs = 183;

%number of hidden states
Nhid = 30;

%max duration
Dmax = 80;

%initial state distribution
A0 = rand(Nhid,1);
A0 = A0./sum(A0); %normalize

%initial duration distribution
D0 = rand(Dmax,1);
D0 = D0./sum(D0); %normalize


%========= transition distribution ==========
%element i,j,k is transition from j to i when d_{t-1} = k: p(i|j,k)

A1 = rand(Nhid);
for i=1:Nhid
    A1(:,i) = A1(:,i)/sum(A1(:,i));
end

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

tmp = ones(1, Nhid);
M = [tmp; zeros(Dmax-1, Nhid)];
D = cat(3, D1, M);

for i=2:Dmax-1
    M = [zeros(i-1, Nhid); tmp; zeros(Dmax-i, Nhid)];
    D = cat(3, D, M);
end

%========= observation distribution =========
%element i,j is observation of symbol i, given we are in state j

O = rand(Nobs, Nhid);
for i=1:Nhid
   O(:,i) = O(:,i)/sum(O(:,i));
end


%==========================================================================
%write to factor graph file

fid = fopen('../libdai/examples/hsmm_factor_graph_init.fg', 'w');
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



%==========================================================================
%write training sequence to file

%directory to read training data from
dir_train = '/Users/igor/Documents/Projects/anomaly/code/HSMM/sim/lpr-mit-live-normal/';

list_train = dir(dir_train);
numSeq_train = length(list_train)-3;

%the result will be written to 
fid_train = fopen('/Users/igor/Documents/Projects/anomaly/code/HSMM/libdai/examples/training.txt', 'w');

names = importdata('/Users/igor/Documents/Projects/anomaly/code/HSMM/sim/lpr-mit-live-callnames.txt');
out = fopen('/Users/igor/Documents/Projects/anomaly/code/HSMM/sim/commands.txt', 'w');


%find indeces of sequences of small lengths
%ind = 1:numSeq_train;

ind = [];

for i=1:numSeq_train
    obs = load(strcat(dir_train, list_train(i+3).name));
    obs = obs(:,2);
    lenObs = length(obs);           
    
    if lenObs <= 160
        ind = [ind i];
    end
end

fprintf(fid_train, '%d\n',length(ind));

for i=1:length(ind)
    
    obs = load(strcat(dir_train, list_train(ind(i)+3).name));
    obs = obs(:,2);
    
    %plot(obs, 'r*')
    
    %sequence = names(obs+ones(size(obs)));
    %fprintf(out, '%s\n', sequence{:});
    
    lenObs = length(obs);            
    fprintf(fid_train, '%d\n', lenObs);
    
    data = [4:3:3*lenObs+1; obs'];
    
    format = repmat('%d\t', 1, lenObs-1);
    fprintf(fid_train, [format,'%d\n'], data');    
end

fclose(fid_train);



%==========================================================================
%write testing sequence to file

%directory to read anomaly data from
dir_test = '/Users/igor/Documents/Projects/anomaly/code/HSMM/sim/lpr-mit-live-anomaly/';

list_test = dir(dir_test);
numSeq_test = length(list_test)-3;

%the result will be written to 
fid_test = fopen('/Users/igor/Documents/Projects/anomaly/code/HSMM/libdai/examples/testing.txt', 'w');

%names = importdata('/Users/igor/Documents/Projects/anomaly/code/HSMM/sim/lpr-mit-live-callnames.txt');
%out = fopen('/Users/igor/Documents/Projects/anomaly/code/HSMM/sim/commands.txt', 'w');


fprintf(fid_test, '%d\n', 20);

%first 10 sequences are normal ones
for i=1:10
    
    obs = load(strcat(dir_train, list_train(i+3).name));
    obs = obs(:,2);
    
    %plot(obs, 'r*')
    
    %sequence = names(obs+ones(size(obs)));
    %fprintf(out, '%s\n', sequence{:});
    
    lenObs = length(obs);            
    fprintf(fid_test, '%d\n', lenObs);
    
    data = [4:3:3*lenObs+1; obs'];
    
    format = repmat('%d\t', 1, lenObs-1);
    fprintf(fid_test, [format,'%d\n'], data');    
end


%second 10 sequences are anomalous ones
%take 5 from top
for i=1:5
    
    obs = load(strcat(dir_test, list_test(i+3).name));
    obs = obs(:,2);
    
    %plot(obs, 'r*')
    
    %sequence = names(obs+ones(size(obs)));
    %fprintf(out, '%s\n', sequence{:});
    
    lenObs = length(obs);            
    fprintf(fid_test, '%d\n', lenObs);
    
    data = [4:3:3*lenObs+1; obs'];
    
    format = repmat('%d\t', 1, lenObs-1);
    fprintf(fid_test, [format,'%d\n'], data');    
end

%take 5 from bottom
for i = numSeq_test:-1:numSeq_test-5
    
    obs = load(strcat(dir_test, list_test(i+3).name));
    obs = obs(:,2);
    
    %plot(obs, 'r*')
    
    %sequence = names(obs+ones(size(obs)));
    %fprintf(out, '%s\n', sequence{:});
    
    lenObs = length(obs);            
    fprintf(fid_test, '%d\n', lenObs);
    
    data = [4:3:3*lenObs+1; obs'];
    
    format = repmat('%d\t', 1, lenObs-1);
    fprintf(fid_test, [format,'%d\n'], data');    
end


fclose(fid_test);













