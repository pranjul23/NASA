function result = testSpectModel_parallel(test, obsDim, numObs, ...
                                          rootTensor,  ...
                                          tailTensor,  ...
                                          obsTensor,  ...
                                          tranTensor, ...
                                          durTensor)
                                   
                                   
                                   
set(findResource(), 'ClusterSize', 8);

L = size(test,1);
result = zeros(L, 1);

%number of cores
C = 6;

assert(L >= C);

%figure out how many jobs per core
extra = rem(L, C);

%spread evenly
numJobs = ((L-extra)/C)*ones(C,1);

%place extra jobs
ind = randsample(1:C, extra, 'false');
numJobs(ind) = numJobs(ind) + 1;

c = 1;

%start all the jobs
for i=1:C
    
    if iscell(test)
        sequence = test(c:c+numJobs(i)-1);
    else
        sequence = test(c:c+numJobs(i)-1,:);
    end

    eval(['job' num2str(i) ' = batch(@testSpectModel, 1, {sequence, obsDim, numObs,'...
                                    ' rootTensor, tailTensor,'...
                                    ' obsTensor, tranTensor, durTensor});']);
    c = c + numJobs(i);
end


%check when all jobs are finished 
for i = 1:C    
    eval(['wait(job' num2str(i) ', ''finished'')']);
end

c = 1;
for i = 1:C    
    eval(['tmp = job' num2str(i) '.getAllOutputArguments;']);    
    result(c:c+numJobs(i)-1) = tmp{1};
    
    c = c + numJobs(i);
end

for i = 1:C    
    eval(['delete(job' num2str(i) ');']);
end





