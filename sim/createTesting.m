function [] = createTesting(T, numSeq, Nobs, ...
                            Nhid_true, Dmax_true, Nhid_anom, Dmax_anom, ...
                            A0_true, D0_true, A_true, D_true, O_true, ...
                            A0_anom, D0_anom, A_anom, D_anom, O_anom)


fidhsmm = fopen('../libdai/examples/HSMMtesting.txt', 'w');
fidhmm = fopen('../libdai/examples/HMMtesting.txt', 'w');

fprintf(fidhsmm, '%d\n',numSeq);
fprintf(fidhmm, '%d\n',numSeq);

%data for Murphyk HMM function
test = cell(numSeq,1);


%generate normal sequences
for i=1:floor(numSeq/2)
    
    prevState = randsample(1:Nhid_true, 1, true, A0_true);
    prevDur = randsample(1:Dmax_true, 1, true, D0_true);

    currTime = 1;
    observations = [];
    states = [];
    
    while(currTime <= T)
        
        %select curr state
        currState = randsample(1:Nhid_true, 1, true, A_true(:, prevState, prevDur))';
        
        %select duration (# of timesteps) until next observation
        dur = randsample(1:Dmax_true, 1, true, D_true(:,currState, prevDur))';
        
        %perform observations in current state
        currObs = randsample(1:Nobs, 1, true, O_true(:,currState))';
        
        %advance time
        currTime = currTime + 1;
        
        prevDur = dur;
        
        %save data
        states = [states, currState];
        observations = [observations, currObs-1];
    end
    
    lenObs = length(observations);
    fprintf(fidhsmm, '%d\n', lenObs);
    fprintf(fidhmm, '%d\n', lenObs);
    
    dataHSMM = [4:3:3*lenObs+1; observations];
    dataHMM = [2:2:2*lenObs; observations];
    
    format = repmat('%d\t', 1, lenObs-1);
    fprintf(fidhsmm,[format,'%d\n'], dataHSMM');
    fprintf(fidhmm,[format,'%d\n'], dataHMM');
    
    %save data for Murphyk HMM function
    test{i} = observations + ones(size(observations));  
end



%generate anomalous sequences
for i=1:ceil(numSeq/2)
    
    prevState = randsample(1:Nhid_anom, 1, true, A0_anom);
    prevDur = randsample(1:Dmax_anom, 1, true, D0_anom);

    currTime = 1;
    observations = [];
    states = [];
    
    while(currTime <= T)
        
        %select curr state
        currState = randsample(1:Nhid_anom, 1, true, A_anom(:, prevState, prevDur))';
        
        %select duration (# of timesteps) until next observation
        dur = randsample(1:Dmax_anom, 1, true, D_anom(:,currState, prevDur))';
        
        %perform dur observations in current state
        currObs = randsample(1:Nobs, 1, true, O_anom(:,currState))';
        
        %advance time
        currTime = currTime + 1;
        
        prevDur = dur;
        
        %save data
        states = [states, currState];
        observations = [observations, currObs-1];
    end
    
    lenObs = length(observations);
    fprintf(fidhsmm, '%d\n', lenObs);
    fprintf(fidhmm, '%d\n', lenObs);
    
    dataHSMM = [4:3:3*lenObs+1; observations];
    dataHMM = [2:2:2*lenObs; observations];
        
    format = repmat('%d\t', 1, lenObs-1);
    
    fprintf(fidhsmm, [format,'%d\n'], dataHSMM');    
    fprintf(fidhmm, [format,'%d\n'], dataHMM');
    
    %save data for Murphyk HMM function
    test{i + floor(numSeq/2)} = observations + ones(size(observations));  
end

save('murphykHMMtestData.mat', 'test');

fclose(fidhsmm);
fclose(fidhmm);












