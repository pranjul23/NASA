function [] = createTraining(T, numSeq, Nobs, Nhid, Dmax, A0, D0, A, D, O)


train = cell(numSeq,1);

fidhsmm = fopen('../libdai/examples/HSMMtraining.txt', 'w');
fidhmm = fopen('../libdai/examples/HMMtraining.txt', 'w');

fprintf(fidhsmm, '%d\n',numSeq);
fprintf(fidhmm, '%d\n',numSeq);


for i=1:numSeq
    
    prevState = randsample(1:Nhid, 1, true, A0);
    prevDur = randsample(1:Dmax, 1, true, D0);

    currTime = 1;
    observations = [];
    states = [];
    
    while(currTime <= T)
        
        %select curr state
        currState = randsample(1:Nhid, 1, true, A(:, prevState, prevDur))';
        
        %select duration (# of timesteps) until next observation
        dur = randsample(1:Dmax, 1, true, D(:,currState, prevDur))';
        
        %perform dur observations in current state
        currObs = randsample(1:Nobs, 1, true, O(:,currState))';
        
        %advance time
        currTime = currTime + 1;        
        prevDur = dur;
        prevState = currState;
        
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
    train{i} = observations + ones(size(observations));    
end

save('murphykHMMtrainData.mat', 'train');

fclose(fidhsmm);
fclose(fidhmm);



















