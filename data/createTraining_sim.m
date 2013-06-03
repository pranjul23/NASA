function [] = createTraining_sim(T, numSeq, Nobs, Nhid, Dmax, A0, D0, A, D, O, ID)


loc = strcat('../libdai/examples/data/HSMMtraining_',num2str(ID),'.txt');
fidhsmm = fopen(loc, 'w');

loc = strcat('../libdai/examples/data/HMMtraining_',num2str(ID),'.txt');
fidhmm = fopen(loc, 'w');

fprintf(fidhsmm, '%d\n',numSeq);
fprintf(fidhmm, '%d\n',numSeq);

%save data for Murphyk HMM function 
train = cell(numSeq,1);
initState = zeros(numSeq,1);
initDur = zeros(numSeq,1);


for i=1:numSeq
    
    prevState = randsample(1:Nhid, 1, true, A0);
    prevDur = randsample(1:Dmax, 1, true, D0);

    initState(i) = prevState;
    initDur(i) = prevDur;
    
    currTime = 1;
    observations = [];
    states = [];
    durations = [];
    
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
        durations = [durations, dur];
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
    train{i} = [durations; states; observations + ones(size(observations))];
    
end

save('murphykHMMtrainData.mat', 'train', 'initState', 'initDur');

fclose(fidhsmm);
fclose(fidhmm);



















