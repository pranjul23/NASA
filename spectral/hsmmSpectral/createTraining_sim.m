function train = createTraining_sim(T, numSeq, Nobs, Nhid, Dmax, A0, A, D, O, ID)

loc = strcat('../../libdai/examples/data/HSMMtraining_',num2str(ID),'.txt');
fidhsmm = fopen(loc, 'w');

fprintf(fidhsmm, '%d\n',numSeq);

train = zeros(numSeq, T);

for i=1:numSeq
    
    prevState = randsample(1:Nhid, 1, true, A0);
    prevDur = 1;
    
    currTime = 1;
    observations = zeros(1, T);    
        
    %initialization
    prevDur = randsample(1:Dmax, 1, true, D(:,prevState, prevDur))';    
    Obs = randsample(1:Nobs, 1, true, O(:,prevState))';    
    observations(1) = Obs;    
    
    currTime = currTime + 1;
            
    while(currTime <= T)
        
        %select curr state
        currState = randsample(1:Nhid, 1, true, A(:, prevState, prevDur))';
        
        %select duration (# of timesteps) until next observation
        dur = randsample(1:Dmax, 1, true, D(:,currState, prevDur))';
        
        %perform dur observations in current state
        currObs = randsample(1:Nobs, 1, true, O(:,currState))';
        
        %save data        
        observations(currTime) = currObs;
        
        %advance time
        currTime = currTime + 1;        
        prevDur = dur;
        prevState = currState;                
    end
    
    lenObs = length(observations);
    fprintf(fidhsmm, '%d\n', lenObs);
    
    dataHSMM = [[2:3:3*(lenObs-1)-1 3*(lenObs-1)+1]; observations-1];
    
    format = repmat('%d\t', 1, lenObs-1);
    
    fprintf( fidhsmm, [format,'%d\n'], dataHSMM' );
    
    train(i,:) = observations;
end

fclose(fidhsmm);


















