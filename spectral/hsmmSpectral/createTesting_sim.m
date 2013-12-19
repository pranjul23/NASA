function test = createTesting_sim(T, numSeq, Nobs, ...
                                  Nhid_true, Dmax_true, Nhid_anom, Dmax_anom, ...
                                  A0_true, A_true, D_true, O_true, ...
                                  A0_anom, A_anom, D_anom, O_anom, ID)


loc = strcat('../../libdai/examples/data/HSMMtesting_',num2str(ID),'.txt');
fidhsmm = fopen(loc, 'w');

fprintf(fidhsmm, '%d\n',numSeq);

N = floor(numSeq/2);
N = numSeq; % <-- only normal sequences

testNorm = zeros(N, T);

%generate normal sequences
for i=1:N
    
    prevState = randsample(1:Nhid_true, 1, true, A0_true);
    prevDur = 1;
    
    currTime = 1;
    durations = [];
    observationsNorm = [];
    states = [];
    
    
    %initialization
    prevDur = randsample(1:Dmax_true, 1, true, D_true(:,prevState, prevDur))';    
    Obs = randsample(1:Nobs, 1, true, O_true(:,prevState))';    
    observationsNorm = [observationsNorm, Obs];    
    
    currTime = currTime + 1;        
    
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
        prevState = currState;
        
        %save data
        durations = [durations, dur];
        states = [states, currState];
        observationsNorm = [observationsNorm, currObs];
    end
    
    lenObs = length(observationsNorm);
    fprintf(fidhsmm, '%d\n', lenObs);
    
    dataHSMM = [[2:3:3*(lenObs-1)-1 3*(lenObs-1)+1]; observationsNorm-1];
    
    format = repmat('%d\t', 1, lenObs-1);
    fprintf(fidhsmm,[format,'%d\n'], dataHSMM');   
    
    testNorm(i,:) = observationsNorm;
end

M = ceil(numSeq/2);
M = 0; % <-- no anomalies
testAnom = zeros(M, T);

%generate anomalous sequences
for i=1:M
    
    prevState = randsample(1:Nhid_anom, 1, true, A0_anom);
    prevDur = 1;
    
    currTime = 1;
    durations = [];
    observationsAnom = [];
    states = [];
    
    
    %initialization
    prevDur = randsample(1:Dmax_true, 1, true, D_true(:,prevState, prevDur))';    
    Obs = randsample(1:Nobs, 1, true, O_true(:,prevState))';    
    observationsAnom = [observationsAnom, Obs];    
    
    currTime = currTime + 1;
    
    
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
        prevState = currState;
        
        %save data
        durations = [durations, dur];
        states = [states, currState];
        observationsAnom = [observationsAnom, currObs];
    end
    
    lenObs = length(observationsAnom);
    fprintf(fidhsmm, '%d\n', lenObs);
    
    dataHSMM = [[2:3:3*(lenObs-1)-1 3*(lenObs-1)+1]; observationsAnom-1];
        
    format = repmat('%d\t', 1, lenObs-1);    
    fprintf(fidhsmm, [format,'%d\n'], dataHSMM');        
    
    testAnom(i,:) = observationsAnom;
end

fclose(fidhsmm);
test= [testNorm; testAnom];










