function [] = createHMMtraining(T, numSeq, Nobs, Nhid, A0, A, O)


data_seq = cell(numSeq,1);

fid = fopen('../libdai/examples/HMMtraining.txt', 'w');

fprintf(fid, '%d\n',numSeq);


for i=1:numSeq
    
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
        
        %perform dur observations in current state
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
    
    fprintf(fidhsmm, [format,'%d\n'], dataHSMM'); 
    fprintf(fidhmm, [format,'%d\n'], dataHMM');
    
    %save data for MATLAB HMM function
    data_seq{i} = observations + ones(size(observations));    
end

save('matlabHMMdata.mat', 'data_seq');

fclose(fid);