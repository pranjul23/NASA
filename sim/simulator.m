%simulator script to generate data for 
%explicit duration HSMM

%OUPUT:
%aray struct data
%data(i).states = sequence of states in run i
%data(i).obs = sequence of observations in run i

clear all;

%simulation time
T = 100;

%number of trials
num = 10;

%number of hidden states
N = 5;

%number of observation symbols
%including NULL
M = 10;

%initial state distribution
pi = rand(N,1);
pi = pi./sum(pi); %normalize


%transit distribution
%self-transitions are not allowed
A = rand(N);
for i=1:N
    A(i,i)=0;
    A(i,:) = A(i,:)/sum(A(i,:));
end

%observation distribution
O = rand(M,N);
for i=1:N
    O(:,i) = O(:,i)/sum(O(:,i));
end

%duration distribution
%it is obtained from duration.m

%allocate memory
data(num,1).obs=[];
data(num,1).states=[];

%do simulation
for i=1:num
    
    currState = randsample(1:N, 1, true, pi);
    currTime = 0;
    observations = [];
    states = currState;
    
    while(currTime < T)
        %select duration (# of timesteps) of current state
        dur = round(duration(currState, N));
        
        %perform dur observation in current state
        currObs = randsample(1:M, dur, true, O(:,currState))';
        observations = [observations; currObs];
        
        %select next state
        currState = randsample(1:N, 1, true, A(currState,:))';
        states = [states; currState];
        
        %advance time
        currTime = currTime + dur;
    end
    
    data(i).obs = observations;
    data(i).states = states;
end


