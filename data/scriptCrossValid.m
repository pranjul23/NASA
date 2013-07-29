%script to compare performance of HMM and HSMM using cross validation

clear;
clc;

%number of observation symbols
Nobs = 12;

%maximum duration in any state
Dmax_init = 100;

%unique ID
ID = 0;


%iterate through different number of hidden states
for Nhid_init = 2:16
    
    %% ================ INIT PARAMETERS =======================================
    [A0_init D0_init A_init Afull D_init O_init] = genHSMMparam_init(Nobs, Nhid_init, Dmax_init);
    
    [A0_init_hmm A_init_hmm O_init_hmm] = genHMMparam(A0_init, Afull, O_init);
    
    
    %% =============== CROSS VALIDATION =======================================
    indices = crossvalind('Kfold', 110, 11);
    for i = 1:11
        test = find(indices == i); train = find(indices ~= i);
        
        %% =============== CREATE FACTOR GRAPH ====================================
        createHSMMfactorGraph(Nobs, Nhid_init, Dmax_init, A0_init, D0_init, A_init, D_init, O_init, ID);
        
        createHMMfactorGraph(Nobs, Nhid_init, A0_init_hmm, A_init_hmm, O_init_hmm, ID);
        
        %% =============== CREATE TRAINING AND TESTING DATA ====================================
        createTrainingData_game(train, ID);
        createTestingData_game(test, ID);
        
        ID = ID + 1;
    end
end
