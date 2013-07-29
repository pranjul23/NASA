%simulator script to generate parameters and data for HSMM/HMM

clear;
clc;

%number of observation steps
Dmax_init = 100;

%number of observation symbols
Nobs = 12;

%unique ID
ID = 0;


%iterate through different number of hidden states
for Nhid_init = 11:15
        
    %different random initializations
    for i=1:4
        %% ================ INIT PARAMETERS =======================================
        
        [A0_init D0_init A_init Afull D_init O_init] = genHSMMparam_init(Nobs, Nhid_init, Dmax_init);
        
        [A0_init_hmm A_init_hmm O_init_hmm] = genHMMparam(A0_init, Afull, O_init);
        
        
        %% =============== CREATE FACTOR GRAPH ====================================
        
        createHSMMfactorGraph(Nobs, Nhid_init, Dmax_init, A0_init, D0_init, A_init, D_init, O_init, ID);
        
        createHMMfactorGraph(Nobs, Nhid_init, A0_init_hmm, A_init_hmm, O_init_hmm, ID);
        
        
        %% =============== CREATE TRAINING SEQUENCE ===============================
        
        createTraining_game(ID);
        
        %% =============== CREATE TESTING SEQUENCE ===============================
        
        createTesting_game(ID);
        
        ID = ID + 1;
    end    
end






