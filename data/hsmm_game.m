%simulator script to generate parameters and data for HSMM/HMM

clear;
clc;

%unique ID
ID = 15;


%number of observation symbols
Nobs = 12;

%% ================ INIT PARAMETERS =======================================

%number of hidden states
Nhid_init = 16;

%number of observation steps
Dmax_init = 100;

[A0_init D0_init A_init Afull D_init O_init] = genHSMMparam_init(Nobs, Nhid_init, Dmax_init);

[A0_init_hmm A_init_hmm O_init_hmm] = genHMMparam(A0_init, Afull, O_init);


%% =============== CREATE FACTOR GRAPH ====================================

createHSMMfactorGraph(Nobs, Nhid_init, Dmax_init, A0_init, D0_init, A_init, D_init, O_init, ID);

createHMMfactorGraph(Nobs, Nhid_init, A0_init_hmm, A_init_hmm, O_init_hmm, ID);


%% =============== CREATE TRAINING SEQUENCE ===============================

createTraining_game(ID);

%% =============== CREATE TESTING SEQUENCE ===============================

createTesting_game(ID);











