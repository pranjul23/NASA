%simulator script to generate parameters and data for HSMM/HMM

clear;
clc;

%unique ID
ID = 0;

%number of observation symbols
Nobs = 12;


%% =============== TRUE PARAMETERS ========================================

%number of hidden states
Nhid_true = 12;

%number of observation steps
Dmax_true = 60;


[A0_true D0_true A_true D_true O_true] = genHSMMparam_true(Nobs, Nhid_true, Dmax_true, ID);


%% ================ INIT PARAMETERS =======================================

%number of hidden states
Nhid_init = 10;

%number of observation steps
Dmax_init = 100;

[A0_init D0_init A_init Afull D_init O_init] = genHSMMparam_init(Nobs, Nhid_init, Dmax_init);

[A0_init_hmm A_init_hmm O_init_hmm] = genHMMparam(A0_init, Afull, O_init);


%% =============== ANOMALY PARAMETERS =====================================

%number of hidden states
Nhid_anom = 4;

%number of observation steps
Dmax_anom = 10;


%min duration
Dmin_anom = Dmax_true;


[A0_anom D0_anom A_anom D_anom O_anom] = genHSMMparam_anom(Nobs, Nhid_anom, Dmax_anom, Dmin_anom, ...
                                                           A0_true, D0_true, A_true, D_true, O_true);


%% =============== CREATE FACTOR GRAPH ====================================

createHSMMfactorGraph(Nobs, Nhid_init, Dmax_init, A0_init, D0_init, A_init, D_init, O_init, ID);
createHSMMfactorGraph_true(Nobs, Nhid_true, Dmax_true, A0_true, D0_true, A_true, D_true, O_true, ID);

createHMMfactorGraph(Nobs, Nhid_init, A0_init_hmm, A_init_hmm, O_init_hmm, ID);


%% =============== CREATE TRAINING SEQUENCE ===============================

%simulation time
T = 100;

%number of obseration sequences
numSeq = 200;

createTraining_sim(T, numSeq, Nobs, Nhid_true, Dmax_true, A0_true, D0_true, A_true, D_true, O_true, ID);

%% =============== CREATE ANOMALY SEQUENCE ================================

%simulation time
T = 100;

%number of obseration sequences
numSeq = 20;

createTesting_sim(T, numSeq, Nobs, ...
                  Nhid_true, Dmax_true, Nhid_anom, Dmax_anom, ...
                  A0_true, D0_true, A_true, D_true, O_true, ...
                  A0_anom, D0_anom, A_anom, D_anom, O_anom, ID);
















