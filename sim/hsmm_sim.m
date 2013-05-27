%simulator script to generate parameters and data for HSMM/HMM

clear;
clc;

%number of observation symbols
Nobs = 7;


%% =============== TRUE PARAMETERS ========================================

%number of hidden states
Nhid_true = 4;

%number of observation steps
Dmax_true = 10;


[A0_true D0_true A_true D_true O_true] = genHSMMparam_true(Nobs, Nhid_true, Dmax_true);


%% ================ INIT PARAMETERS =======================================

%number of hidden states
Nhid_init = 4;

%number of observation steps
Dmax_init = 10;

[A0_init D0_init A_init D_init O_init] = genHSMMparam_init(Nobs, Nhid_init, Dmax_init);

A0_init = A0_true;
D0_init = D0_true;
A_init = A_true;
D_init = D_true;
O_init = O_true;


[A0_init_hmm A_init_hmm O_init_hmm] = genHMMparam(Nobs, Nhid_init);


%% =============== ANOMALY PARAMETERS =====================================

%number of hidden states
Nhid_anom = 4;

%the larger the Nhid the larger the spread likelihood values of data
%the change of Nhid_anom wrt true Nhid_true does not help in detecting
%anomaly


%number of observation steps
Dmax_anom = 50;


%min duration
Dmin_anom = Dmax_true;


[A0_anom D0_anom A_anom D_anom O_anom] = genHSMMparam_anom(Nobs, Nhid_anom, Dmax_anom, Dmin_anom, ...
                                                           A0_true, D0_true, A_true, D_true, O_true);


%% =============== CREATE FACTOR GRAPH ====================================

createHSMMfactorGraph(Nobs, Nhid_init, Dmax_init, A0_init, D0_init, A_init, D_init, O_init);

createHMMfactorGraph(Nobs, Nhid_init, A0_init_hmm, A_init_hmm, O_init_hmm);


%% =============== CREATE TRAINING SEQUENCE ===============================

%simulation time
T = 100;

%number of obseration sequences
numSeq = 100;

createTraining(T, numSeq, Nobs, Nhid_true, Dmax_true, A0_true, D0_true, A_true, D_true, O_true);


%% =============== CREATE ANOMALY SEQUENCE ================================

%simulation time
T = 200;

%number of obseration sequences
numSeq = 1000;

createTesting(T, numSeq, Nobs, ...
              Nhid_true, Dmax_true, Nhid_anom, Dmax_anom, ...
              A0_true, D0_true, A_true, D_true, O_true, ...
              A0_anom, D0_anom, A_anom, D_anom, O_anom);
















