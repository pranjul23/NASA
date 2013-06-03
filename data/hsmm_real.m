%script to generate structure and data for
%explicit duration HSMM

clear;
clc;


%unique ID
ID = 0;


%% ================ INIT PARAMETERS =======================================

%number of observation symbols
%including NULL
Nobs = 183;

%number of hidden states
Nhid = 5;

%max duration
Dmax = 50;

[A0 D0 A D O] = genHSMMparam_init(Nobs, Nhid, Dmax);

[A0_hmm A_hmm O_hmm] = genHMMparam(A0, A, O);



%% =============== CREATE FACTOR GRAPH ====================================

createHSMMfactorGraph(Nobs, Nhid, Dmax, A0, D0, A, D, O, ID);

createHMMfactorGraph(Nobs, Nhid, A0_hmm, A_hmm, O_hmm, ID);



%% =============== CREATE TRAINING SEQUENCE ===============================

createTraining_real(ID);




%% =============== CREATE TESTING SEQUENCE ================================

createTesting_real(ID);