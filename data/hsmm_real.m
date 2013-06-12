%script to generate structure and data for
%explicit duration HSMM

clear;
clc;


%unique ID
ID = 15;


%% ================ INIT PARAMETERS =======================================

%number of observation symbols
%including NULL
Nobs = 41;

%number of hidden states
Nhid = 12;

%max duration
Dmax = 150;

[A0 D0 A D O] = genHSMMparam_init(Nobs, Nhid, Dmax);

[A0_hmm A_hmm O_hmm] = genHMMparam(A0, A, O);



%% =============== CREATE FACTOR GRAPH ====================================

createHSMMfactorGraph(Nobs, Nhid, Dmax, A0, D0, A, D, O, ID);

createHMMfactorGraph(Nobs, Nhid, A0_hmm, A_hmm, O_hmm, ID);



%% =============== CREATE TRAINING SEQUENCE ===============================

trainInd = createTraining_real(ID);




%% =============== CREATE TESTING SEQUENCE ================================

createTesting_real(ID, trainInd);