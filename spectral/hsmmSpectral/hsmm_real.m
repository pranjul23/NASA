%script to assemble data for HSMM testing

clear variables;
clc;

%unique ID
ID = 29;

%number of observation symbols
Nobs = 5;


%% ================ INIT PARAMETERS =======================================

%number of hidden states
Nhid_init = 3;

%max duration
Dmax_init = 6;

%min duration
Dmin_init = 1;

[A1_init A_init Afull D_init O_init] = genHSMMparam_init(Nobs, Nhid_init, Dmin_init, Dmax_init);


%% =============== CREATE FACTOR GRAPH ====================================

createHSMMfactorGraph(Nobs, Nhid_init, Dmax_init, A1_init, A_init, D_init, O_init, ID, 'init');



%% =============== CREATE TRAINING SEQUENCE ===============================

train = createTrainingAir_real(ID);


%% =============== CREATE TESTING SEQUENCE ================================

test = createTestingAir_real(ID);


[rootTensor  ...
 tailTensor  ...
 obsTensor  ...
 tranTensor ...
 durTensor ] = learnSpectModel2(train, ...
                                numObs, ...
                                span, ...
                                skip, ...
                                Nobs, ...
                                Nhid_true, ...
                                Dmin_true, ...
                                Dmax_true, [], [], []);


result = testSpectModel2(test, Nobs, numObs, ...
                         rootTensor,  ...
                         tailTensor,  ...
                         obsTensor,  ...
                         tranTensor, ...
                         durTensor);
