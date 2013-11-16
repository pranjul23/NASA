%simulator script to generate parameters and data for HSMM

clear variables;
clc;

%unique ID
ID = 29;

%number of observation symbols
Nobs = 5;


%% =============== TRUE PARAMETERS ========================================

%number of hidden states
Nhid_true = 3;

%maximum possible duration
Dmax_true = 6;

%minimum possible duration
Dmin_true = 1;

% [A1_true A_true D_true O_true] = genHSMMparam_true(Nobs, Nhid_true, Dmin_true, Dmax_true, ID);
[A1_true A_true D_true O_true] = genHSMMparam_true_dense(Nobs, Nhid_true, Dmin_true, Dmax_true, ID);

if any(isnan(A1_true)) || any(any(isnan(A_true(:,:,1)))) || any(any(isnan(D_true(:,:,1)))) || any(any(isnan(O_true(:,:))))
    error('true param: encountered NaN');
end


%% ================ INIT PARAMETERS =======================================

%number of hidden states
Nhid_init = 3;

%max duration
Dmax_init = 6;

%min duration
Dmin_init = 1;

[A1_init A_init Afull D_init O_init] = genHSMMparam_init(Nobs, Nhid_init, Dmin_init, Dmax_init);


%% =============== ANOMALY PARAMETERS =====================================

%number of hidden states
Nhid_anom = 3;

%number of observation steps
Dmax_anom = 6;

%min duration
%Dmin_anom = Dmax_true;
Dmin_anom = 1;


[A0_anom A_anom D_anom O_anom] = genHSMMparam_anom(Nobs, Nhid_anom, Dmax_anom, Dmin_anom, ...
                                                           A1_true, A_true, D_true, O_true);

                                                       
if any(isnan(A0_anom)) || any(any(isnan(A_anom(:,:,1)))) || any(any(isnan(D_anom(:,:,1)))) || any(any(isnan(O_anom(:,:))))
    error('anomaly param: encountered NaN');
end                


%% =============== CREATE FACTOR GRAPH ====================================

createHSMMfactorGraph(Nobs, Nhid_init, Dmax_init, A1_init, A_init, D_init, O_init, ID, 'init');

createHSMMfactorGraph(Nobs, Nhid_true, Dmax_true, A1_true, A_true, D_true, O_true, ID, 'true');

%% =============== CREATE TRAINING SEQUENCE ===============================

%simulation time
Ttrain = 40;

%number of obseration sequences
numSeq = 1000;

train = createTraining_sim(Ttrain, numSeq, Nobs, Nhid_true, Dmax_true, A1_true, A_true, D_true, O_true, ID);

%% =============== CREATE ANOMALY SEQUENCE ================================


%minimum number of observations needed for spectral method
%(Sec. 7 of Ankur's paper)
numObs = ceil( (log(Nhid_true)+log(Dmax_true)) / log(Nhid_true) );
assert(numObs >= 1);

%span of these observations
span = Nhid_true^(numObs-1);

%number of observations to skip between two neighboring observations
skip = floor(span/(numObs-1)) - 1;
assert(skip >= 1);

%find the length of testing sequence for spectral method
Ttest = Ttrain - span;

%number of obseration sequences
numSeq = 100;

test = createTesting_sim(Ttest, numSeq, Nobs, ...
                         Nhid_true, Dmax_true, Nhid_anom, Dmax_anom, ...
                         A1_true, A_true, D_true, O_true, ...
                         A0_anom, A_anom, D_anom, O_anom, ID);

                     
%% ============= LEARN SPECTRAL MODEL =====================================

[rootTensor  ...
 tailTensor  ...
 obsTensors  ...
 tranTensors ...
 durTensors ] = learnSpectModel(train(1:end, :), numObs, span, skip, Nobs, Nhid_true, Dmin_true, Dmax_true, A1_true, A_true, D_true, O_true);

%index for spectral method
stop_ind = span + 2;
result = testSpectModel(test, Nobs, stop_ind, numObs, ...
                        rootTensor,  ...
                        tailTensor,  ...
                        obsTensors,  ...
                        tranTensors, ...
                        durTensors);
                    
result_eff = testSpectModel_eff(test, Nobs, stop_ind, numObs, ...
                                rootTensor,  ...
                                tailTensor,  ...
                                obsTensors,  ...
                                tranTensors, ...
                                durTensors);                    


[rootTensor  ...
 tailTensor  ...
 obsTensor  ...
 tranTensor ...
 durTensor ] = learnSpectModel2(train(1:end, :), numObs, span, skip, Nobs, Nhid_true, Dmin_true, Dmax_true, A1_true, A_true, D_true, O_true);

result2 = testSpectModel2(test, Nobs, numObs, ...
                          rootTensor,  ...
                          tailTensor,  ...
                          obsTensor,  ...
                          tranTensor, ...
                          durTensor);
                      
[rootTensor  ...
 tailTensor  ...
 obsTensor  ...
 tranTensor ...
 durTensor ] = learnSpectModel3(train(1:end, :), numObs, span, skip, Nobs, Nhid_true, Dmin_true, Dmax_true, A1_true, A_true, D_true, O_true);                      




result3 = testSpectModel2(test, Nobs, numObs, ...
                          rootTensor,  ...
                          tailTensor,  ...
                          obsTensor,  ...
                          tranTensor, ...
                          durTensor);






