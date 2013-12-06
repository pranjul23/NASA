%simulator script to generate parameters and data for HSMM

clear variables;
clc;

%unique ID
ID = 29;

%number of observation symbols
Nobs = 5;


%% =============== TRUE PARAMETERS ========================================

%number of hidden states
Nhid_true = 4;

%maximum possible duration
Dmax_true = 3;

%minimum possible duration
Dmin_true = 1;

% [A1_true A_true D_true O_true] = genHSMMparam_true(Nobs, Nhid_true, Dmin_true, Dmax_true, ID);
[A1_true A_true D_true O_true] = genHSMMparam_true_dense(Nobs, Nhid_true, Dmin_true, Dmax_true, ID);

if any(isnan(A1_true)) || any(any(isnan(A_true(:,:,1)))) || any(any(isnan(D_true(:,:,1)))) || any(any(isnan(O_true(:,:))))
    error('true param: encountered NaN');
end


%% ================ INIT PARAMETERS =======================================

%number of hidden states
Nhid_init = 4;

%max duration
Dmax_init = 3;

%min duration
Dmin_init = 1;

[A1_init A_init Afull D_init O_init] = genHSMMparam_init(Nobs, Nhid_init, Dmin_init, Dmax_init);


%% =============== ANOMALY PARAMETERS =====================================

%number of hidden states
Nhid_anom = 4;

%number of observation steps
Dmax_anom = 3;

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
Ttrain = 50;

%number of obseration sequences
numSeq = 100000;

train = createTraining_sim(Ttrain, numSeq, Nobs, Nhid_true, Dmax_true, A1_true, A_true, D_true, O_true, ID);

%% =============== CREATE ANOMALY SEQUENCE ================================

%span of these observations
span =  Dmax_true;

%find the length of testing sequence for spectral method
Ttest = Ttrain - span;

%number of obseration sequences
numSeq = 50;

test = createTesting_sim(Ttest, numSeq, Nobs, ...
                         Nhid_true, Dmax_true, Nhid_anom, Dmax_anom, ...
                         A1_true, A_true, D_true, O_true, ...
                         A0_anom, A_anom, D_anom, O_anom, ID);

                     
%% ============= LEARN SPECTRAL MODEL =====================================

% [rootTensor  ...
%  tailTensor  ...
%  obsTensors  ...
%  tranTensors ...
%  durTensors ] = learnSpectModel(train(1:end, :), numObs, span, skip, Nobs, Nhid_true, Dmin_true, Dmax_true, A1_true, A_true, D_true, O_true);
% 
% %index for spectral method
% stop_ind = span + 2;
% result = testSpectModel(test, Nobs, stop_ind, numObs, ...
%                         rootTensor,  ...
%                         tailTensor,  ...
%                         obsTensors,  ...
%                         tranTensors, ...
%                         durTensors);
%                     
% result_eff = testSpectModel_eff(test, Nobs, stop_ind, numObs, ...
%                                 rootTensor,  ...
%                                 tailTensor,  ...
%                                 obsTensors,  ...
%                                 tranTensors, ...
%                                 durTensors);                    


[rootTensor  ...
 tailTensor  ...
 obsTensor  ...
 tranTensor ...
 durTensor ] = learnSpectModel2(train, Nobs, Nhid_true, Dmin_true, Dmax_true, A1_true, A_true, D_true, O_true);

numObs = ceil( (log(Nhid_true)+log(Dmax_true)) / log(Nhid_true) );

result3 = testSpectModel2(test, Nobs, numObs, ...
                          rootTensor,  ...
                          tailTensor,  ...
                          obsTensor,  ...
                          tranTensor, ...
                          durTensor);





