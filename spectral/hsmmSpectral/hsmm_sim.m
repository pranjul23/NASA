%simulator script to generate parameters and data for HSMM

clear variables;
clc;

%unique ID
ID = 30;

%number of observation symbols
Nobs = 5;


%% =============== TRUE PARAMETERS ========================================

%number of hidden states
Nhid_true = 4;

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
Nhid_init = 4;

%max duration
Dmax_init = 6;

%min duration
Dmin_init = 1;

[A1_init A_init Afull D_init O_init] = genHSMMparam_init(Nobs, Nhid_init, Dmin_init, Dmax_init);


%% =============== ANOMALY PARAMETERS =====================================

%number of hidden states
Nhid_anom = 4;

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
Ttrain = 100;

%number of obseration sequences
numSeq = 1000000;

train = createTraining_sim(Ttrain, numSeq, Nobs, Nhid_true, Dmax_true, A1_true, A_true, D_true, O_true, ID);

%% =============== CREATE ANOMALY SEQUENCE ================================

%span of these observations
span =  Dmax_true;

%find the length of testing sequence for spectral method
Ttest = 100;

%number of obseration sequences
numSeq = 1000;

disp('Create testing start');
test = createTesting_sim(Ttest, numSeq, Nobs, ...
                         Nhid_true, Dmax_true, Nhid_anom, Dmax_anom, ...
                         A1_true, A_true, D_true, O_true, ...
                         A0_anom, A_anom, D_anom, O_anom, ID);

disp('Create testing end');                     
%% ============= SAVE DATA ================================================

save(strcat('SpectralTainData', num2str(ID),'.mat'), ...
     'train', ...
     'test', ...
     'Nobs', ...
     'Nhid_true', ...
     'Dmin_true', ...
     'Dmax_true', ...
     'A1_true', ...
     'A_true', ...
     'D_true', ...
     'O_true');

 
%% =============== RUN HSMM EM ============================================

Ntrain = [500, 1000, 5000, 10000, 100000, 1000000];

origin = pwd;

flg = 0;

cd '../tensor_toolbox';
addpath(pwd);
cd(origin);

for k = 1:length(Ntrain)
    
  cd '../../libdai/examples/'    
    system(['./hsmm_spectest ', num2str(ID), ' 70 ', num2str(Ntrain(k))]);    
  cd(origin);
    
    N = 7;
    
    %true joint likelihood
    loc = strcat('../../libdai/examples/data/HSMMlikelihood_true_', num2str(ID), '.txt');
    tru = load(loc);
    
    %em-computed joint likelihood 
    em = zeros(numSeq, N+1);
    for i=0:N
        loc = strcat('../../libdai/examples/data/HSMMlikelihood_test_', num2str(ID), '-', num2str(i), '.txt');
        em(:,i+1) = load(loc);
    end
    
    numTrain = Ntrain(k);
    hsmm;
    
    save(['emTestResults', num2str(ID), '-', num2str(k), '.mat'], 'tru', 'em', 'sp');        
end








