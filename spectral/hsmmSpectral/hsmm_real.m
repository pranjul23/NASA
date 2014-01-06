%script to assemble data for HSMM testing

clear variables;
clc;

%unique ID
ID = 29;

%number of observation symbols
Nobs = 10;


%% ================ INIT PARAMETERS =======================================

%number of hidden states
Nhid = 7;

%max duration
Dmax = 45;

%min duration
Dmin = 1;

[A1_init A_init Afull D_init O_init] = genHSMMparam_init(Nobs, Nhid, Dmin, Dmax);


%% =============== CREATE FACTOR GRAPH ====================================

createHSMMfactorGraph(Nobs, Nhid, Dmax, A1_init, A_init, D_init, O_init, ID, 'init');


%% =============== CREATE TRAINING SEQUENCE ===============================

train = createTrainingAir_real(ID);


%% =============== CREATE TESTING SEQUENCE ================================

test = createTestingAir_real(ID);


%% =============== RUN HSMM and EM ========================================

origin = pwd;

% cd '../tensor_toolbox';
% addpath(pwd);
% cd(origin);
% 
% cd '../../libdai/examples/'    
%    system(['./hsmm ', num2str(ID), ' 70']);    
% cd(origin);
% 
% loc = strcat('../../libdai/examples/data/HSMMlikelihood_test_', num2str(ID), '.txt');
% em = load(loc);

sp = hsmm(train, test, Nobs, Nhid, Dmin, Dmax, [], [], [], [], 0);

% save(['RealResults', num2str(ID), '.mat'], 'em', 'sp');