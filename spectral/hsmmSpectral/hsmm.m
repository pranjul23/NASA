%function to train and test spectral model
function [sp len] = hsmm(train, test, Nobs, Nhid_true, Dmin_true, Dmax_true, A1_true, A_true, D_true, O_true, flg)


disp('Spectral Training ...');

[rootTensor  ...
 tailTensor  ...
 obsTensor  ...
 tranTensor ...
 durTensor ] = learnSpectModel2(train, Nobs, Nhid_true, Dmin_true, Dmax_true, A1_true, A_true, D_true, O_true, flg);

disp('Spectral Training ... Done!');

numObs = ceil( (log(Nhid_true)+log(Dmax_true)) / log(Nhid_true) );

% sp = testSpectModel_parallel(test, Nobs, numObs, ...
%                              rootTensor,  ...
%                              tailTensor,  ...
%                              obsTensor,  ...
%                              tranTensor, ...
%                              durTensor);
            
[sp len] = testSpectModel2(test, Nobs, numObs, ...
                     rootTensor,  ...
                     tailTensor,  ...
                     obsTensor,  ...
                     tranTensor, ...
                     durTensor);                           
                     
                     