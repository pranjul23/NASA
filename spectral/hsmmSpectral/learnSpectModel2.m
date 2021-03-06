function [rootTensor ...
          tailTensor ...
          obsTensor ...
          tranTensor ...
          durTensor] = learnSpectModel2(train, ...
                                        obsDim, ...
                                        stateDim, ...
                                        durMin, ...
                                        durMax, A1_true, A_true, D_true, O_true, flg)




                                    
%check to make sure we have enough data to learn model                                    
[L, N] = size(train);

%assumption 
assert(obsDim >= stateDim);

%minimum number of observations needed for spectral method
%(Sec. 7 of Ankur's paper)
numObs = ceil( (log(stateDim)+log(durMax)) / log(stateDim) );
assert(numObs >= 1);

%span of these observations
span =  durMax;


%estimate root tensor
rootTensor = estimateRootTensor(train, ...
                                numObs,...
                                obsDim,...
                                stateDim, ...                                                                
                                durMax, A1_true, A_true, D_true, O_true, flg);
                     
                       
%estimate observation tensor
obsTensor = estObsTensor(train, ...
                         obsDim, ...
                         stateDim, A1_true, A_true, D_true, O_true, flg);

%estimate transition tensor 
tranTensor = estTranTensor(train, ...
                           numObs, ...
                           span, ...
                           obsDim, ...
                           stateDim, ...
                           durMin, ...
                           durMax, A1_true, A_true, D_true, O_true, flg);

%estimate duration tensor
durTensor = estDurTensor(train, ...
                         numObs, ...
                         span, ...
                         obsDim, ...
                         stateDim, ...
                         durMin, ...
                         durMax, A1_true, A_true, D_true, O_true, flg);
                     
                     

%estimate tail tensor
tmp = tensor(ones(obsDim,1));

%marginalize modes 2:numObs which do not exist in tailTensor
tailTensor = squeeze(ttt(tranTensor, tmp, 2, 1));
for i = 2:numObs
    tailTensor = squeeze(ttt(tailTensor, tmp, 2, 1));
end
