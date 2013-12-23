function rootTensor = estimateRootTensor(train, ...
                                         numObs, ...
                                         obsDim, ...
                                         stateDim, ...
                                         durMax, A1_true, A_true, D_true, O_true, flg)

%sequences size
[L, N] = size(train);

%compute linear indeces for efficient matrix access
m = 2 + numObs;
P = ones(m,1);
for k = 2:m
    P(k) = P(k-1)*obsDim;
end

tens = zeros(obsDim*(ones(1, numObs+2)));

indices = [1 2 getRightInd(3, stateDim, durMax, numObs)];

ind = train(:, indices);
ind = [ind(:,1)  ind(:, 2:m)-ones(size(ind(:, 2:m)))];
ind = ind*P;

%compute empirical estimate of probabilities
uniq_ind = unique(ind);
count = histc(ind, uniq_ind);
tens(uniq_ind) = count;

rootTensor = tensor(tens/L);
% rootTensor.tensor = sptensor(tensor(computeO1O2O3O4(A1_true, A_true, D_true(:,:,1), D_true, O_true)));

if flg
    L = getLeftCondProb(O_true, D_true, A_true, A1_true, [1 2 2], 0);
    R = getRightCondProb(O_true, D_true, A_true, [2  getRightInd(3, stateDim, durMax, numObs)], 0);
    
    rootTensor = ttt(L, R, [3 4], [length(R.size)-1 length(R.size)]);
end


% M = getRightCondProb(O_true, D_true, A_true, [2 getRightInd(3, stateDim, durMax, numObs)], 1);
% eA = tensor(embedA1(A_true, D_true(:,:,1), A1_true));
% 
% res2 = ttt(eA, M, [1 4], [length(M.size)-1 length(M.size)]);
% res2 = ttt(tensor(O_true), res2, 2, 1);
% res2 = ttt(tensor(O_true), res2, 2, 2);
% 
% rootTensor = res2;



