function rootTensor = estimateRootTensor(train, numObs, skip, span, obsDim, durMax, A1_true, A_true, D_true, O_true)

%sequences size
[L, N] = size(train);

%compute linear indeces for efficient matrix access
m = 2 + numObs;
P = ones(m,1);
for k = 2:m
    P(k) = P(k-1)*obsDim;
end


tens = zeros(obsDim*(ones(1, numObs+2)));

indices = [1 2 3:skip:(3+span-1)];

ind = train(:, indices);
ind = [ind(:,1)  ind(:, 2:m)-ones(size(ind(:, 2:m)))];
ind = ind*P;

%compute empirical estimate of probabilities
uniq_ind = unique(ind);
count = histc(ind, uniq_ind);
tens(uniq_ind) = count;

rootTensor.tensor = sptensor(tensor(tens/L));
% rootTensor.tensor = sptensor(tensor(computeO1O2O3O4(A1_true, A_true, D_true(:,:,1), D_true, O_true)));

% M = computeO3O4D1X2(A1_true, A_true, D_true(:,:,1), D_true, O_true);
% M = tensor(M);
% eA = tensor(embedA1(A_true, D_true(:,:,1), A1_true));
% res2 = ttt(eA, M, [4 1], [3 4]);
% res2 = ttt(tensor(O_true), res2, 2, 1);
% res2 = ttt(tensor(O_true), res2, 2, 2);
% 
% rootTensor.tensor = sptensor(res2);

rootTensor.var_ind = indices;