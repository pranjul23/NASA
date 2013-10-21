function rootTensor = estimateRootTensor(train, numObs, obsDim, durMax)

%sequences size
[L, N] = size(train);

%compute linear indeces for efficient matrix access
m = 2 + numObs;
P = ones(m,1);
for k = 2:m
    P(k) = P(k-1)*obsDim;
end


tens = zeros(obsDim*(ones(1, numObs+2)));

ind = [1 2 2+durMax : 2+durMax+numObs-1];

ind = train(:, ind);
ind = [ind(:,1)  ind(:, 2:m)-ones(size(ind(:, 2:m)))];
ind = ind*P;

%compute empirical estimate of probabilities
count = histc(ind, unique(ind));
tens(unique(ind)) = count;

rootTensor.tensor = sptensor(tensor(tens/L));
rootTensor.var_ind = [1 2 2+durMax : 2+durMax+numObs-1];