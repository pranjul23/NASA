function durTensors = estimateDurTensors(train, iter_ind, numObs, obsDim, stateDim, durMin, durMax)

%sequences size
[L, N] = size(train);


%compute linear indeces for efficient matrix access
m = 2*numObs;
P = ones(m,1);
for k = 2:m
    P(k) = P(k-1)*obsDim;
end

for k = 1:length(iter_ind)
    
    i = iter_ind(k);
    
    tens = zeros(obsDim*(ones(1, 2*numObs)));
    
    left_start = i - durMin - numObs;
    left_end = i - 1 - durMin;
    
    right_start = i - 1 + durMax;
    right_end = i - 1 + durMax + numObs - 1;
    
    ind = [left_start:left_end  right_start:right_end];    
    
    ind = train(:, ind);
    ind = [ind(:,1)  ind(:, 2:m)-ones(size(ind(:, 2:m)))];    
    ind = ind*P;
    
    %compute empirical estimate of probabilities    
    count = histc(ind, unique(ind));
    tens(unique(ind)) = count;
    
    
    %scale the computed multdim array and convert in to tensor class
    scaled_tensor = sptensor(tens/L);
    scaled_tensor_inv = tensor(tens/L);
    
    %tensor's modes - 1:numObs,1:numObs
    %eliminate first set of numObs modes
    
    %modes to be eliminated are put as columns
    
    mat_tensor = tenmat(scaled_tensor_inv, 1:numObs, 't');
    
    %do svd by exracting only the K largest values/vectors
    [U S V] = svds(sparse(mat_tensor.data), durMax*stateDim);
            
    %drop singular values with small values
    num = nnz(diag(S) >= 1e-3);
    assert(num >= 1);
    
    U = U(:,1:num);
    S = S(1:num, 1:num);
    V = V(:,1:num);
    
    %U corresponds to second set of numObs modes
    %V corresponds to first set of numObs modes
    
    %So, we need USV' * VS^U' = UU'
        
    mat_tensor(:) = (V*diag(1./diag(S))*U')';
    %mat_tensor(:) = (V*inv(U'*mat_tensor.data*V)*U')';
        
        
    res = ttt(scaled_tensor, sptensor(tensor(mat_tensor)), 1:numObs);    
    
    %% TODO: make it symmetric !!
    durTensors(iter_ind(k)).tensor = res;
    durTensors(iter_ind(k)).var_ind = [right_start:right_end right_start:right_end];    
end



