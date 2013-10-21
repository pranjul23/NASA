function obsTensors = estimateObsTensors(train, M, obsDim, stateDim)

%sequence length
[L, N] = size(train);

%compute linear indeces for efficient matrix access
m = 2;
P = ones(m,1);
for k = 2:m
    P(k) = P(k-1)*obsDim;
end

for i=1:M
    
    tensr = zeros(obsDim,obsDim); %modes i i+1
    
    %compute indeces
    ind = train(:, [i i+1]);
    ind = [ind(:,1)  ind(:, 2:m)-ones(size(ind(:, 2:m)))];
    
    ind = ind*P;
    
    %compute empirical estimate of probabilities    
    count = histc(ind, unique(ind));
    tensr(unique(ind)) = count;
            
    scaled_tensor = tensr/L;
    
    %do svd by exracting only the K largest values/vectors
    [U S V] = svds(sparse(scaled_tensor), stateDim);
    
    %drop singular values with small values
    num = nnz(diag(S)>=1e-3);
    assert(num>=1);
    
    U = U(:,1:num);
    S = S(1:num, 1:num);
    V = V(:,1:num);
    
    %eliminate mode i+1
    
    %U corresponds to i
    %V corresponds to i+1
    
    %So, we need USV' * VS^U' = UU' (like i over i)
    
    %obsTensors{i} = scaled_tensor*(V*inv(U'*scaled_tensor*V)*U');
    
    res = scaled_tensor*(V*diag(1./diag(S))*U');
    
    obsTensors(i).tensor = sptensor(tensor(res));
    obsTensors(i).var_ind = [i i];
    
    %obsTensors{i} = scaled_tensor*((V/(U'*scaled_tensor*V))*U');        
end