function obsTensors = estObsTensor2(train, M, obsDim, stateDim, A1_true, A_true, D_true, O_true)

%sequence length
[L, N] = size(train);

%compute linear indeces for efficient matrix access
m = 2;
P = ones(m,1);
for k = 2:m
    P(k) = P(k-1)*obsDim;
end

res_acc = zeros(obsDim,obsDim);


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
%     [U S V] = svd(scaled_tensor);
    
    %drop singular values with small values
%     num = nnz(diag(S)>=1e-5);
    num = nnz(diag(S)>=1e-20);
    
    assert(num>=1);
%     
    U = U(:,1:num);
    S = S(1:num, 1:num);
    V = V(:,1:num);
    
    %eliminate mode i+1
    
    %U corresponds to i
    %V corresponds to i+1
    
    %So, we need USV' * VS^U' = UU' (like i over i)
    
    %obsTensors{i} = scaled_tensor*(V*inv(U'*scaled_tensor*V)*U');
    
    res = scaled_tensor*(V*diag(1./diag(S))*U');
        
    res_acc = res_acc + res;
end

obsTensors =  sptensor(tensor(res_acc/M));
    
 
