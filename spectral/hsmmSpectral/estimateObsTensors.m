function obsTensors = estimateObsTensors(train, M, obsDim, stateDim, A1_true, A_true, D_true, O_true)

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
        
    obsTensors(i).tensor = sptensor(tensor(res));
    obsTensors(i).var_ind = [i -i];
    
    %obsTensors{i} = scaled_tensor*((V/(U'*scaled_tensor*V))*U');        
end
 
% P = computeO1O2(A1_true, A_true, D_true(:,:,1), D_true, O_true);
% % [U S V] = svd(P);
% [U S V] = svds(sparse(P), stateDim);
% 
% % num = nnz(diag(S)>=1e-30);
% % U = U(:,1:num);
% % S = S(1:num, 1:num);
% % V = V(:,1:num);
% 
% res = P*(V*diag(1./diag(S))*U');
% obsTensors(1).tensor = sptensor(tensor(res));
% obsTensors(1).var_ind = [1 -1];
% 
% 
% P = computeO2O3(A1_true, A_true, D_true(:,:,1), D_true, O_true);
% % [U S V] = svd(P);
% [U S V] = svds(sparse(P), stateDim);
% 
% % num = nnz(diag(S)>=1e-30);
% % U = U(:,1:num);
% % S = S(1:num, 1:num);
% % V = V(:,1:num);
% 
% res = P*(V*diag(1./diag(S))*U');
% obsTensors(2).tensor = sptensor(tensor(res));
% obsTensors(2).var_ind = [2 -2];
% 
% 
% P = computeO3O4(A1_true, A_true, D_true(:,:,1), D_true, O_true);
% % [U S V] = svd(P);
% [U S V] = svds(sparse(P), stateDim);
% 
% % num = nnz(diag(S)>=1e-30);
% % U = U(:,1:num);
% % S = S(1:num, 1:num);
% % V = V(:,1:num);
% 
% res = P*(V*diag(1./diag(S))*U');
% obsTensors(3).tensor = sptensor(tensor(res));
% obsTensors(3).var_ind = [3 -3];
% 
% 
% P = computeO4O5(A1_true, A_true, D_true(:,:,1), D_true, O_true);
% % [U S V] = svd(P);
% [U S V] = svds(sparse(P), stateDim);
% 
% % num = nnz(diag(S)>=1e-30);
% % U = U(:,1:num);
% % S = S(1:num, 1:num);
% % V = V(:,1:num);
% 
% res = P*(V*diag(1./diag(S))*U');
% obsTensors(4).tensor = sptensor(tensor(res));
% obsTensors(4).var_ind = [4 -4];
% 
% 
% 
% P = computeO5O6(A1_true, A_true, D_true(:,:,1), D_true, O_true);
% % [U S V] = svd(P);
% [U S V] = svds(sparse(P), stateDim);
% 
% % num = nnz(diag(S)>=1e-30);
% % U = U(:,1:num);
% % S = S(1:num, 1:num);
% % V = V(:,1:num);
% 
% res = P*(V*diag(1./diag(S))*U');
% obsTensors(5).tensor = sptensor(tensor(res));
% obsTensors(5).var_ind = [5 -5];

 
 
% P = computeO6O7(A1_true, A_true, D_true(:,:,1), D_true, O_true);
% % [U S V] = svd(P);
% [U S V] = svds(sparse(P), stateDim);
% 
% % num = nnz(diag(S)>=1e-30);
% % U = U(:,1:num);
% % S = S(1:num, 1:num);
% % V = V(:,1:num);
% 
% res = P*(V*diag(1./diag(S))*U');
% obsTensors(6).tensor = sptensor(tensor(res));
% obsTensors(6).var_ind = [6 6];






