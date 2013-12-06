function durTensors = estimateDurTensors(train, ...
                                         iter_ind, ...
                                         numObs,...
                                         skip, ...
                                         span, ...
                                         obsDim, ...
                                         stateDim, ...
                                         durMin, ...
                                         durMax, A1_true, A_true, D_true, O_true)

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
    
    %compute indeces
    left_start = i - 1 - (span-1);
    left_end = i - 1;
    
    right_start = i + 2;
    right_end = i + span + 1;
    
    ind = [left_start:skip:left_end  right_start:skip:right_end];    
    
    ind_res = right_start:skip:right_end;
    
    ind = train(:, ind);
    ind = [ind(:,1)  ind(:, 2:m)-ones(size(ind(:, 2:m)))];    
    ind = ind*P;
    
    %compute empirical estimate of probabilities    
    count = histc(ind, unique(ind));
    tens(unique(ind)) = count;
    
    
    %scale the computed multdim array and convert in to tensor class
    scaled_tensor = sptensor(tens/L);
    
    if k==1
        scaled_tensor = sptensor(computeO1O2O5O6(A1_true, A_true, D_true(:,:,1), D_true, O_true));
    else
        scaled_tensor = sptensor(computeO2O3O6O7(A1_true, A_true, D_true(:,:,1), D_true, O_true));
    end

    
    
    %compute indeces for inverse tensor
    left_start = i - 1 - (span-1);
    left_end = i - 1;
    
    right_start = i + 1;
    right_end = i + span;
    
    ind = [left_start:skip:left_end  right_start:skip:right_end];    
    
    ind_res = [ind_res right_start:skip:right_end];
    
    ind = train(:, ind);
    ind = [ind(:,1)  ind(:, 2:m)-ones(size(ind(:, 2:m)))];    
    ind = ind*P;
    
    %compute empirical estimate of probabilities    
    count = histc(ind, unique(ind));
    tens(unique(ind)) = count;
    
    
    scaled_tensor_inv = tensor(tens/L);
    
    if k==1
        scaled_tensor_inv = sptensor(computeO1O2O4O5(A1_true, A_true, D_true(:,:,1), D_true, O_true));
    else
        scaled_tensor_inv = sptensor(computeO2O3O5O6(A1_true, A_true, D_true(:,:,1), D_true, O_true));
    end
    

    %tensor's modes - 1:numObs,1:numObs
    %eliminate first set of numObs modes
    
    %modes to be eliminated are put as columns
    
    mat_tensor = tenmat(scaled_tensor_inv, 1:numObs, 't');
    
    %do svd by exracting only the K largest values/vectors
    [U S V] = svds(sparse(mat_tensor.data), durMax*stateDim);
%     [U S V] = svd(mat_tensor.data);
    
    %drop singular values with small values
%     num = nnz(diag(S) >= 1e-5);
    num = nnz(diag(S) >= 1e-30);
    assert(num >= 1);
    
    U = U(:,1:num);
    S = S(1:num, 1:num);
    V = V(:,1:num);
%     
    %U corresponds to second set of numObs modes
    %V corresponds to first set of numObs modes
    
    %So, we need USV' * VS^U' = UU'
        
    mat_tensor(:) = (V*diag(1./diag(S))*U')';
    %mat_tensor(:) = (V*inv(U'*mat_tensor.data*V)*U')';
        
        
    res = ttt(scaled_tensor, sptensor(tensor(mat_tensor)), 1:numObs);    
    
    durTensors(iter_ind(k)).tensor = res;
    durTensors(iter_ind(k)).var_ind = ind_res; 
end


iter_ind = [2 3 4];

for k = 1:length(iter_ind)

    
%     if k==1
%         M = computeO3O4D1X2(A1_true, A_true, D_true(:,:,1), D_true, O_true);
%         M = tensor(M);
%         Mm = tenmat(M,[1 2]);
%         invMm = pinv(Mm.data);
%         Mm(:) = invMm';
%         M = tensor(Mm);
%         eD = tensor(embedD(D_true));
%         V = computeO4O5D2X2(A1_true, A_true, D_true(:,:,1), D_true, O_true);
%         V = tensor(V);
%         res2 = ttt(M, eD, [3 4], [4 3]);
%         res2 = ttt(res2, V, [3 4], [3 4]);        
%         
%         res2m = tenmat(res2, [1 2]);
%         res2m(:) = res2m.data';
%         res = sptensor(tensor(res2m));
%         
%         ind_res = [4 5 3 4];
    if k==1
        M = computeO4O5D2X3(A1_true, A_true, D_true(:,:,1), D_true, O_true);
        M = tensor(M);
        Mm = tenmat(M,[1 2]);
        invMm = pinv(Mm.data);
        Mm(:) = invMm';
        M = tensor(Mm);
        eD = tensor(embedD(D_true));
        V = computeO5O6D3X3(A1_true, A_true, D_true(:,:,1), D_true, O_true);
        V = tensor(V);
        res2 = ttt(M, eD, [3 4], [4 3]);
        res2 = ttt(res2, V, [3 4], [3 4]);
         
        res2m = tenmat(res2, [1 2]);
        res2m(:) = res2m.data';
        res = sptensor(tensor(res2m));
        
        ind_res = [5 6 4 5];
    elseif k==2
        M = computeO5O6D3X4(A1_true, A_true, D_true(:,:,1), D_true, O_true);
        M = tensor(M);
        Mm = tenmat(M,[1 2]);
        invMm = pinv(Mm.data);
        Mm(:) = invMm';
        M = tensor(Mm);
        eD = tensor(embedD(D_true));
        V = computeO6O7D4X4(A1_true, A_true, D_true(:,:,1), D_true, O_true);
        V = tensor(V);
        res2 = ttt(M, eD, [3 4], [4 3]);
        res2 = ttt(res2, V, [3 4], [3 4]);  
        
        res2m = tenmat(res2, [1 2]);
        res2m(:) = res2m.data';
        res = sptensor(tensor(res2m));
        
        ind_res = [6 7 5 6];
    end
    
    
    durTensors(iter_ind(k)).tensor = res;
    durTensors(iter_ind(k)).var_ind = ind_res; 

end



















