function durTensor = estDurTensor(train, ...
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

scaled_tens = zeros(obsDim*(ones(1, 2*numObs)));
scaled_tens_inv = zeros(obsDim*(ones(1, 2*numObs)));

for k = 1:length(iter_ind)
    
    i = iter_ind(k);
    
    tens = zeros(obsDim*(ones(1, 2*numObs)));
    tens_inv = zeros(obsDim*(ones(1, 2*numObs)));
    
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
            
    
    %compute indeces for inverse tensor
    left_start = i - 1 - (span-1);
    left_end = i - 1;
    
    right_start = i + 1;
    right_end = i + span;
    
    ind = [left_start:skip:left_end  right_start:skip:right_end];    
    
    ind = train(:, ind);
    ind = [ind(:,1)  ind(:, 2:m)-ones(size(ind(:, 2:m)))];    
    ind = ind*P;
    
    %compute empirical estimate of probabilities    
    count = histc(ind, unique(ind));
    tens_inv(unique(ind)) = count;
    
    %scale the computed multdim array and convert in to tensor class
    scaled_tens = scaled_tens + tens/L;
    scaled_tens_inv = scaled_tens_inv + tens_inv/L;
end


scaled_tens = scaled_tens/length(iter_ind);
scaled_tens_inv = scaled_tens_inv/length(iter_ind);


%scale the computed multdim array and convert in to tensor class
scaled_tensor = sptensor(scaled_tens);
scaled_tensor_inv = tensor(scaled_tens_inv);

%tensor's modes - 1:numObs,1:numObs
%eliminate first set of numObs modes

%modes to be eliminated are put as columns

mat_tensor = tenmat(scaled_tensor_inv, 1:numObs, 't');

%do svd by exracting only the K largest values/vectors
[U S V] = svds(sparse(mat_tensor.data), durMax*stateDim);

%drop singular values with small values
num = nnz(diag(S) >= 1e-30);
assert(num >= 1);

U = U(:,1:num);
S = S(1:num, 1:num);
V = V(:,1:num);

%U corresponds to second set of numObs modes
%V corresponds to first set of numObs modes

%So, we need USV' * VS^U' = UU'

mat_tensor(:) = (V*diag(1./diag(S))*U')';

res = ttt(scaled_tensor, sptensor(tensor(mat_tensor)), 1:numObs);

durTensor = res;




