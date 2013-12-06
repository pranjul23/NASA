function durTensor = estDurTensor(train, ...
                                  numObs,...
                                  span, ...
                                  obsDim, ...
                                  stateDim, ...
                                  durMin, ...
                                  durMax, A1_true, A_true, D_true, O_true)

%sequences size
[L, N] = size(train);

%duration clique should not go beyond transition clique (which is last)
dur_start_ind = span + 1;
dur_end_ind = N - span - 1;

iter_ind = dur_start_ind : dur_end_ind;

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
    ind = [getLeftInd(i-1, stateDim, durMax, numObs)  getRightInd(i+2, stateDim, durMax, numObs)];   
    
    ind = train(:, ind);
    ind = [ind(:,1)  ind(:, 2:m)-ones(size(ind(:, 2:m)))];    
    ind = ind*P;
    
    %compute empirical estimate of probabilities    
    count = histc(ind, unique(ind));
    tens(unique(ind)) = count;
            
    
    %compute indeces for inverse tensor
    ind = [getLeftInd(i-1, stateDim, durMax, numObs)  getRightInd(i+1, stateDim, durMax, numObs)];    
                
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

% scaled_tensor = sptensor(computeO1O2O5O6(A1_true, A_true, D_true(:,:,1), D_true, O_true));
% scaled_tensor_inv = sptensor(computeO1O2O4O5(A1_true, A_true, D_true(:,:,1), D_true, O_true));

%tensor's modes - 1:numObs,1:numObs
%eliminate first set of numObs modes

%modes to be eliminated are put as columns

mat_tensor = tenmat(scaled_tensor_inv, 1:numObs, 't');

%do svd by exracting only the K largest values/vectors
[U S V] = svds(sparse(mat_tensor.data), durMax*stateDim);

%drop singular values with small values
num = nnz(diag(S) >= 1e-7);
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




indM = getRightInd(4, stateDim, durMax, numObs);
M = getCondProb(O_true, D_true, A_true, [3 indM], 1); %O..O|xt,dt-1

Mm = tenmat(M, 1:length(indM));
[U S V] = svds(Mm.data, durMax*stateDim);
invMm = (V*diag(1./diag(S))*U');
%invMm = pinv(Mm.data);

Mm(:) = invMm';
M = tensor(Mm);

eD = tensor(embedD(D_true)); %d_t x_t x_t d_t-1

indK = getRightInd(4, stateDim, durMax, numObs);
K = getCondProb(O_true, D_true, A_true, [2 indK], 0); %O..O|xt,dt

res2 = ttt(M, eD, [length(M.size)-1 length(M.size)], [3 4]);
res2 = ttt(res2, K, [length(indM)+2 length(indM)+1], [length(K.size)-1 length(K.size)]);

res2m = tenmat(res2, 1:length(indM));
res2m(:) = res2m.data';

% durTensor = tensor(res2m);









