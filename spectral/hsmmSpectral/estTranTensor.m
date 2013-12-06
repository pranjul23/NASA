function tranTensor = estTranTensor(train,...
                                    numObs, ...
                                    span, ...
                                    obsDim, ...
                                    stateDim, ...
                                    durMin, ...
                                    durMax, A1_true, A_true, D_true, O_true)

%sequences size
[L, N] = size(train);


%find start and end indeces in the data sequence, so that tensors can be
%estimated from enough data
tran_start_ind = span + 2; %+2 because clique dixixi+1 spans 2 observations and Oi+1 is the one we reference against
tran_end_ind = N - span;

iter_ind = tran_start_ind : tran_end_ind;


%compute linear indeces for efficient matrix access
m = 2*numObs+1;
P = ones(m,1);
for k = 2:m
    P(k) = P(k-1)*obsDim;
end


%compute linear indeces for efficient matrix access
m_inv = 2*numObs;
P_inv = ones(m_inv,1);
for k = 2:m_inv
    P_inv(k) = P_inv(k-1)*obsDim;
end

scaled_tens = zeros(obsDim*(ones(1, 2*numObs+1)));
scaled_tens_inv = zeros(obsDim*(ones(1, 2*numObs)));

for k = 1:length(iter_ind)
    
    i = iter_ind(k);
    
    tens = zeros(obsDim*(ones(1, 2*numObs+1)));
    tens_inv = zeros(obsDim*(ones(1, 2*numObs)));
        
    %compute indeces    
    ind = [getLeftInd(i-2, stateDim, durMax, numObs) i getRightInd(i+1, stateDim, durMax, numObs) ];
    
    ind = train(:, ind);
    ind = [ind(:,1)  ind(:, 2:m)-ones(size(ind(:, 2:m)))];
    ind = ind*P;
        
    %compute empirical estimate of probabilities
    count = histc(ind, unique(ind));
    tens(unique(ind)) = count;
    
        
    %compute indeces for inverse tensor
    ind = [getLeftInd(i-2, stateDim, durMax, numObs)  getRightInd(i+1, stateDim, durMax, numObs)];
    
    ind = train(:, ind);
    ind = [ind(:,1)  ind(:, 2:m_inv)-ones(size(ind(:, 2:m_inv)))];
    ind = ind*P_inv;
    
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

% scaled_tensor = sptensor(computeO1O2O4O5O6(A1_true, A_true, D_true(:,:,1), D_true, O_true));
% scaled_tensor_inv = tensor(computeO1O2O5O6(A1_true, A_true, D_true(:,:,1), D_true, O_true));

%tensor_inv modes - 1:numObs, 1:numObs
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
%mat_tensor(:) = (V*inv(U'*mat_tensor.data*V)*U')';

%tranTensors{c} = ttt(scaled_tensor, tensor(mat_tensor), 1:numObs);
res = tensor(ttt(scaled_tensor, sptensor(tensor(mat_tensor)), 1:numObs));

tranTensor = res;




indM = getRightInd(4, stateDim, durMax, numObs);
M = getCondProb(O_true, D_true, A_true, [2 indM], 0); %O..O|xt,dt

Mm = tenmat(M, 1:length(indM));
[U S V] = svds(Mm.data, durMax*stateDim);
invMm = (V*diag(1./diag(S))*U');
%invMm = pinv(Mm.data);

Mm(:) = invMm';
M = tensor(Mm);

eA = tensor(embedA(A_true)); %x_t x_t x_t-1 d_t-1 d_t-1

indK = getRightInd(4, stateDim, durMax, numObs);
K = getCondProb(O_true, D_true, A_true, [3 indK], 1); %O..O|xt,dt-1

res2 = ttt(K, eA, [length(K.size)-1 length(K.size)], [2 4]);
res2 = ttt(res2, M, [length(indK)+2 length(indK)+3], [length(M.size)-1 length(M.size)]);

% tranTensor = ttt(tensor(O_true), res2, 2, length(indK)+1);












