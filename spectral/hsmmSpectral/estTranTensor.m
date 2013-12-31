function tranTensor = estTranTensor(train,...
                                    numObs, ...
                                    span, ...
                                    obsDim, ...
                                    stateDim, ...
                                    durMin, ...
                                    durMax, A1_true, A_true, D_true, O_true, flg)

%sequences size
[L, N] = size(train);

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


%find start and end indeces in the data sequence, so that tensors can be
%estimated from enough data
tran_start_ind = span + 1;
tran_end_ind = N - span;

iter_ind = tran_start_ind : tran_end_ind;


for k = 1:length(iter_ind)
    
    i = iter_ind(k);
    
    tens = zeros(obsDim*(ones(1, 2*numObs+1)));
    tens_inv = zeros(obsDim*(ones(1, 2*numObs)));
        
    %compute indeces    
    ind = [getLeftInd(i-1, stateDim, durMax, numObs)  i  getRightInd(i+1, stateDim, durMax, numObs)];
    
    ind = train(:, ind);
    
    %for sequences of unequal length
    %rows which contain -1 are invalid, need to delete whole row
    inv_ind = (sum(ind<0,2)>0);
    ind(inv_ind,:)=[];    
    
    ind = [ind(:,1)  ind(:, 2:m)-ones(size(ind(:, 2:m)))];
    ind = ind*P;
        
    %compute empirical estimate of probabilities
    count = histc(ind, unique(ind));
    tens(unique(ind)) = count/L;
    
    
    %%%
    if flg
        Le = getLeftCondProb(O_true, D_true, A_true, A1_true, [getLeftInd(i-1, stateDim, durMax, numObs) i], 0);
        Ri = getRightCondProb(O_true, D_true, A_true, [i i getRightInd(i+1, stateDim, durMax, numObs)], 0);
        tens = ttt(Le, Ri, [length(Le.size)-1 length(Le.size)], [length(Ri.size)-1 length(Ri.size)]);
    end
    %%%
    
    
        
    %compute indeces for inverse tensor
    ind = [getLeftInd(i-1, stateDim, durMax, numObs)  getRightInd(i+1, stateDim, durMax, numObs)];
    
    ind = train(:, ind);
    
    %for sequences of unequal length
    %rows which contain -1 are invalid, need to delete whole row
    inv_ind = (sum(ind<0,2)>0);
    ind(inv_ind,:)=[];
        
    ind = [ind(:,1)  ind(:, 2:m_inv)-ones(size(ind(:, 2:m_inv)))];
    ind = ind*P_inv;
    
    %compute empirical estimate of probabilities
    count = histc(ind, unique(ind));
    tens_inv(unique(ind)) = count/L;
    
    
    %%%
    if flg
        Le = getLeftCondProb(O_true, D_true, A_true, A1_true, [getLeftInd(i-1, stateDim, durMax, numObs) i], 0);
        Ri = getRightCondProb(O_true, D_true, A_true, [i getRightInd(i+1, stateDim, durMax, numObs)], 0);
        tens_inv = ttt(Le, Ri, [length(Le.size)-1 length(Le.size)], [length(Ri.size)-1 length(Ri.size)]);
    end
    %%%
    
    
    %scale the computed multdim array and convert in to tensor class
    scaled_tens = scaled_tens + tens;
    scaled_tens_inv = scaled_tens_inv + tens_inv; 
end

scaled_tens = scaled_tens/length(iter_ind);
scaled_tens_inv = scaled_tens_inv/length(iter_ind);

% N = length(iter_ind);
% scaled_tens = scaled_tens/((N^2+N)/2);
% scaled_tens_inv = scaled_tens_inv/((N^2+N)/2);

%scale the computed multdim array and convert in to tensor class
scaled_tensor = tensor(scaled_tens);
scaled_tensor_inv = tensor(scaled_tens_inv);

% scaled_tensor = sptensor(computeO1O2O4O5O6(A1_true, A_true, D_true(:,:,1), D_true, O_true));
% scaled_tensor_inv = tensor(computeO1O2O5O6(A1_true, A_true, D_true(:,:,1), D_true, O_true));

%tensor_inv modes - 1:numObs, 1:numObs
%eliminate first set of numObs modes

%modes to be eliminated are put as columns

mat_tensor = tenmat(scaled_tensor_inv, 1:numObs, 't');


%do svd by exracting only the K largest values/vectors
[U S V] = svd(mat_tensor.data);

%drop singular values with small values and try to keep at least durMax*stateDim
num = min(nnz(diag(S) >= eps), durMax*stateDim);
% num = 30;
assert(num >= 1);

U = U(:,1:num);
S = S(1:num, 1:num);
V = V(:,1:num);

%U corresponds to second set of numObs modes
%V corresponds to first set of numObs modes

%So, we need USV' * VS^U' = UU'

mat_tensor(:) = (V*diag(1./diag(S))*U')';

tranTensor = tensor(ttt(scaled_tensor, tensor(mat_tensor), 1:numObs));



% indM = getRightInd(4, stateDim, durMax, numObs);
% M = getRightCondProb(O_true, D_true, A_true, [2 indM], 0); %O..O|xt,dt
% 
% Mm = tenmat(M, 1:length(indM));
% [U S V] = svd(Mm.data, 0);
% 
% num = min(nnz(diag(S) >= eps), durMax*stateDim);
% assert(num == durMax*stateDim);
% 
% invMm = (V*diag(1./diag(S))*U');
% %invMm = pinv(Mm.data);
% 
% Mm(:) = invMm';
% M = tensor(Mm);
% 
% eA = tensor(embedA(A_true)); %x_t x_t x_t-1 d_t-1 d_t-1
% 
% indK = getRightInd(4, stateDim, durMax, numObs);
% K = getRightCondProb(O_true, D_true, A_true, [3 indK], 1); %O..O|xt,dt-1
% 
% res2 = ttt(K, eA, [length(K.size)-1 length(K.size)], [2 4]);
% res2 = ttt(res2, M, [length(indK)+2 length(indK)+3], [length(M.size)-1 length(M.size)]);
% 
% tranTensor = ttt(tensor(O_true), res2, 2, length(indK)+1);






