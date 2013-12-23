function result = testSpectModel2(test, obsDim, numObs, ...
                                  rootTensor,  ...
                                  tailTensor,  ...
                                  obsTensor,  ...
                                  tranTensor, ...
                                  durTensor)
                                   
                                   
                                   
[L, N] = size(test);

result = zeros(L,1);


% rootTensor.tensor = tensor(rootTensor.tensor);
% tailTensor = tensor(tailTensor);
% obsTensor = tensor(obsTensor);
% tranTensor = tensor(tranTensor);
% durTensor = tensor(durTensor);

for i=1:L
    sequence = test(i,:);
    
    last_ind = length(sequence);
    val = sequence(last_ind);
    
    delta = zeros(obsDim,1);
    delta(val) = 1;    
    tdelta = tensor(delta);
                
    tail_obs = squeeze(ttt(obsTensor, tdelta, 1, 1));            
    
    res = ttt(tailTensor, tail_obs, 1, 1);
    
    scale = 1;%10^(round(length(sequence)/2));    
    res = res * scale;
    
    %3 is the last index  before rootTensor starts
    for k = last_ind-1 : -1 : 3
        
        res = ttt(durTensor, res, 1:numObs, 1:numObs);        
        
        val = sequence(k);
        delta = zeros(obsDim,1);
        delta(val) = 1;    
        tdelta = tensor(delta);    
            
        obs = squeeze(ttt(obsTensor, tdelta, 1, 1));
        
        leaf = ttt(tranTensor, obs, 1, 1);
        
        res = ttt(leaf, res, 1:numObs, 1:numObs);        
    end
    
    
    res = ttt(durTensor, res, 1:numObs, 1:numObs);
    
    val = sequence(1);
    delta = zeros(obsDim,1);
    delta(val) = 1;    
    tdelta = tensor(delta);
        
    obs1 = squeeze(ttt(obsTensor, tdelta, 1, 1));
    
    
    val = sequence(2);
    delta = zeros(obsDim,1);
    delta(val) = 1;    
    tdelta = tensor(delta);
    
    obs2 = squeeze(ttt(obsTensor, tdelta, 1, 1));
        
    
    root = ttt(rootTensor, obs1, 1, 1);
    
    root = ttt(root, obs2, 1, 1);    
    
    res = ttt(root, res, 1:numObs, 1:numObs); 
    result(i) = res;
    i
end

