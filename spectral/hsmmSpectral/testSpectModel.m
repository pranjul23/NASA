function result = testSpectModel(test, obsDim, stop_ind, numObs, ...
                                 rootTensor,  ...
                                 tailTensor,  ...
                                 obsTensors,  ...
                                 tranTensors, ...
                                 durTensors)
                                   
                                   
                                   
[L, N] = size(test);

result = zeros(L,1);

for i=1:L
    sequence = test(i,:);
    
    last_ind = length(sequence);
    val = sequence(last_ind);
    
    delta = zeros(obsDim,1);
    delta(val) = 1;    
    tdelta = sptensor(delta);
        
    mode = find(tailTensor.var_ind == last_ind);
    
    %select the slice corresponding to actual observation (tensor)
    assert(tailTensor.seq_ind == last_ind);    
    tail_tran = squeeze(ttt(tailTensor.tensor, tdelta, mode, 1));
        
    %select the slice corresponding to actual observation (scalar)
    tail_obs = obsTensors(last_ind).tensor(val, val);
    
    %SCALE: myltiply by some big number to avoid underflow
    scale = 1;%10^(round(length(sequence)/2));
    
    res = tail_tran * tail_obs * scale;
    
    %the indeces of the variables in the result tensor
    res_ind = tailTensor.var_ind;
    res_ind(mode) = [];
        
    %assert(all(durTensors(last_ind).ind(numObs+1:2*numObs) == 1+1:numObs));    
        
    
    for k = last_ind-1 : -1 : stop_ind
        
        assert(durTensors(k+1).var_ind(1:numObs) == res_ind(1:numObs));        
        res = ttt(durTensors(k+1).tensor, res, 1:numObs, 1:numObs);        
        %res_ind remains the same after the previous step
        
        val = sequence(k);
        delta = zeros(obsDim,1);
        delta(val) = 1;    
        tdelta = sptensor(delta);    
    
        mode = find(tranTensors(k).var_ind == k);
        tran = squeeze(ttt(tranTensors(k).tensor, tdelta, mode, 1));
        
        leaf_ind = tranTensors(k).var_ind;
        leaf_ind(mode) = [];
        
        obs = obsTensors(k).tensor(val, val);
    
        leaf = tran * obs;
    
        assert(leaf_ind(1:numObs) == res_ind(1:numObs));
        res = ttt(leaf, res, 1:numObs, 1:numObs);
        
        res_ind = leaf_ind(numObs+1:2*numObs);
    end
    
    
    for k = stop_ind-1: -1 : 3
        
        res = ttt(durTensors(stop_ind).tensor, res, 1:numObs, 1:numObs);
    
        val = sequence(k);
        delta = zeros(obsDim,1);
        delta(val) = 1;    
        tdelta = sptensor(delta);
        
        mode = find(tranTensors(stop_ind).var_ind == stop_ind);
        tran = squeeze(ttt(tranTensors(stop_ind).tensor, tdelta, mode, 1));
    
        obs = obsTensors(k).tensor(val, val);
    
        leaf = tran * obs;
    
        res = ttt(leaf, res, numObs+1:2*numObs, 1:numObs);        
    end
    
    
    res = ttt(durTensors(stop_ind).tensor, res, 1:numObs, 1:numObs);
    
    val = sequence(1);
    delta = zeros(obsDim,1);
    delta(val) = 1;    
    tdelta = sptensor(delta);
    
    mode = find(rootTensor.var_ind == 1);
    root = ttt(rootTensor.tensor, tdelta, mode, 1);
    
    obs1 = obsTensors(1).tensor(val, val);
    
    
    val = sequence(2);
    delta = zeros(obsDim,1);
    delta(val) = 1;    
    tdelta = sptensor(delta);
    
    mode = find(rootTensor.var_ind == 1);
    root = squeeze(ttt(root, tdelta, mode, 1));
    
    obs2 = obsTensors(2).tensor(val, val);
    
    root = root * obs1;
    root = root * obs2;
    
    res = ttt(root, res, 1:numObs, 1:numObs); 
    result(i) = res;
end


