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
                
    %select the slice corresponding to actual observation (tensor)
    mode = find(obsTensors(last_ind).var_ind == last_ind);
    tail_obs = squeeze(ttt(obsTensors(last_ind).tensor, tdelta, mode, 1));            
    
    assert(tailTensor.seq_ind == last_ind);    
    
    mode = find(tailTensor.var_ind == last_ind);
    res = ttt(tailTensor.tensor, tail_obs, mode, 1);
    
    %SCALE: myltiply by some big number to avoid underflow
    scale = 1;%10^(round(length(sequence)/2));    
    res = res * scale;
    
    %the indeces of the variables in the result tensor
    res_ind = tailTensor.var_ind;
    res_ind(mode) = [];
        
    %assert(all(durTensors(last_ind).ind(numObs+1:2*numObs) == 1+1:numObs));    
        
    
    for k = last_ind-1 : -1 : stop_ind
        
        assert(all(durTensors(k).var_ind(1:numObs) == res_ind(1:numObs)));        
        res = ttt(durTensors(k).tensor, res, 1:numObs, 1:numObs);        
        
        res_ind = durTensors(k).var_ind(numObs+1:2*numObs);
        
        val = sequence(k);
        delta = zeros(obsDim,1);
        delta(val) = 1;    
        tdelta = sptensor(delta);    
            
        mode = find(obsTensors(k).var_ind == k);
        obs = squeeze(ttt(obsTensors(k).tensor, tdelta, mode, 1));
        
        mode = find(tranTensors(k).var_ind == k);
        leaf = ttt(tranTensors(k).tensor, obs, mode, 1);
        
        leaf_ind = tranTensors(k).var_ind;
        leaf_ind(mode) = [];
        
        assert(all(leaf_ind(1:numObs) == res_ind(1:numObs)));
        res = ttt(leaf, res, 1:numObs, 1:numObs);
        
        res_ind = leaf_ind(numObs+1:2*numObs);
    end
    
    
    for k = stop_ind-1: -1 : 3
        
        res = ttt(durTensors(stop_ind-1).tensor, res, 1:numObs, 1:numObs);
    
        val = sequence(k);
        delta = zeros(obsDim,1);
        delta(val) = 1;    
        tdelta = sptensor(delta);
        
        mode = find(obsTensors(k).var_ind == k);
        obs = squeeze(ttt(obsTensors(k).tensor, tdelta, mode, 1));
            
        mode = find(tranTensors(stop_ind).var_ind == stop_ind);
        leaf = ttt(tranTensors(stop_ind).tensor, obs, mode, 1);
                
        res = ttt(leaf, res, 1:numObs, 1:numObs);        
    end
    
    
    res = ttt(durTensors(stop_ind-1).tensor, res, 1:numObs, 1:numObs);
    
    val = sequence(1);
    delta = zeros(obsDim,1);
    delta(val) = 1;    
    tdelta = sptensor(delta);
        
    mode = find(obsTensors(1).var_ind == 1);
    obs1 = squeeze(ttt(obsTensors(1).tensor, tdelta, mode, 1));
    
    res_ind = rootTensor.var_ind;
    res_ind(mode) = [];
        
        
    
    val = sequence(2);
    delta = zeros(obsDim,1);
    delta(val) = 1;    
    tdelta = sptensor(delta);
    
    mode = find(obsTensors(2).var_ind == 2);
    obs2 = squeeze(ttt(obsTensors(2).tensor, tdelta, mode, 1));
        
    
    root = ttt(rootTensor.tensor, obs1, 1, 1);
    
    root = ttt(root, obs2, 1, 1);    
    
    res = ttt(root, res, 1:numObs, 1:numObs); 
    result(i) = res;
end


