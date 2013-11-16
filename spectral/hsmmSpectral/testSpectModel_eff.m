function result = testSpectModel_eff(test, obsDim, stop_ind, numObs, ...
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
                
    tail_obs = squeeze(ttt(obsTensors(last_ind).tensor, tdelta, 1, 1));            
    
    res = ttt(tailTensor.tensor, tail_obs, 1, 1);
    
    scale = 1;%10^(round(length(sequence)/2));    
    res = res * scale;
    
    for k = last_ind-1 : -1 : 3
        
        res = ttt(durTensors(last_ind-1).tensor, res, 1:numObs, 1:numObs);        
        
        val = sequence(k);
        delta = zeros(obsDim,1);
        delta(val) = 1;    
        tdelta = sptensor(delta);    
            
        obs = squeeze(ttt(obsTensors(last_ind).tensor, tdelta, 1, 1));
        
        leaf = ttt(tranTensors(last_ind-1).tensor, obs, 1, 1);
        
        res = ttt(leaf, res, 1:numObs, 1:numObs);        
    end
    
    
    res = ttt(durTensors(last_ind-1).tensor, res, 1:numObs, 1:numObs);
    
    val = sequence(1);
    delta = zeros(obsDim,1);
    delta(val) = 1;    
    tdelta = sptensor(delta);
        
    obs1 = squeeze(ttt(obsTensors(last_ind).tensor, tdelta, 1, 1));
    
    
    val = sequence(2);
    delta = zeros(obsDim,1);
    delta(val) = 1;    
    tdelta = sptensor(delta);
    
    obs2 = squeeze(ttt(obsTensors(last_ind).tensor, tdelta, 1, 1));
        
    
    root = ttt(rootTensor.tensor, obs1, 1, 1);
    
    root = ttt(root, obs2, 1, 1);    
    
    res = ttt(root, res, 1:numObs, 1:numObs); 
    result(i) = res;
end


