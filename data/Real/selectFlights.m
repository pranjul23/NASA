%function to select flights that have similar rate of change per variable
function ind_res = selectFlights(rates) 

[K, N] = size(rates);

ind = true(K, N);

for i=1:K
    
    D = rates(i,:); % ITS BETTER OT USE obj = gmdistribution.fit(D',2);
    
    edges = 0:100;
    step = 1;

    %max allowable difference between mean and mode
    max_delta = range(D)/4;
    
    %find a mode of histogram
    [n, bin] = histc(D, edges);    
    hist_mode = mode(bin);
    
    %make sure data is unimodal
    if abs( mean(D) -  hist_mode) > max_delta
        
        %keep increasing window around mode until diff 
        %between mode and mean is smaller than max_delta
        shift= step;
        while(1)
           mi = max([hist_mode-shift 0]);
           ma = min([hist_mode+shift 100]);
           
           if abs( mean(D(D<=ma & D>=mi)) -  hist_mode) > max_delta
               break;
           end
           
           shift = shift + step;           
        end
        
        ind_tmp = D<=ma & D>=mi;
        D_tmp = D(ind_tmp);
        
        %make sure data is not spread out too much
        %if std(D_tmp) > max_std
        in = (D_tmp >= prctile(D_tmp,2)) & (D_tmp <= prctile(D_tmp,98));
        %end
        
        c = find(ind_tmp);
        c = c(in);
        
        ind(i,:) = false;
        ind(i,c) = true;
        
    else
        
        %make sure data is not spread out too much
        %if std(D) > max_std
        ind(i,:) = (D >= prctile(D,5)) & (D <= prctile(D,95));
        %end
                
    end
end

%compute final indeces
ind = sum(ind,1);
ind_res = (ind == K);


