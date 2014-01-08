%function to select flights that have similar rate of change per variable
function ind = selectAnomFlights(rates)

[K, N] = size(rates);

ind = true(K, N);
edges = 0:100;


for i=1:K
    
    D = rates(i,:);
    
%     %find a mode of histogram
%     n = histc(D, edges);
%     
%     bar(edges, n, 'histc');
%     
%     mi = input('Specify min: ');
%     
%     if mi >= 0
%         
%         ma = input('Specify max: ');
%         
%         ind_tmp = D<=ma & D>=mi;
%         D_tmp = D(ind_tmp);
%         
%         %make sure data is not spread out too much
%         in = (D_tmp >= prctile(D_tmp,2)) & (D_tmp <= prctile(D_tmp,98));
%         
%         c = find(ind_tmp);
%         c = c(in);
%         
%         ind(i,:) = false;
%         ind(i,c) = true;
%         
%     else
        %make sure data is not spread out too much
        ind(i,:) = D > prctile(D,99);
%     end
end
