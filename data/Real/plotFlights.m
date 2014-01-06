%script to plot flight trajectories from NASA dataset

endPoints = zeros(15000,3);
fileNames = cell(15000,1);

origin = pwd;


folder = 'Tail_6';

c1 = 1;
c2 = 1;

files = 71:73;

%% ========================================================================
% determine take-off and landing positions    

for f = 1:length(files)
    
    k = 1;
    
    while(1)        
                
        loc = [folder num2str(files(f)) '_' num2str(k)];
        
        if ~exist(loc, 'dir')
           break; 
        end
        
        list = dir(loc);
        cd(loc);
        
        for i=3:length(list)
            
            load(list(i).name);
            
            if(length(unique(PH.data)) < 6)
                fprintf('skipped %d: corrupt PH - missing symbols ', i);
                disp(unique(PH.data)');
                continue;
            end
            
            if(any(unique(PH.data)>7))
                fprintf('skipped %d: corrupt PH - extra symbols ', i);
                disp(unique(PH.data)');
                continue;
            end
            
            if( length(unique(LATP.data)) < 10 || length(unique(LONP.data)) < 10 )
                fprintf('skipped %d: corrupt LATP or LONP\n', i);
                continue;
            end
            
            
            %select part of flight from takeoff to landing
            assert(PH.Rate == 1);
            ind = ((PH.data ~= 2)&(PH.data ~= 1)&(PH.data ~= 0));
            
            assert(length(unique(PH.data(ind))) <= 5);
            
            
            %altitude
            assert(RALT.Rate == 8);
            A = RALT.data(1:8:end);
            
            if(max(A)<2000)
                fprintf('skipped %d: currupt RALT ', i);
                disp(max(A));
                continue;
            end
            
            
            %lat lon alt [deg deg m]
            assert(LONP.Rate == 1); assert(LATP.Rate == 1);
            LLA = [LATP.data(ind) LONP.data(ind) 0.3048*A(ind)];
            
            %transform to ECEF coordinate system
            xyz = lla2ecef(LLA);
            
            if any(any(xyz == 0))
                fprintf('skipped %d: corrupt trajectory \n', i);
                continue;
            end
            
            if norm(xyz(1,:)-xyz(end,:)) < 10000
                fprintf('skipped %d: too short trajectory \n', i);
                continue;
            end
            
            plot3(xyz(:,1), xyz(:,2), xyz(:,3), '-c');
            hold on
            plot3(xyz(end,1), xyz(end,2), xyz(end,3), 'o', 'markersize', 5, 'MarkerFaceColor','b','MarkerEdgeColor','b');
            hold on
            plot3(xyz(1,1), xyz(1,2), xyz(1,3), 'o', 'markersize', 5, 'MarkerFaceColor','r','MarkerEdgeColor','r');
            
            endPoints(c1,:) = xyz(1,:); %odd - departure airport
            endPoints(c1+1,:) = xyz(end,:); %even - arrival airport
            
            fileNames{c2} = [loc '/' list(i).name];
            
            c1 = c1 + 2;
            c2 = c2 + 1;
        end
        
        %  grid on
        %  axis equal
        
        cd(origin);
        
        k = k + 1;
    end    
end

% %quick visualization of locations of all the airports
% plot3(endPoints(1:2:end,1), endPoints(1:2:end,2), endPoints(1:2:end,3), '.r'); %departure
% hold on;
% plot3(endPoints(2:2:end,1), endPoints(2:2:end,2), endPoints(2:2:end,3), '.b'); %arrival
% grid on
% axis equal

%% ========================================================================
% cluster

%first, remove zero rows (extra left after preallocaiton of space)
endPoints(c1:end, :) = [];
fileNames(c2:end) = [];


%compute all pairwise distances
D = dist(endPoints');

%points that are close are from same airport
L = D<=50000;

%find # of clusters and cluster membership
[S, C] = graphconncomp(sparse(double(L)), 'Directed', false);

[~, largest_clust_ID] = sort((histc(C, unique(C))), 'descend');

%% ========================================================================
% select flights 

% SINGLE ARR_DEP PAIR
% %cluster ID for departure airport
% dep = 88;
% 
% %select odd indeces - they correspond to departure airports
% ind = find(C==dep);
% 
% ind_odd_dep = ind(logical(mod(ind,2)));
% ind_even = ind_odd_dep + 1;
% 
% %find arrival airport such that number of depart-arrival is maximized
% arriv = C(ind_even);
% 
% arriv_uniq = unique(arriv);
% M = -inf;
% M_val = 0;
% 
% for i=1:length(arriv_uniq)
%     tmp = nnz(arriv == arriv_uniq(i));
%     if tmp > M
%         M = tmp;
%         M_val = arriv_uniq(i);
%     end
% end
% 
% ind = find(C==M_val);
% ind_odd_arr = ind(~logical(mod(ind,2)))-1;
% 
% ind = intersect(ind_odd_dep, ind_odd_arr);
% 
% 
% %transform to indeces for fileNames
% ind = (ind+1)/2;


% ALL FLIGHTS COMMING TO ARR
arr = largest_clust_ID(1);

%select all aiplanes arriving to the same airport
ind = find(C==arr);

ind_even_arr = ind(~logical(mod(ind,2)));

% %transform to indeces for fileNames
ind = ind_even_arr/2;



%pick selected file names
data = fileNames(ind);

%% ========================================================================
%visualize resutls

% N = length(data);
% for i=1:N
%     
%     load(data{i});
%     
%     %select part of flight from takeoff to landing
%     ind = ((PH.data ~= 2)&(PH.data ~= 1)&(PH.data ~= 0));
%     
%     %altitude
%     A = RALT.data(1:8:end);
%     
%     %lat lon alt [deg deg m]
%     assert(LONP.Rate == 1); assert(LATP.Rate == 1);
%     LLA = [LATP.data(ind) LONP.data(ind) 0.3048*A(ind)];
%     
%     %transform to ECEF coordinate system
%     xyz = lla2ecef(LLA);
%     
%     plot3(xyz(:,1), xyz(:,2), xyz(:,3), '-g');
%     hold on
%     plot3(xyz(end,1), xyz(end,2), xyz(end,3), 'o', 'markersize', 5, 'MarkerFaceColor','b','MarkerEdgeColor','b');
%     
%     hold on
%     plot3(xyz(1,1), xyz(1,2), xyz(1,3), 'o', 'markersize', 5, 'MarkerFaceColor','r','MarkerEdgeColor','r');
% end
% 
% grid on
% axis equal

%% ========================================================================
% exctract pilot actions

[actions rates] = extractPilotActions(data);

save('../../spectral/hsmmSpectral/data.mat', 'actions', 'rates');


















