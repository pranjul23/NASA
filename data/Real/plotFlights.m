%script to plot flight trajectories from NASA dataset
close all
load '../../spectral/hsmmSpectral/dataNorm.mat'

N = length(data);

%visualize trajectories
f1 = figure;

for i=1:N
    
    load(data{i});
    
    %select part of flight from takeoff to landing
    ind = zeros(size(PH.data))';
    ind(startInd(i):endInd(i))=1;
    ind=logical(ind);
    
    %altitude
    A = RALT.data(1:8:end);
    
    %lat lon alt [deg deg m]
    assert(LONP.Rate == 1); assert(LATP.Rate == 1);
    LLA = [LATP.data(ind) LONP.data(ind) 0.3048*A(ind)];
    
    %transform to ECEF coordinate system
    xyz = lla2ecef(LLA);
    
    plot3(xyz(:,1), xyz(:,2), xyz(:,3), '-g');
    
    hold on
    plot3(xyz(end,1), xyz(end,2), xyz(end,3), 'o', 'markersize', 5, 'MarkerFaceColor','b','MarkerEdgeColor','b');
    
    hold on
    plot3(xyz(1,1), xyz(1,2), xyz(1,3), 'o', 'markersize', 5, 'MarkerFaceColor','c','MarkerEdgeColor','c');
end

grid on
axis equal





%number of actions
A = 10;

%pilot actions
act_names = cell(A,1);
k=1;

act_names{k} = 'ABRK'; %airbrake position
k=k+1;
act_names{k} = 'ALTS'; %selected altitude
k=k+1;
act_names{k} = 'CASS'; %selected airspeed
k=k+1;
% act_names{k} = 'CCPC'; %control column position capt
% k=k+1;
% act_names{k} = 'CWPC'; %control wheel position capt
% k=k+1;
act_names{k} = 'FLAP'; %flap position
k=k+1;
act_names{k} = 'HDGS'; %selected heading
k=k+1;
act_names{k} = 'MNS'; %selected mach
k=k+1;
act_names{k} = 'N1C'; %N1 command
k=k+1;
% act_names{k} = 'RUDP'; %rudder pedal position
% k=k+1;
act_names{k} = 'VSPS'; %selected vertical speed
k=k+1;
act_names{k} = 'TMODE'; %thrust mode
k=k+1;
act_names{k} = 'VMODE'; %vertical engage mode
k=k+1;

rates = [1 1 1 1 1 1 4 1 1 1];

%visualize actions
f2 = figure;

for i = 1:N
    
    load(data{i});        
        
    %select part of flight from takeoff to landing
    ind = zeros(size(PH.data))';
    ind(startInd(i):endInd(i))=1;
    ind=logical(ind);    
    
    for k = 1:A
        
        variab = genvarname(act_names{k});
       
        if rates(k) == 1
            D = eval([variab '.data(ind)']);
        else
            D = eval([variab '.data']);
            D = D(1:rates(k):end);
            D = D(ind);
        end
        
        subplot(A+1,1,k);
        plot(D);
        set(gca,'XTickLabel',[]);
        ylabel(act_names{k});
    end    
    
    subplot(A+1,1,k+1);
    plot(PH.data(ind));
    set(gca,'XTickLabel',[])
    ylabel('PH');
    
    
    AA = RALT.data(1:8:end);
    
    %lat lon alt [deg deg m]
    LLA = [LATP.data(ind) LONP.data(ind) 0.3048*AA(ind)];
    
    %transform to ECEF coordinate system
    xyz = lla2ecef(LLA);
    
    figure(f1);
    o1 = plot3(xyz(:,1), xyz(:,2), xyz(:,3), '-r', 'linewidth', 3);
    o2 = plot3(xyz(1,1), xyz(1,2), xyz(1,3), 'o', 'markersize', 10, 'MarkerFaceColor','k','MarkerEdgeColor','k');
    
    pause(5);    
    
    delete(o1);
    delete(o2);
    figure(f2);    
end










