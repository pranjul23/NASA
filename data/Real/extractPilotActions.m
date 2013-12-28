%function to extract pilot actions from NASA data
function res = extractPilotActions(data)


N = length(data);
tmp = cell(N,1);
res = cell(N,1);

plotting = 0;

CH = zeros(13, N);

for i=1:N
    
    load(data{i});
    
    assert(PH.Rate == 1);
    
    % select specific part of flight
    %-----------------------------------------------------------------------
    % select whole flight
    ind = ((PH.data ~= 2)&(PH.data ~= 1)&(PH.data ~= 0));
    
    %find index of last time stamp
    d = abs(diff(ind));
    endPoint = max(find(d==max(d)));
    
    %find index of "initial descend" time stamp
    ind = PH.data == 6;
    d = abs(diff(ind));
    
    mi = min(find(d==max(d)));
    ma = max(find(d==max(d)));
    
    startPoint = mi + round((ma-mi)/2);
    
    %make sure this portion no longer than 15 minutes
    if (endPoint-startPoint) > 1000
        startPoint = endPoint - 1000;
    end
    
    ind = zeros(size(PH.data))';
    ind(startPoint:endPoint)=1;
    ind = logical(ind);
    %-----------------------------------------------------------------------
    
    
    L = nnz(ind)-1;
    
    if plotting
        subplot(14,1,1);
        plot(PH.data(ind));axis off;        
    end
    
    actions = zeros(13, L);
    
    % airbrake position
    assert(ABRK.Rate == 1);
    D = ABRK.data(ind);
    actions(1, :) = double( abs(diff(D)) > 0 )';
        
    
    if plotting
        subplot(14,1,2);
        plot(D);
        axis off;
        fprintf('%.1f ',nnz(actions(1,:))/L*100);        
    end    
    CH(1, i) = nnz(actions(1,:))/L*100;
    
    % selected altitude
    assert(ALTS.Rate == 1);
    D = ALTS.data(ind);
    actions(2, :) = double( abs(diff(D)) > 0 )';
    
    if plotting
        subplot(14,1,3);
        plot(D);axis off;
        fprintf('%.1f ',nnz(actions(2,:))/L*100);
    end
    CH(2, i) = nnz(actions(2,:))/L*100;
    
    % selected airspeed
    assert(CASS.Rate == 1);
    D = CASS.data(ind);
    actions(3, :) = double( abs(diff(D)) > 0 )';
    
    if plotting
        subplot(17,1,4);
        plot(D);axis off;
        fprintf('%.1f ',nnz(actions(3,:))/L*100);
    end
    CH(3, i) = nnz(actions(3,:))/L*100;
    
    
    % control column position capt
    assert(CCPC.Rate == 2);
    D = CCPC.data;
    D = D(1:2:end);
    D = D(ind);
    actions(4, :) = double( abs(diff(D)) > 0 )';
    
    if plotting
        subplot(14,1,5);
        plot(D);axis off;
        fprintf('%.1f ',nnz(actions(4,:))/L*100);
    end
    CH(4, i) = nnz(actions(4,:))/L*100;

    
    % control wheel position capt
    assert(CWPC.Rate == 2);
    D = CWPC.data;
    D = D(1:2:end);
    D = D(ind);
    actions(5, :) = double( abs(diff(D)) > 0 )';
    
    
    if plotting
        subplot(14,1,6);
        plot(D);axis off;
        fprintf('%.1f ',nnz(actions(5,:))/L*100);
    end
    CH(5, i) = nnz(actions(5,:))/L*100;
    
    
    % flap position
    assert(FLAP.Rate == 1);
    D = FLAP.data(ind);
    actions(6, :) = double( abs(diff(D)) > 0 )';
    
    if plotting
        subplot(14,1,7);
        plot(D);axis off;
        fprintf('%.1f ',nnz(actions(6,:))/L*100);
    end
    CH(6, i) = nnz(actions(6,:))/L*100;
    
    % selected heading
    assert(HDGS.Rate == 1);
    D = HDGS.data(ind);
    actions(7, :) = double( abs(diff(D)) > 0 )';
    
    if plotting
        subplot(14,1,8);
        plot(D);axis off;
        fprintf('%.1f ',nnz(actions(7,:))/L*100);
    end
    CH(7, i) = nnz(actions(7,:))/L*100;
    
    % selected mach
    assert(MNS.Rate == 1);
    D = MNS.data(ind);
    actions(8, :) = double( abs(diff(D)) > 0 )';
    
    if plotting
        subplot(14,1,9);
        plot(D);axis off;
        fprintf('%.1f ',nnz(actions(8,:))/L*100);
    end
    CH(8, i) = nnz(actions(8,:))/L*100;
    
    % N1 command
    assert(N1C.Rate == 4);
    D = N1C.data;
    D = D(1:4:end);
    D = D(ind);
    actions(9, :) = double( abs(diff(D)) > 0 )';
    
    if plotting
        subplot(14,1,10);
        plot(D);axis off;
        fprintf('%.1f ',nnz(actions(9,:))/L*100);
    end
    CH(9, i) = nnz(actions(9,:))/L*100;
    
    % rudder pedal position
    assert(RUDP.Rate == 2);
    D = RUDP.data;
    D = D(1:2:end);
    D = D(ind);
    actions(10, :) = double( abs(diff(D)) > 0 )';
    
    if plotting
        subplot(14,1,11);
        plot(D);axis off;
        fprintf('%.1f ',nnz(actions(10,:))/L*100);
    end
    CH(10, i) = nnz(actions(10,:))/L*100;
    
    % selected vertical speed
    assert(VSPS.Rate == 1);
    D = VSPS.data(ind);
    actions(11, :) = double( abs(diff(D)) > 0 )';
    
    if plotting
        subplot(14,1,12);
        plot(D);axis off;
        fprintf('%.1f ',nnz(actions(11,:))/L*100);
    end
    CH(11, i) = nnz(actions(11,:))/L*100;
    
    % thrust mode
    assert(TMODE.Rate == 1);
    D = TMODE.data(ind);
    actions(12, :) = double( abs(diff(D)) > 0 )';
    
    if plotting
        subplot(14,1,13);
        plot(D);axis off;
        fprintf('%.1f ',nnz(actions(12,:))/L*100);
    end
    CH(12, i) = nnz(actions(12,:))/L*100;
    
    % vertical engage mode
    assert(VMODE.Rate == 1);
    D = VMODE.data(ind);
    actions(13, :) = double( abs(diff(D)) > 0 )';
    
    if plotting
        subplot(14,1,14);
        plot(D); axis off;
        fprintf('%.1f \n',nnz(actions(13,:))/L*100);
    end
    CH(13, i) = nnz(actions(13,:))/L*100;
    
%     pause(1);
    
    tmp{i} = actions;
    
    %flatten out the data
    R = [1:13]';    
    actions = actions.*(R*ones(1,L));
    
    %random permutations at each time stamp    
    for j=1:L
        actions(:,j) = actions(randperm(13)',j);
    end
    actions = reshape(actions, numel(actions), 1)';
    actions = actions(actions ~= 0);
    
    res{i} = actions;
end








