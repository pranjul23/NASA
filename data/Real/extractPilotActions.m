%function to extract pilot actions from NASA data
function [res CH data startInd endInd ind_bad] = extractPilotActions(data)


N = length(data);
tmp = cell(N,1);
res = cell(N,1);

plotting = 0;

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

%sampling rates of the above actions
% rates = [1 1 1 2 2 1 1 1 4 2 1 1 1]; %original
rates = [1 1 1 1 1 1 4 1 1 1];

CH = zeros(A, N);

startInd = zeros(N, 1);
endInd = zeros(N, 1);

counter = 1;
flag = 0;

%indeces to use to extract sequences from data which are used
ind_to_keep = [];

for i=1:N
    
    load(data{i});
    
    assert(PH.Rate == 1);
    
    % select specific part of flight
    %----------------------------------------------------------------------
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
    
    %make sure this portion is no longer than 15 minutes
    if (endPoint-startPoint) > 1000
        startPoint = endPoint - 1000;
    end
    
    ind = zeros(size(PH.data))';
    ind(startPoint:endPoint)=1;
    ind = logical(ind);
    %----------------------------------------------------------------------
    
    
    L = nnz(ind)-1;
    
    actions = zeros(A, L);
    
    for k = 1:A
        
        variab = genvarname(act_names{k});
        assert(eval([variab '.Rate']) == rates(k));            
        
        if rates(k) == 1
            D = eval([variab '.data(ind)']);
        else
            D = eval([variab '.data']);
            D = D(1:rates(k):end);
            D = D(ind);
        end
                                
        actions(k, :) = double( abs(diff(D)) > 0 )';
        
        
        assert(RALT.Rate == 8);
        AL = RALT.data(1:8:end);
        AL = AL(ind);
        
        %eliminate flights which have big jumps in altitude
        if max(abs(diff(AL))) > 1000
            flag=1;
            break;
        end
        
        %eliminate flights which start at low altitude
        if AL(1) < 4000
            flag=1;
            break;
        end
        
        %eliminate flights which are landing too fast
        if length(AL) < 400
            flag=1;
            break;
        end
                       
        if plotting
            subplot(A+1,1,k);
            plot(D);
            set(gca,'XTickLabel',[]);
            ylabel(act_names{k});
            %fprintf('%.1f ',nnz(actions(k,:))/L*100);
        end
        CH(k, counter) = nnz(actions(k,:))/L*100;
    end
    
    if flag == 1
        flag = 0;
        continue;
    end
    
    if plotting
        subplot(A+1,1,k+1);
        plot(PH.data(ind));
        set(gca,'XTickLabel',[])
        ylabel('PH');
    end
    
%      pause(0.01);
    
     
    tmp{counter} = actions;
    
    %flatten out the data
    R = [1:A]';
    actions = actions.*(R*ones(1,L));
    
    %random permutations at each time stamp
    for j=1:L
        actions(:,j) = actions(randperm(A)',j);
    end
    actions = reshape(actions, numel(actions), 1)';
    actions = actions(actions ~= 0);
    
    res{counter} = actions;  
    
    startInd(counter) = startPoint;
    endInd(counter) = endPoint;
    
    counter = counter + 1;
        
    ind_to_keep = [ind_to_keep; i];        
    
    i
end

%clean unused parts
res(counter:end)=[];
CH(:,counter:end)=[];

startInd(counter:end) = [];
endInd(counter:end) = [];

data = data(ind_to_keep);

ind_bad = 1:N;
ind_bad(ind_to_keep)=[];


% ind = selectFlights_manual(CH);
% 
% CH = CH(:,ind);
% 
% %plot how variabes change from flight to flight
% if 1%plotting
%     figure;
%     
%     for k=1:A
%         subplot(A,1,k);
%         plot(CH(k,:));
%         set(gca,'XTickLabel',[]);
%         ylabel(act_names{k});
%     end
%     
%     figure;
%     
%     for k=1:A
%         subplot(A,1,k);
%         hist(CH(k,:));
%         ylabel(act_names{k});
%     end
%     
%     median(CH,2)
% end




