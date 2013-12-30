%script to diplay results of testing EM and Spectral HSMM

close all;

ID = 31;

N = 4;

% f1 = figure;
f2 = figure;

col = hsv(5);

err2_em_end = zeros(4,1);
err2_sp_end = zeros(4,1);

for i=1:N
    
    name = ['emTestResults', num2str(ID), '-', num2str(i), '.mat'];    
    load(name);
    
    name = ['emTestResultsSP', num2str(ID), '-', num2str(i), '.mat'];
    load(name);
    
    %remove logarithm
    tru = exp(tru);
    
    %spectral error
    e = (sp-tru)./tru;
    
    %Mean absolute deviation
    err1_sp = mean(abs(e));
    
    %Root mean squared error
    err2_sp = sqrt(mean(e.*e));
    
    M = 8;
    
    err1_em = zeros(M,1);
    err2_em = zeros(M,1);
    
    %remove logarithm
    em = exp(em);
    
    %EM error
    e = (em - tru*ones(1,M))./(tru*ones(1,M));
    
    for j = 1:M
        
        %Mean absolute deviation
        err1_em(j) = mean(abs(e(:,j)));
        
        %Root mean squared error
        err2_em(j) = sqrt(mean(e(:,j).*e(:,j)));
        
    end
    
    err2_em_end(i) = err2_em(end); 
    err2_sp_end(i) = err2_sp;
    
%     figure(f1);
%     plot(err1_em, '-+', 'linewidth', 2, 'color', col(i,:));
%     title('MAE');
%     
%     hold on
%     plot(err1_sp*ones(size(err1_em)), '-', 'linewidth', 2, 'color', col(i,:));
%     
% %     legend('EM', 'Spectral');
%     xlabel('# of EM iterations', 'fontsize', 15);
%     ylabel('Error', 'fontsize', 15);
        
    
    figure(f2);
    plot(err2_em, '-+', 'linewidth', 2, 'color', col(i,:));
    title('RMSE');
    
    hold on
    plot(err2_sp*ones(size(err2_em)), '-', 'linewidth', 2, 'color', col(i,:));
    
%     legend('EM', 'Spectral');
    xlabel('# of EM iterations', 'fontsize', 15);
    ylabel('Error', 'fontsize', 15);            
end


i=N+1;



name = ['emTestResults', num2str(ID), '-', num2str(i), '.mat'];    
load(name);


%spectral error
e = (sp-tru)./tru;

%Mean absolute deviation
err1_sp = mean(abs(e));

%Root mean squared error
err2_sp = sqrt(mean(e.*e));

err2_sp_end(i) = err2_sp;

% figure(f1);
% plot(err1_sp*ones(size(err1_em)), '-', 'linewidth', 2, 'color', col(i,:));

figure(f2);
plot(err2_sp*ones(size(err2_em)), '-', 'linewidth', 2, 'color', col(i,:));

figure;
plot(err2_sp_end,'-+r','linewidth', 2);
hold on
plot(err2_em_end,'-+b','linewidth', 2);

legend('Spectral', 'EM');
xlabel('Training Data Size','fontsize',15);
ylabel('Error','fontsize',15);








