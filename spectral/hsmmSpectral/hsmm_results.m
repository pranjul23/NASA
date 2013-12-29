%script to diplay results of testing EM and Spectral HSMM

close all;

ID = 30;

N = 4;

for i=1:N
    
    name = ['emTestResults', num2str(ID), '-', num2str(i), '.mat'];
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
    
    
    figure
    plot(err1_em, '-+r', 'linewidth', 2);
    title('MAE');
    
    hold on
    plot(err1_sp*ones(size(err1_em)), '-b', 'linewidth', 2);
    
    legend('EM', 'Spectral');
    xlabel('# of EM iterations', 'fontsize', 15);
    ylabel('Error', 'fontsize', 15);
        
    
    
    figure
    plot(err2_em, '-+r', 'linewidth', 2);
    title('RMSE');
    
    hold on
    plot(err2_sp*ones(size(err2_em)), '-b', 'linewidth', 2);
    
    legend('EM', 'Spectral');
    xlabel('# of EM iterations', 'fontsize', 15);
    ylabel('Error', 'fontsize', 15);
    
    g=9;
    
end
