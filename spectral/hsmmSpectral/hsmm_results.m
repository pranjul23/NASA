%script to diplay results of testing EM and Spectral HSMM

close all;

ID = 29;

%true joint llikelihood
loc = strcat('../../libdai/examples/data/HSMMlikelihood_true_', num2str(ID), '.txt');
tru = load(loc);

%remove logarithm
tru = exp(tru);

%spectral error
e = (result-tru)./tru;

%Mean absolute deviation
err1_sp = mean(abs(e))

%Root mean squared error
err2_sp = sqrt(mean(e.*e))


N = 11;

err1_em = zeros(N,1);
err2_em = zeros(N,1);

for i = 0:N
    loc = strcat('../../libdai/examples/data/HSMMlikelihood_test_', num2str(ID), '-', num2str(i), '.txt');
    em = load(loc);
    
    %remove logarithm
    em = exp(em);
    
    %EM error
    e = (em-tru)./tru;
    
    %Mean absolute deviation
    err1_em(i+1) = mean(abs(e));
    
    %Root mean squared error
    err2_em(i+1) = sqrt(mean(e.*e));
    
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
