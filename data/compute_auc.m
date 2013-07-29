%script to e get AUC scores for classifier
close all

%ID = [40:59 0:39 60:79];
ID = 0:164;

N = length(ID);

%loc = '../libdai/examples/data/resultSim/';
%loc = '../libdai/examples/data/resultChandola/';
loc = '../libdai/examples/data/resultGame2/';

HSMMscores = [];
for i=1:N
    file = strcat(loc, 'HSMMlikelihood_test_', num2str(ID(i)), '.txt');
    %file = strcat(loc, 'stideScores', num2str(ID(i)));
    HSMMscores = [HSMMscores getScores(file)];
end


HMMscores = [];
for i=1:N
    file = strcat(loc, 'HMMlikelihood_test_', num2str(ID(i)), '.txt');
    HMMscores = [HMMscores getScores(file)];
end



%plot results
S = 15;
A = 5;

aucHMM = zeros(A,S);
devHMM = zeros(A,S);

aucHSMM = zeros(A,S);
devHSMM = zeros(A,S);

for i=1:A
    for j=1:S
        aucHMM(i,j) = mean(HMMscores(i, (j-1)*11+1 : j*11));
        devHMM(i,j) = std(HMMscores(i, (j-1)*11+1 : j*11));
        
        aucHSMM(i,j) = mean(HSMMscores(i, (j-1)*11+1 : j*11));
        devHSMM(i,j) = std(HSMMscores(i, (j-1)*11+1 : j*11));
    end

    figure;    
    h=errorbar(2:2+S-1, aucHMM(i,:), devHMM(i,:), 'linewidth', 1.1);
    hc=get(h, 'Children');
    set(hc(2),'Linewidth',1);
    
    title(['Anomaly type ', num2str(i)], 'fontsize', 15);
    xlabel('Number of Hidden States', 'fontsize', 15);
    ylabel('AUC', 'fontsize', 15);
    ylim([0 1]);
    
    hold on;
    h=errorbar(2:2+S-1, aucHSMM(i,:), devHSMM(i,:), 'r-o', 'MarkerFaceColor','r', 'linewidth', 1.1);    
    hc=get(h, 'Children');
    set(hc(2),'Linewidth',1);
    
    xlabel('Number of Hidden States', 'fontsize', 15);
    ylabel('AUC', 'fontsize', 15);    
    ylim([0 1]);
    
    legend('HMM', 'HSMM', 'Location', 'SouthEast');
end



