%script to diplay results of testing EM and Spectral HSMM

% close all;

ID = 32;

N = 5;

% f1 = figure;
f2 = figure;

col = lines(N);
ltype = {'--', '-.', ':', '--', '-'};

err2_em_end = zeros(N,1);
err2_sp_end = zeros(N,1);

% folder = 'orinjadeResults';
folder = 'dioResults';
%folder = 'maximusResults';

iter = 0:10:70;
dataSizes = [500, 1000, 5000, 10^4, 10^5];
% Ldata = [100, 1000, 10^4, 10^5];

for i=1:N
    
    name = [folder '/emTestResults', num2str(ID), '-', num2str(i), '.mat'];
%     name = [folder '/simTestResults', num2str(ID), '-', num2str(i), '.mat'];
    load(name);
    
%un-comment for maximusResults
%     if i<=4
%         name = [folder '/emTestResultsSP', num2str(ID), '-', num2str(i), '.mat'];
%         load(name);
%     else
%         name = [folder '/emTestResultsExtra', num2str(ID), '-', num2str(1), '.mat'];
%         load(name);
%     end
    
    %remove logarithm
    tru = exp(tru);
%     tru = tru/100;
    
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
%     em = em/100;
    
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
    plot(iter, err2_sp*ones(size(err2_em)), ltype{i}, 'linewidth', 2, 'color', col(i,:));
    
    hold on    
    h = plot(iter, err2_em, [ltype{i} 'o'], 'linewidth', 2, 'color', col(i,:),'MarkerFaceColor',col(i,:));
    
    set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
    
    xlabel('EM iterations', 'fontsize', 17);
    ylabel('Error', 'fontsize', 17);       
    
%     set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15);
    set(gca,'YTick',get(gca,'YTick'),'fontsize',15);

end

h = legend('500', '1000', '5000', '10^4', '10^5');
set(h, 'FontSize',15);


figure;
semilogx(dataSizes, err2_sp_end, '-', 'linewidth', 2, 'MarkerFaceColor', 'b');
hold on
h = semilogx(dataSizes, err2_sp_end, '-o', 'linewidth', 2, 'MarkerFaceColor', 'b');
set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');

semilogx(dataSizes, err2_em_end, 'r-.', 'linewidth', 2, 'MarkerFaceColor', 'r');
h = semilogx(dataSizes, err2_em_end, 'r-.o', 'linewidth', 2, 'MarkerFaceColor', 'r');
set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');

xlim([500 10^5]);

xlabel('Log (Num. train. samples)', 'fontsize', 17);
ylabel('Error', 'fontsize', 17);

h = legend('Spectral' ,'EM');
set(h, 'FontSize',15);

% set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15);
set(gca,'YTick',get(gca,'YTick'),'fontsize',15);
% i=N+1;
% 
% 
% 
% name = ['emTestResults', num2str(ID), '-', num2str(i), '.mat'];    
% load(name);
% 
% 
% %spectral error
% e = (sp-tru)./tru;
% 
% %Mean absolute deviation
% err1_sp = mean(abs(e));
% 
% %Root mean squared error
% err2_sp = sqrt(mean(e.*e));
% 
% err2_sp_end(i) = err2_sp;
% 
% % figure(f1);
% % plot(err1_sp*ones(size(err1_em)), '-', 'linewidth', 2, 'color', col(i,:));
% 
% figure(f2);
% plot(err2_sp*ones(size(err2_em)), '-', 'linewidth', 2, 'color', col(i,:));
% 
% figure;
% plot(err2_sp_end,'-+r','linewidth', 2);
% hold on
% plot(err2_em_end,'-+b','linewidth', 2);
% 
% legend('Spectral', 'EM');
% xlabel('Training Data Size','fontsize',15);
% ylabel('Error','fontsize',15);


st30=[0.5176 0.6201 2.0596 3.3634 38.3607];
st32=[0.2514 0.2720 2.1452 3.3960 38.4607];

et30=[315+313+312+312+312+312+313 635+640+629+624+624+635+641 3169+3127+3149+3133+3137+3193+3179 6290+6311+6358+6312+6355+6365 380000];
et32=[263+264+264+263+265+263+264 530+529+530+529+530+529+529 2588+2649+2647+2608+2645+2623+2624 5228+5236+5217+5235+5226+5239+5215 366000];

figure;

loglog(dataSizes, st30, '-', 'linewidth', 2, 'MarkerFaceColor', 'b');
hold on
h=loglog(dataSizes, st30, '-o', 'linewidth', 2, 'MarkerFaceColor', 'b');
set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');

loglog(dataSizes, et30, 'r-.', 'linewidth', 2, 'MarkerFaceColor', 'r');
hold on
h=loglog(dataSizes, et30, 'r-.o', 'linewidth', 2, 'MarkerFaceColor', 'r');
set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');

xlim([500 10^5]);

xlabel('Log (Num. train. samples)', 'fontsize', 17);
ylabel('Log (Runtime), sec', 'fontsize', 17);

h = legend('Spectral' ,'EM');
set(h, 'FontSize',15);

set(gca,'YTick',get(gca,'YTick'),'fontsize',15);
    
figure;

loglog(dataSizes, st32, '-', 'linewidth', 2, 'MarkerFaceColor', 'b');
hold on
h=loglog(dataSizes, st32, '-o', 'linewidth', 2, 'MarkerFaceColor', 'b');
set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');

loglog(dataSizes, et32, 'r-.', 'linewidth', 2, 'MarkerFaceColor', 'r');
hold on
h=loglog(dataSizes, et32, 'r-.o', 'linewidth', 2, 'MarkerFaceColor', 'r');
set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');

xlim([500 10^5]);

xlabel('Log (Num. train. samples)', 'fontsize', 17);
ylabel('Log (Runtime), sec', 'fontsize', 17);

h = legend('Spectral' ,'EM');
set(h, 'FontSize',15);

set(gca,'YTick',get(gca,'YTick'),'fontsize',15);




