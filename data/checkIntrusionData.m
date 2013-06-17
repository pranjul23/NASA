
%directory to read training data from
dir_train = 'lpr-mit-live-normal/';

list_train = dir(dir_train);
numSeqTrain = length(list_train)-3;


%directory to read anomaly data from
dir_test = 'lpr-mit-live-anomaly/';

list_test = dir(dir_test);
numSeqTest = length(list_test)-3;

A = [];
mapObj = containers.Map('KeyType','int32','ValueType','int64');

for i=1:numSeqTrain
    obs = load(strcat(dir_train, list_train(i+3).name));
    obs = obs(:,2);
    
    for k=1:length(obs)
        if mapObj.isKey(obs(k))
            mapObj(obs(k)) = mapObj(obs(k)) + 1;
        else
            mapObj(obs(k)) = 1;
        end
    end
    
    A = union(A,unique(obs)); 
    i
end


for i=1:numSeqTest
    obs = load(strcat(dir_test, list_test(i+3).name));
    obs = obs(:,2);
    
    for k=1:length(obs)
        if mapObj.isKey(obs(k))
            mapObj(obs(k)) = mapObj(obs(k)) + 1;
        else
            mapObj(obs(k)) = 1;
        end
    end
    
    A = union(A,unique(obs));    
    i
end


A = A + 1;
map = (-1)*ones(183,1);

for i=1:length(A)
   map(A(i)) = i-1;
end

% %remove repetitive test sequences
% A = cell(0);
% ind = 1;
% good_ind = [];
%     
% 
% for i=1:numSeqTest
%     obs = load(strcat(dir_test, list_test(i+3).name));
%     obs = obs(:,2)';
%         
%     found = 0;
%     
%     if(isempty(A))
%         found = 0;
%     end
%     
%     for j=1:length(A)
%         if length(obs)==length(A{j})
%             if(all(A{j} == obs))
%                 found = 1;  
%             end
%         end
%     end
%     
%     if ~found
%         good_ind = [good_ind; i];
%         A{ind} = obs;
%         ind = ind + 1;
%     end
% end

% 
%                 11098
%                  3675
%               1283705
%               1262187
%                 48898
%                 76315
%                  6526
%                  2624
%                  7793
%                  9306
%                 13110
%                  5525
%                  3930
%                  8588
%                    12
%                  3703
%                 19154
%                    16
%                 13051
%                 17581
%                  3702
%                 32685
%                  3702
%                  3685
%                  6699
%                  5740
%                 12822
%                  3564
%                 13346
%                 78406
%                 39203
%                 39202
%                  9423
%                 10210
%                  5737
%                  7381
%                  5740
%                  3703
%                  3851
%                  3687
%                  6092



