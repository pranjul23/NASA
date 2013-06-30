%given testing results, compute performance characteristic of classifier

function getScore(train, testNorm, testAnom)

perc = [0 5 10 15 20 25 30 35];

N = length(perc);

%get different thresholds which are based on percentiles from train data
thresh = prctile(train, perc);

TP = zeros(N,1);
FP = zeros(N,1);
TN = zeros(N,1);
FN = zeros(N,1);
ACC = zeros(N,1);

for i=1:N
   %detect anomalies in test data
   norm = testNorm < thresh(i); 
   anom = testAnom < thresh(i);
   
   %correctly detected anomalies
   TP(i) = numel(find(anom == 1));
   
   %incorrectly detected anomalies
   FP(i) = numel(find(norm == 1));
   
   %correctly detected normals
   TN(i) = numel(find(norm == 0));
   
   %incorrectly detected normals
   FN(i) = numel(find(anom == 0)); 
   
   %accuracy of prediction
   ACC(i) = (TP(i)+TN(i)) / (TP(i)+TN(i)+FP(i)+FN(i));
end

