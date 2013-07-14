%given testing results, compute performance characteristic (AUC) of classifier

function res = getScores(file)

data = load(file);

%remove infinities
tmp = data;
tmp(tmp == -inf) = [];
data(data == -inf) = min(tmp) - 0.5;

N = 5;
L = 10;
res = zeros(N,1);

%prepare labels
% [anomalies; normals]
labels = [ones(L,1);zeros(L,1)];


for i=1:N
    %create scores
    %[anomalies; normals]
    scores = [data((i-1)*L+1 : i*L); data(N*L+1 : (N+1)*L)];
    
    [~, ~, ~, AUC] = perfcurve(labels, scores, 0);
    res(i) = AUC;
end