function ind = getLeftInd(start, Nhid, Dmax, numObs)


ind = zeros(1, numObs);
ind(1) = start - (Dmax-1);
ind(end) = start;

p = 1;
for i = 2:numObs-1
    ind(i) = ind(i-1) + (Nhid^p - Nhid^(p-1));
    p=p+1;
end
