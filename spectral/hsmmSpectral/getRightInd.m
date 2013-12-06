function ind = getRightInd(start, Nhid, Dmax, numObs)


ind = zeros(1, numObs);
ind(end) = start+Dmax-1;
ind(1) = start;

p = 1;
for i = numObs-1 : -1 : 2
    ind(i) = ind(i+1) - (Nhid^p - Nhid^(p-1));
    p=p+1;
end
