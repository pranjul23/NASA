function res = computeTrueProbTensor(test, A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);



Tr1 = tensor(embedA1(A, D1, A1));
Dur = tensor(embedD(D));
Tr = tensor(embedA(A));
Trend = tensor(A);

e = 1e-0;


I = ttt(tensor(pinv(O)'), tensor(O), 1, 1);
I = tensor(eye(x));
P = ttt(tensor(O), I, 2, 1);

P = ttt(Trend, P, 1, 2);
%% =============================================

O6O7 = tensor(computeO6O7D4X4(A1, A, D1, D, O));
O6O7m = tenmat(O6O7, [1 2]);
O6O7m(:) = pinv(O6O7m.data)';
O6O7inv = tensor(O6O7m);

I = ttt(tensor(O6O7), O6O7inv, [1 2], [1 2]);
M = tensor(zeros(d,x,d,x));
Mm = tenmat(M, [1 2]);
Mm(:) = eye(x*d) + e*rand(x*d);
I = tensor(Mm);
P = ttt(I, P, [1 2], [2 1]);

P = ttt(Dur, P, [1 2], [1 2]);

O5O6 = computeO5O6D3X4(A1, A, D1, D, O);
O5O6m = tenmat(O5O6, [1 2]);
O5O6m(:) = pinv(O5O6m.data)';
O5O6inv = tensor(O5O6m);

I = ttt(tensor(O5O6), O5O6inv, [1 2], [1 2]);
M = tensor(zeros(d,x,d,x));
Mm = tenmat(M, [1 2]);
Mm(:) = eye(x*d) + e*rand(x*d);
I = tensor(Mm);
P = ttt(I, P, [1 2], [2 1]);

P = ttt(Tr, P, [4 2], [1 2]);

I = ttt(tensor(pinv(O)'), tensor(O), 1, 1);
I = tensor(eye(x));
tmp = ttt(tensor(O), I, 2, 1);

P = ttt(P, tmp, [1], [2]);
%% =============================================

O5O6 = computeO5O6D3X3(A1, A, D1, D, O);
O5O6m = tenmat(O5O6, [1 2]);
O5O6m(:) = pinv(O5O6m.data)';
O5O6inv = tensor(O5O6m);

I = ttt(tensor(O5O6), O5O6inv, [1 2], [1 2]);
M = tensor(zeros(d,x,d,x));
Mm = tenmat(M, [1 2]);
Mm(:) = eye(x*d) + e*rand(x*d);
I = tensor(Mm);
P = ttt(I, P, [1 2], [2 1]);

P = ttt(Dur, P, [1 2], [1 2]);

O4O5 = computeO4O5D2X3(A1, A, D1, D, O);
O4O5m = tenmat(O4O5, [1 2]);
O4O5m(:) = pinv(O4O5m.data)';
O4O5inv = tensor(O4O5m);

I = ttt(tensor(O4O5), O4O5inv, [1 2], [1 2]);
M = tensor(zeros(d,x,d,x));
Mm = tenmat(M, [1 2]);
Mm(:) = eye(x*d) + e*rand(x*d);
I = tensor(Mm);
P = ttt(I, P, [1 2], [2 1]);

P = ttt(Tr, P, [4 2], [1 2]);

I = ttt(tensor(pinv(O)'), tensor(O), 1, 1);
I = tensor(eye(x));
tmp = ttt(tensor(O), I, 2, 1);

P = ttt(P, tmp, [1], [2]);
%% =============================================

O4O5 = computeO4O5D2X2(A1, A, D1, D, O);
O4O5m = tenmat(O4O5, [1 2]);
O4O5m(:) = pinv(O4O5m.data)';
O4O5inv = tensor(O4O5m);

I = ttt(tensor(O4O5), O4O5inv, [1 2], [1 2]);
M = tensor(zeros(d,x,d,x));
Mm = tenmat(M, [1 2]);
Mm(:) = eye(x*d) + e*rand(x*d);
I = tensor(Mm);
P = ttt(I, P, [1 2], [2 1]);

P = ttt(Dur, P, [1 2], [1 2]);

O3O4 = computeO3O4D1X2(A1, A, D1, D, O);
O3O4m = tenmat(O3O4, [1 2]);
O3O4m(:) = pinv(O3O4m.data)';
O3O4inv = tensor(O3O4m);

I = ttt(tensor(O3O4), O3O4inv, [1 2], [1 2]);
M = tensor(zeros(d,x,d,x));
Mm = tenmat(M, [1 2]);
Mm(:) = eye(x*d) + e*rand(x*d);
I = tensor(Mm);
P = ttt(I, P, [1 2], [2 1]);

P = ttt(Tr1, P, [4 2], [1 2]);

I = ttt(tensor(pinv(O)'), tensor(O), 1, 1);
I = tensor(eye(x));
tmp = ttt(tensor(O), I, 2, 1);

P = ttt(P, tmp, [1], [2]);
%% =============================================

I = ttt(tensor(pinv(O)'), tensor(O), 1, 1);
I = tensor(eye(x));
tmp = ttt(tensor(O), I, 2, 1);

P = ttt(P, tmp, [1], [2]);
%% =============================================

N = size(test,1);

res = zeros(N,1);

for i = 1:N
    
    res(i) = P(test(i,5), test(i,4), test(i,3), test(i,2), test(i,1));
    
end