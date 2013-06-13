function [A0 A O] = genHMMparam(A0_hsmm, A_hsmm, O_hsmm)

A0 = A0_hsmm;
A =  A_hsmm;
O = O_hsmm;

%data for running murphyk HMM code
save('murphykHMMinit', 'A0', 'A', 'O');