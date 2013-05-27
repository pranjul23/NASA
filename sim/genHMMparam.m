function [A0 A O] = genHMMparam(A0_hsmm, A_hsmm, O_hsmm)

A0 = A0_hsmm;
A = A_hsmm(:,:,1);
O = O_hsmm;