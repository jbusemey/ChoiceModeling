% build simMDF_MEX

clear
clc
na = 3;
G3 = 1.0.*eye(na);
C3 = -1/(na-1)*ones(na,na);   % C matrix for na options
C3 = C3 - diag(diag(C3)) + eye(na);
M3 = [ 1 2 ; 3 4; 5 6]./6;
nd = size(M3,2);
w = ones(nd,1)./nd;
theta1 = 1.0;
sig2 = 1.0;
Ns = 10000;
codegen simMDF.m -args {G3,C3,M3,w,theta1,sig2,Ns}

