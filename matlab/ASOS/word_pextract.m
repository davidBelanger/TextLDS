function [A,B,C,D,Q,invR,insize,outsize,hidsize,pi_1,V_1] = word_pextract( params )

A = params.A;
B = params.B;
C = params.C;
D = params.D;
Q = params.Q;
invR = params.invR;
insize = size(B,2);
outsize = size(C,1);
hidsize = size(A,1);
pi_1 = params.pi_1;
V_1 = params.V_1;
end
