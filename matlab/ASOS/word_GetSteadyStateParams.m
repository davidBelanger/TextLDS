function params = word_GetSteadyStateParams(params,ape)


[A,B,C,D,Q,invR,insize,outsize,hidsize,pi_1,V_1] = word_pextract( params );
[V_0_0, V_0_1, V_0_T, V_1_T, K, J, invS] = word_KalmanDoubling( A, Q, C, invR, 50,false);


params.V_0_0 = V_0_0;
params.V_0_1 = V_0_1;
params.V_0_T = V_0_T;
params.V_1_T = V_1_T;
params.K = K;
params.J = J;
params.invS = invS;

