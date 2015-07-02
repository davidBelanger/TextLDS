function params = RandomInitialParams(ape,hidsize)

A = randn(hidsize,hidsize);
params.A = A/max(abs(eig(A)) + .4);
  
params.Q = eye(hidsize) * 1.0;

load('embeddings');
params.C = embeddings*0.0001;


Gamma_0 = ape.y_y_0_with_counts_scaled;


params.Q = eye(hidsize);
Sigma_0  = lyap(params.A,params.Q);
assert(norm(imag(Sigma_0)) < 0.0000001);
Sigma_0 = real(Sigma_0);
assertPSD(Sigma_0); %%
assert(all(isreal(params.C)));
assert(all(isreal(params.A)));
assert(all(isreal(Sigma_0)));


Gamma_0 = ape.y_y_0_with_counts_scaled;%%this is actually the identity, since we assume the data has been whitened
%account for the mean shifting of Sigma_0 in a low rank fashion
%%todo: separate the meanShift from C (b/c we need to get the maxeig
%of C sigma C

%todo: make sure meanShift is sqrt(mu).
n = size(ape.y_y_{1},1);

[Ur Sr Vr] = SparsePlusLRSVD(sparse(n,n),params.C,Sigma_0,params.C,1); %%this is just using MATLAB's API for power iteration

leading_eigenvalue = Sr(1,1);

if(n < 100)
    m = params.C*Sigma_0*params.C';
    evs = eig(m);
    assert(abs(max(evs) - leading_eigenvalue)< 0.00000001);
end
if(leading_eigenvalue > 1)
    disp('WARNING: subspace id produced a noise covariance that is not psd on the dim d-1 subspace orthogonal to mu...manually correcting');
    fprintf('the max eig before manually shrinking was: %f\n',leading_eigenvalue);
    cscale = 0.9/(leading_eigenvalue);
else
    cscale = 1;
end

% C_w_shift = [params.C ape.meanShift];
% P_w_shift = [cscale*Sigma_0 zeros(hidsize,1); zeros(1,hidsize) 1];
% invR =  mat2MIL(Gamma_0,C_w_shift,-cscale*P_w_shift,C_w_shift');

invR =  mat2MIL(Gamma_0,params.C,-cscale*Sigma_0,params.C');

params.invR = invR;

%%these are just some things to make it compatible with downstream code
params.B = randn(hidsize,0);
params.D = randn(n,0);
[m v] = eig(params.A);
params.pi_1 = real(m(:,1));
params.V_1 = Sigma_0;

