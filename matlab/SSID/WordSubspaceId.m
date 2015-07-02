function params = WordSubspaceId(ape,hidsize,useTropp,oversampling,numPowerIters,numObs,order,simpleMethodForC,transformForSVD)
%%%NOTE: numerous aspects of this algorithm expect that the data has been
%%%whitened when passed in

n = size(ape.y_y_{1},1);

if(~isempty(transformForSVD))
    unTransform = inv(transformForSVD);
    transform = transformForSVD;
    blockTransform =  kron(speye(order),transform);
else
    transform = speye(n);
    unTransform = speye(n);
    blockTransform = speye(n*order);
end

transformedMeanShift = (transform*ape.meanShift);

%this is the term for the rank one whitened mean-shifting of the Hankel Matrix
block_meanShift = repmat(transformedMeanShift,[order 1]);

AA = ConstructHankel(ape,order,transform,numObs);

disp('doing SVD');
tic;
usePowerMethod = true; %alternative is to use Krylov Subspaces
    if(useTropp)
        [U, S, V] = SparsePlusLRTroppSVD(AA,block_meanShift,-1,block_meanShift,hidsize,numPowerIters,usePowerMethod,oversampling);
    else
        [U, S, V] = SparsePlusLRSVD(AA,block_meanShift,-1,block_meanShift,hidsize);
    end


fprintf('time to do svd: %f\n',toc);
clear AA;


sq_S = sqrt(S);
Delta = sq_S*(V')*blockTransform';
G = Delta(:,(end -n+1):end);

D1 = Delta(:,1:(end -n));
D2 = Delta(:,(n+1) : end);
params.A = D1*pinv(D2);


%%
sq_S = sqrt(S);
Delta = sq_S*(V');
D1 = Delta(:,1:(end -n));
D2 = Delta(:,(n+1) : end);
params.A = D1*pinv(D2);
%%

Lambda = blockTransform*U*sq_S;

if(simpleMethodForC)
    params.C = Lambda(1:n,:);
else
    Abar = [];
    for i = 0:(order-1)
        Abar = [Abar (params.A^i)'];
    end
    Lambdabar = cell2mat(cellfun(@(x) x',mat2cell(Lambda,n*ones(1,order),hidsize),'UniformOutput',false))';
    params.C = Lambdabar * pinv(Abar);
end
toc;

maxEig = 0.98;
e = max(abs(eig(params.A)));
fprintf('max eig of A = %f\n',e);

%assert(e < maxEig,num2str(e)); %% we should be alarmed if this fails
%[params.A, denom] = scaleSpectrum(params.A,maxEig);
%params.C = params.C*denom;
params.A = projectSpectrum(params.A,maxEig,'A-SSID');

params.Q = eye(hidsize);
Sigma_0  = lyap(params.A,params.Q);
assert(norm(imag(Sigma_0)) < 0.00001);
Sigma_0 = real(Sigma_0);
assertPSD(Sigma_0);
assert(all(isreal(params.C)));
assert(all(isreal(params.A)));
assert(all(isreal(Sigma_0)));

Gamma_0 = ape.y_y_0_with_counts_scaled;

    whitener = spdiags(1./sqrt(spdiags(Gamma_0)),0,n,n);
    [Ur Sr Vr] = SparsePlusLRSVD(sparse(n,n),whitener*params.C,Sigma_0,whitener*params.C,1); %%this is just using MATLAB's API for power iteration
    leading_eigenvalue = Sr(1,1);
    
    if(leading_eigenvalue > 1)
        disp('WARNING: subspace id produced a noise covariance that is not psd on the dim d-1 subspace orthogonal to mu...manually correcting');
        fprintf('the max eig before manually shrinking was: %f\n',leading_eigenvalue);
        cscale = maxEig/(leading_eigenvalue);
    else
        cscale = 1;
    end
    invR =  mat2MIL(Gamma_0,params.C,-cscale*Sigma_0,params.C');
%     [u s] = MIL_eigs(invR,3);
%     assert(all(diag(s) > 0));
%     fprintf('minimal eig of invR = %f\n',min(diag(s).^(-1)));
    params.invR = invR;

%%these are just some things to make it compatible with downstream code
params.B = randn(hidsize,0);
params.D = randn(n,0);
[m v] = eig(params.A);
params.pi_1 = real(m(:,1));
params.V_1 = Sigma_0;

