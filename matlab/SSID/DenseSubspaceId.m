function params = DenseSubspaceId(ape,hidsize,useTropp,oversampling,numPowerIters,numObs,order,simpleMethodForC,whitenForSVD)

%%todo: this is still using block whitening from the left to whiten the
%%right. we could be something more exact
n = size(ape.y_y_{1},1);
assert(~whitenForSVD,'probably should not be whitening, since the input data is likely whitened in word space');
if(whitenForSVD)
    whiten = ape.whiten;
else
    whiten = speye(n);
end
unWhiten = inv(whiten);
blockWhitening =  kron(speye(order),whiten);
blockUnWhitening =  kron(speye(order),unWhiten);

AA = ConstructHankel(ape,order,whiten,numObs);

disp('doing SVD');
tic;
usePowerMethod = true; %alternative is to use Krylov Subspaces

if(useTropp)
    [U, S, V] = fsvd(AA, hidsize, numPowerIters, usePowerMethod,oversampling);
else
    [U, S, V] = svds(AA,hidsize);
end

fprintf('time to do svd: %f\n',toc);
%clear AA;

sq_S = sqrt(S);
Delta = sq_S*(V')*blockWhitening; %todo: or should this be un-whitening?
D1 = Delta(:,1:(end -n));
D2 = Delta(:,(n+1) : end);
params.A = D1*pinv(D2);


Lambda = blockWhitening*U*sq_S;%todo: or should this be un-whitening?

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
assert(norm(imag(Sigma_0),'inf') < 0.0000001);
Sigma_0 = real(Sigma_0);
assertPSD(Sigma_0);
assert(all(isreal(params.C)));
assert(all(isreal(params.A)));
assert(all(isreal(Sigma_0)));

Gamma_0 = ape.y_y_0_scaled;


R = Gamma_0 - params.C*Sigma_0*params.C';
minEig = 0.0001;
params.R = projectPSD(R,minEig);
assertPSD(params.R);

%%these are just some things to make it compatible with downstream code
params.B = randn(hidsize,0);
params.D = randn(n,0);
[m v] = eig(params.A);
params.pi_1 = real(m(:,1));
params.V_1 = Sigma_0;

