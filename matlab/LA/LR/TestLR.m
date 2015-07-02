clear;
  

for i = 1:50
    
    numRows = 750;
    numCols = 750;
    numObs = 7500;
    is = randi(numRows,numObs,1);
    js = randi(numCols,numObs,1);
    values = randn(numObs,1);
    Asparse = sparse(is,js,values,numRows,numCols);
    
    b = randn(numRows,numCols);
    k = 3;
    [U S V] = svds(b,k);
    A = Asparse + U*S*V';
    mu = rand(numRows,1);
    mu = mu/sum(mu);
    mus = sqrt(mu);
    mus'*mus;
   % projection = @(x) x - mus*(mus'*x);

%    projection = @(x) x;
    
%     z = randn(numCols,1);
%            z0 = projection(z);
% 
%     for t = 1:10
%        zi = projection(z);
%        assert(norm(z0 - zi) < 0.000000001); 
%     end
    [Ur Sr Vr] = SparsePlusLRSVD(Asparse,U,S,V,k);

    powerIters = 50;
    usePowerMethod = true;
    oversampling = 10;
    
    [Ur2 Sr2 Vr2] = SparsePlusLRTroppSVD(Asparse,U,S,V,k,powerIters,usePowerMethod,oversampling);
    
%     
%     assert(all(mus'*Ur < 0.0000001));
%     assert(all(mus'*Ur2 < 0.0000001));
%     
%      
%     [U S V] = svds(A,k);
     x = Ur*Sr*Vr';

%     x = Ur*Sr*Vr';
     x2 = Ur2*Sr2*Vr2';
     norm(x - x2)/norm(x)


end