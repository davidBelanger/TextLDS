function  [Ur Sr Vr] = SparsePlusLRTroppSVD(sparsePart,U,S,V,k,powerIters,usePowerMethod,oversampling,projection)
%%todo: really we should be able to handle left and right projections
if(nargin == 8)
    projection = @(x) x;
end


    multiplyByA = @(x) projection((sparsePart*x +  U*S*(V'*x)));
    multiplyByAtranspose = @(x) projection((sparsePart'*x +  V*S*(U'*x)));

    numCols = size(sparsePart,2);
    %todo: shouldn't need to pass projection in
    [Ur,Sr,Vr] = fsvdFromLambda(multiplyByA,multiplyByAtranspose,numCols, k, powerIters, usePowerMethod,oversampling,projection);
    
    
    Ur = projection(Ur);
    Vr = projection(Vr);
    
    