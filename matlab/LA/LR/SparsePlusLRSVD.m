function  [Ur Sr Vr] = SparsePlusLRSVD(sparsePart,U,S,V,k,projection)

if(nargin == 5)
    projection = @(x) x;
end

multiplyByA = @(x) projection((sparsePart*x +  U*S*(V'*x)));
multiplyByAtranspose = @(x) projection((sparsePart'*x +  V*S*(U'*x)));

multiplyBy_A_Atranspose = @(x) multiplyByA(multiplyByAtranspose(x));

numRows = size(sparsePart,1);
[Ur Sr Vr] = SVDFromLambda(multiplyBy_A_Atranspose,multiplyByAtranspose,numRows,k);


