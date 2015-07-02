function [U S V] = MIL_svd(invR,k)

if(~isstruct(invR))
   [U S V] = svds(invR,k); 
else
    numRows = size(invR.firstTerm,1);
    multiplyByAtranspose = @(x) MIL_leftMultiply(invR,x)';
    multiplyByA = @(x) MIL_rightMultiply(invR,x);
    multiplyBy_A_Atranspose = @(x) multiplyByA(multiplyByAtranspose(x));
    [U S V] = SVDFromLambda(multiplyBy_A_Atranspose,multiplyByAtranspose,numRows,k);
end
