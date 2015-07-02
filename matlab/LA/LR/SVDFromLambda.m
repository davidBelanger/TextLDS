function [Ur Sr Vr] = SVDFromLambda(multiplyBy_A_Atranspose,multiplyByAtranspose,numRows,k)
  
    
    [Ur, sv, flag] = eigs(multiplyBy_A_Atranspose,numRows, k, 'lm');
    Sr = sqrt(sv);
    
    Vr = multiplyByAtranspose(Ur*spdiags(1./diag(Sr),0,k,k));
    