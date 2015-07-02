function [Ur sv] = MIL_eigs(invR,k)

if(~isstruct(invR))
   [Ur sv] = eigs(invR,k); 
else
    numRows = size(invR.firstFactor,1);
%    multiplyByAtranspose = @(x) MIL_leftMultiply(invR,x')';
    multiplyByA = @(x) MIL_rightMultiply(invR,x);
 %   multiplyBy_A_Atranspose = @(x) multiplyByA(multiplyByAtranspose(x));
    [Ur, sv, flag] = eigs(multiplyByA,numRows, k, 'lm');
end

