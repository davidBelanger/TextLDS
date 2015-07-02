function s = mat2MIL(A ,U ,S,V,firstTermInverted)
if(nargin == 4)
    firstTermInverted = false;
end
if(firstTermInverted)
    firstTerm = A;
else
    firstTerm = inv(A);
end

%firstFactor = firstTerm*U;
firstFactor = MIL_rightMultiply(firstTerm,U);
bottleneckTerm = S;
%termInInverse = inv(inv(S) + V*firstTerm*U);
termInInverse = inv(inv(S) + MIL_innerMultiply(firstTerm,V,U));

%secondFactor = termInInverse*V*firstTerm;
secondFactor = MIL_leftMultiply(firstTerm,termInInverse*V);

inputs.A = A;
inputs.U = U;
inputs.S = S;
inputs.V = V;
s = MIL_struct(firstTerm,firstFactor,termInInverse,secondFactor,bottleneckTerm,inputs);
