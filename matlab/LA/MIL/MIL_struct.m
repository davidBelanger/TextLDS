 function invS = MIL_struct(firstTerm,firstFactor,termInInverse,secondFactor,bottleneckTerm,origInputs)
if(nargin == 5)
    origInputs = [];
end
 
        invS = struct;
        
        invS.firstTerm = firstTerm;
        invS.firstFactor = firstFactor;
        invS.termInInverse = termInInverse;
        invS.secondFactor = secondFactor;
        invS.bottleneckTerm = bottleneckTerm;
        invS.origInputs = origInputs;
        
        