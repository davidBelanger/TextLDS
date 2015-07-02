function O = MIL_innerMultiply(invS,A,B)

if(~isstruct(invS))
    O = A*invS*B;
else
    
    if(isstruct(invS.firstTerm))
        firstProduct = MIL_innerMultiply(invS.firstTerm,A,B);
    else
        firstProduct = A*invS.firstTerm*B;
    end
    
    O =  firstProduct - (A*invS.firstFactor)*(invS.secondFactor*B);
end