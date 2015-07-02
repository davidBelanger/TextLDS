function O = MIL_leftMultiply(invS,A)

if(~isstruct(invS))
    O = A*invS;
else
    if(isstruct(invS.firstTerm))
        firstProduct = MIL_leftMultiply(invS.firstTerm,A);
    else
        firstProduct = A*invS.firstTerm;
    end
    O =   firstProduct - (A*invS.firstFactor)*invS.secondFactor;
end
