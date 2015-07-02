function O = MIL_rightMultiply(invS,A)
if(~isstruct(invS))
    O = invS*A;
else
    if(isstruct(invS.firstTerm))
        firstProduct = MIL_rightMultiply(invS.firstTerm,A);
    else
        firstProduct = invS.firstTerm*A;
    end
    O =   firstProduct - invS.firstFactor*(invS.secondFactor*A);
end
