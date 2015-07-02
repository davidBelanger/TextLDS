function O = MIL_traceOfProduct(invS,A)
 
if(~isstruct(invS))
    
    O = matrixProductTrace(invS,A);
    
else

    if(isstruct(invS.firstTerm))
       firstTrace = MIL_traceOfProduct(invS.firstTerm,A);
    else
       firstTrace = matrixProductTrace(invS.firstTerm,A); 
    end

    O =   firstTrace - matrixProductTrace(invS.firstFactor,(invS.secondFactor*A));

end