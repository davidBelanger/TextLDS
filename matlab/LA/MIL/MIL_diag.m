function O = MIL_diag(invS)
 
if(~isstruct(invS))
    
    O = diag(invS);
    
else
    
    O = MIL_diag(invS.firstTerm) - diagOfMatrixProduct(invS.firstFactor,invS.secondFactor);

end

