function d = MIL_logDet(invS)

%%this is an implementation of the matrix determinant lemma (from Wikipedia)
%%det(A + UVW') = det(inv(A)) + V'inv(A)U)det(W)det(A)
%%invS represents inv(A + UVW'), so logdet(invS) = -logdet(A + UVW')

if(~isstruct(invS))
    if(issparse(invS))
        d = sparseDiagLogdet(invS); %we're assuming that the only sparse matrices that get passed in are diagonal ones
    else
        d = logdet(invS);
    end
else
    if(isfield(invS,'logDet'))
        d= invS.logDet;
    else
        d = -(logdet(invS.bottleneckTerm) - MIL_logDet(invS.firstTerm)  - logdet(invS.termInInverse));
    end
end
d = real(d);