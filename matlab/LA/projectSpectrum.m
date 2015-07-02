function [A denom] = projectSpectrum(A,maxEig,msg)

if(any(abs(eig(A))> maxEig))
    [M v] = eig(A);
    
    evs = diag(v);
    fprintf('%s: max abs(ev) before scaling: %f\n',msg,max(abs(evs)));
    norms = abs(evs);
    badInds = find(norms > maxEig);
    evs(badInds) = maxEig*evs(badInds)./norms(badInds);
    A = real(M*diag(evs)*inv(M));
    denom = max(norms(badInds)/maxEig);
else
    denom = 1;
end
