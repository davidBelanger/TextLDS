function [A denom] = scaleSpectrum(A,maxEig)

if(any(eig(A)> maxEig))
    [M v] = eig(A);
    evs = diag(v);
    denom = (max(abs(evs)))/maxEig;
    fprintf('max eig is %f. dividing through by %f\n',max(abs(evs)),denom);
    A = A / denom;
else
    denom = 1;
end
