function assertPSD(M)
assertSymmetric(M);
e = eig(M);
if(~(all(e> 0)))
    e(find(e <= 0))
end
assert(all(e > 0),num2str(min(e)));
