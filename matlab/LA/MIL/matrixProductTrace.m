function s = matrixProductTrace(A,B)

[n m] = size(A);
[n2 m2] = size(B);

assert(m == n2);
if(n < m)
    s = sum(sum(A'.*B,2));
else
    s = sum(sum(B.*A',2));
end
