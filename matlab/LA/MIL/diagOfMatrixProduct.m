function d = diagOfMatrixProduct(A,B)
assert(all(size(A) == size(B')));
d = sum(A.*B',2);
