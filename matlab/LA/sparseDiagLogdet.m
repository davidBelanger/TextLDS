function d  = sparseDiagLogdet(A)

a = diag(A);
d = sum(log(a));