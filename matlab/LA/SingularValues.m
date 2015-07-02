function v = SingularValues(A)

[U S V] = svd(A);
v = diag(S)';
