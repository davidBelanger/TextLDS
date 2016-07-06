function U = orthogonalizeColumns(C)
[U S V] = svd(C,'econ');
