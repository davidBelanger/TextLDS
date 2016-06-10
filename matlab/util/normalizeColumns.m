function n = normalizeColumns(matrix)
nr = size(matrix,1);
norms = sqrt(sum(matrix.^2,1));
n = matrix./repmat(norms,[nr 1]);
