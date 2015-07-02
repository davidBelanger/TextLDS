function assertSymmetric(M)

assert(norm(M - M')/norm(M) < 0.0000001);