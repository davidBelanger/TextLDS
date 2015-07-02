clear all;
rng(0);
for T = 1:100
    n = 250;
    A = randn(n,n);
    A = A*A';
    [U , S, V] = svds(A,10);
    S = diag(abs(diag(S)));
    z = randn(1,n)';
    A = diag(z);
    %A = 10*spdiags(z,0,n,n);
    
    X  = A +  U*S*V';
    Xi = inv(X);

    s = mat2MIL(A,U,S,V');

    assertEqual = @(a,b) assert(norm(a(:) - b(:),'inf')/norm(a) < 0.000001);
    assert(norm(Xi - MIL2mat(s),'inf') < 0.00001);
    logdet(Xi)
    MIL_logDet(s)
    d= real(logdet(Xi));
    assert(abs(d - MIL_logDet(s)) < 0.0001);
% 
%     for t = 1:10
%         A = randn(n,10);
%         assertEqual(Xi*A,MIL_rightMultiply(s,A));
% 
%         B = randn(10,n);
% 
%         assertEqual(B*Xi, MIL_leftMultiply(s,B));
% 
%         assertEqual(B*Xi*A , MIL_innerMultiply(s,B,A));
% 
%         prod1 = MIL_innerMultiply(s,B,B');
%         assertEqual(prod1,prod1');
%         
%         C = randn(n,n);
%         assertEqual(trace(Xi*C), MIL_traceOfProduct(s,C));
% 
%     end
%     assertEqual(diag(Xi) , MIL_diag(s));
% 
%     A = randn(10,n);
%     B = randn(n,10);
%     assert(norm(diagOfMatrixProduct(A,B) - diag(A*B)) < 0.0000001);
%     assert(abs(trace(A*B) - matrixProductTrace(A,B))< 0.000001);
%     assert(abs(trace(A*B) - matrixProductTrace(B',A'))< 0.000001);

    
end



