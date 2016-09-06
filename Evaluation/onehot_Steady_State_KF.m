function [LL, x1, xT] = onehot_Steady_State_KF( y_,  params , computeLL)
%%x1 are the filters and xT are the smooths

edgesize = 0;

computeSmooths = nargout > 2;
T = size(y_,2);
fparams = params.filteringParams;

pi_1 = fparams.pi_1;
A = fparams.A;
K_s = fparams.K;
J_s = fparams.J;
invS_logDet = fparams.invS_logDet;
SinvC = fparams.SinvC;
Ctranspose_Sinv_C = fparams.Ctranspose_Sinv_C;
likTerm = fparams.likTerm;
A_KCA = fparams.A_KCA;

outsize = size(fparams.K,2);
hidsize = size(A,1);

x0 = zeros( hidsize, T );
x1 = zeros( hidsize, T );

LL = 0;

x1(:,1) = pi_1;

%x1 is the posterior @ time t not including w_t. x0 includes w_t
for t = 2:T
    yi = y_(t);
    x1(:,t) = A*x0(:,t-1);% + B*u_(:,t-1);
    x0(:,t) = A_KCA*x0(:,t-1) + K_s(:,yi);     
    if(computeLL)
            term1 = likTerm(yi);
            term2 = -2*SinvC(yi,:)*x1(:,t);            
            term3 = x1(:,t)'*Ctranspose_Sinv_C*x1(:,t);
            LL = LL - (term1 + term2 + term3)/2;
    end  
end

if(computeLL)
    LL = LL - (outsize*log(2*pi)/2 - invS_logDet/2)*(T-edgesize);
end


if(computeSmooths)
    xT = zeros( hidsize, T );
    xT(:,T) = x0(:,T);
    for t = (T-1):-1:2
        xT(:,t-1) = x0(:,t-1) + J_s*(xT(:,t) - x1(:,t));
    end
end

end
