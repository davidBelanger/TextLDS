function [V_0_0, V_0_1, V_0_T, V_1_T, K, J, invS] = word_KalmanDoubling( A, Q, C, invR, T ,computeInvS)

%A is the state dynamics matrix
%Q is the dynamics noise covariance

%C is the prediction matrix
%invR is the inverse prediction noise covariance

numhid = size(A,1);

%% Forwards phase
Phi = A';
Psi = MIL_innerMultiply(invR,C',C);
assertPSD(Psi); 
Theta = Q;
for t = 1:T
    mat = inv(eye(numhid) + Psi*Theta);
    incr = Phi'*Theta*mat*Phi;
    Theta = Theta + incr;
    if norm(vec(incr)) < 1e-12*norm(vec(Theta))
        break;
    end
    Psi = Psi + Phi*mat*Psi*Phi';
    Phi = Phi*mat*Phi;
end
%t

%% Asymptotic quantities from forwards phase
V_0_1 = Theta;
assertPSD(V_0_1);
invV = inv(V_0_1);
%S = (D + C*V_0_1*C')

invR_withoutShift = invR;%mat2MIL(invR,meanShift,1,meanShift',true);

invS = mat2MIL(invR_withoutShift,C,V_0_1,C',true);

invS.logDet = MIL_logDet(invS); %%todo: fix this

%a = MIL2mat(invS);
%assertEqual(logdet(a),invS.logDet);
%%todo: these should probably be w.r.t invR_withoutShift
Ctranspose_invR_C =  MIL_innerMultiply(invR_withoutShift,C',C);
Ctranspose_invR = MIL_leftMultiply(invR_withoutShift,C');
invR_C = MIL_rightMultiply(invR_withoutShift,C);

invS.Ctranspose_invR_C =  Ctranspose_invR_C;
invS.Ctranspose_invR = Ctranspose_invR;
invS.invR_C = invR_C;

%for testing
logDetR = -MIL_logDet(invR_withoutShift);
logDetS = logdet(V_0_1) + logDetR  + logdet(inv(V_0_1) + Ctranspose_invR_C); %%using the Matrix Determinant Lemma
%assert(abs((real(logDetS )- (-MIL_logDet(invS)))/logDetS) < 0.0001,['diff = ' num2str(abs((logDetS - (-MIL_logDet(invS)))/logDetS))]);


%K = V_0_1*C'*invS; %%this is the formula. 
Ctranspose_invS = MIL_leftMultiply(invS,C');
K = V_0_1*Ctranspose_invS;


%%this is just for testing
if(computeInvS)
    invS_t = inv(C*V_0_1*C' + inv(MIL2mat(invR)));
    Kdb = V_0_1*C'*invS_t;
    if(~(norm(K - Kdb,'inf')/norm(K) < 0.0001))
        disp('here');
    end
    assert(norm(K - Kdb,'inf')/norm(K) < 0.0001);
    
    assert(abs(logdet(invS_t) - invS.logDet) < 0.00001);
    invS = invS_t;
end
V_0_0 = V_0_1 - K*C*V_0_1; %%it seems like we could also solve for V_0_0 with doubling (i.e. fixed point iteration), but this is what Martens does

J = V_0_0*A'*inv(V_0_1);

D = V_0_0 - J*V_0_1*J';


%% Backwards phase
Phi = J';
Theta = D;

for t = 1:T
    incr = Phi'*Theta*Phi;
    Theta = Theta + incr;
    if norm(vec(incr)) < 1e-8*norm(vec(Theta))
        break;
    end
    Phi = Phi^2;
end
%t
V_0_T = Theta;

assert(all(isreal(J(:))));

%% Cross phase
V_1_T = V_0_T*J';

%assertSymmetric(V_1_T);
%assertSymmetric(V_0_T);
