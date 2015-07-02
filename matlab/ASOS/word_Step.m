function [o, expt, err, LL] = word_Step( o, params)

%this is the function that implements the ASOS idea.  Note that the large
%commented blocks of code that appear throughout are derivations of the
%formulae that follow them

testingMode = size(params.C, 1) < 200;
computeLikelihood = true;
meanShift = o.meanShift;

[A,B,C,D,Q,invR,insize,outsize,hidsize,pi_1,V_1] = word_pextract( params );

klim = o.klim;
klag = o.klag;
edgesize = o.edgesize;

y_y_ = o.y_y_;
u_u_ = o.u_u_;
u_y_ = o.u_y_;
y_u_ = o.y_u_;

y_lead_ = o.y_lead_;
u_lead_ = o.u_lead_;
y_trail_ = o.y_trail_;
u_trail_ = o.u_trail_;

y_pre_ = o.y_pre_;
u_pre_ = o.u_pre_;
y_post_ = o.y_post_;
u_post_ = o.u_post_;

T = o.T;

F = 1;
total = klim+1+1;

x_x_   = zeros( hidsize, hidsize, total );
x_y_   = zeros( hidsize, outsize, total );
xT_xT_ = zeros( hidsize, hidsize, total );
xT_x_  = zeros( hidsize, hidsize, total );
xT_y_  = zeros( hidsize, outsize, total );
x_u_   = zeros( hidsize, insize, total );
xT_u_  = zeros( hidsize, insize, total );

%secondary versions:
u_x_ = zeros( insize, hidsize, total );
x_xT_ = zeros( hidsize, hidsize, 2 );
u_xT_ = zeros( insize, hidsize, 2 );

[V_0_0, V_0_1, V_0_T, V_1_T, K, J, invS] = word_KalmanDoubling( A, Q, C, invR, 150, false);

if(testingMode)
    [V_0_0_test, V_0_1_test, V_0_T_test, V_1_T_test, K_test, J_test, invS_test] = word_KalmanDoubling( A, Q, C, invR, 50 , true);
    assert(norm(K_test - K,'inf') < 0.0001);
    assert(norm(J_test - J,'inf') < 0.0001);
    R = inv(MIL2mat(invR));
end


fprintf('max ev of J: %f\n',max(abs(eig(J)))); %this will also mess things up if evalues are greater than 1 (smoothing will blow up variances). 
%This will lead to unstable estimation of A in the subsequent iter.


%x_t = K y_t + H x_t-1 + L u_t-1 - KD u_t
%xT_t = J xT_t+1 + P x_t - JB u_t

CB = C*B;
CA = C*A;
H = A - K*CA;
fprintf('max eig of H = A - KCA: %f\n',max(abs(eig(H)))); %note: if H has an e-values greater than 1, everything will diverge. 

L = B - K*CB;
P = eye(hidsize) - J*A;
nKD = -K*D;
nJB = -J*B;



LL = 0;
err = 0;

%% lead and pre 1st order stats

x0_lead_ = zeros(hidsize, edgesize);
V0_lead_ = zeros(hidsize, hidsize, edgesize);
K_lead_ = zeros(hidsize, outsize, edgesize);
J_lead_ = zeros(hidsize, hidsize, edgesize);
x1_lead_ = zeros(hidsize, edgesize+1);
xT_lead_ = zeros(hidsize, edgesize+1);
V1_lead_ = zeros(hidsize, hidsize, edgesize+1);
VT_lead_ = zeros(hidsize, hidsize, edgesize+1);
Vc_lead_ = zeros(hidsize, hidsize, edgesize+1);

x1_pre_ = zeros(hidsize, klag);
x0_pre_ = zeros(hidsize, klag);
xT_pre_ = zeros(hidsize, klag);


x1_lead_(:,1) = pi_1;
V1_lead_(:,:,1) = V_1;
for t = 1:edgesize
    invS_t = inv(C*V1_lead_(:,:,t)*C' + R);
    K_lead_(:,:,t) = V1_lead_(:,:,t)*C'*invS_t;
    V0_lead_(:,:,t) = V1_lead_(:,:,t) - K_lead_(:,:,t)*C*V1_lead_(:,:,t);
    
    y_inov = y_lead_(:,t) - C*x1_lead_(:,t) - D*u_lead_(:,t);
    LL = LL - y_inov'*invS_t*y_inov / 2;
    LL = LL - outsize*log(2*pi)/2 + logdet(invS_t)/2;    
    err = err + y_inov'*y_inov;
    
    x0_lead_(:,t) = x1_lead_(:,t) + K_lead_(:,:,t)*y_inov;
    
    x1_lead_(:,t+1) = A*x0_lead_(:,t) + B*u_lead_(:,t);
    V1_lead_(:,:,t+1) = A*V0_lead_(:,:,t)*A' + Q;
    
end
x1_pre_(:,1) = x1_lead_(:,edgesize+1);

for t = 1:klag
    if t > 1
        x1_pre_(:,t) = A*x0_pre_(:,t-1) + B*u_pre_(:,t-1);
    end
    x0_pre_(:,t) = x1_pre_(:,t) + K*(y_pre_(:,t) - C*x1_pre_(:,t) - D*u_pre_(:,t));
end
xT_pre_(:,klag) = x0_pre_(:,klag);
for t = klag:-1:2
    xT_pre_(:,t-1) = x0_pre_(:,t-1) + J*(xT_pre_(:,t) - x1_pre_(:,t));
end

xT_lead_(:,edgesize+1) = xT_pre_(:,1);
VT_lead_(:,:,edgesize+1) = V_0_T;
for t = (edgesize+1):-1:2
    J_lead_(:,:,t-1) = V0_lead_(:,:,t-1)*A'*inv(V1_lead_(:,:,t));
    xT_lead_(:,t-1) = x0_lead_(:,t-1) + J_lead_(:,:,t-1)*(xT_lead_(:,t) - x1_lead_(:,t));
    VT_lead_(:,:,t-1) = V0_lead_(:,:,t-1) + J_lead_(:,:,t-1)*(VT_lead_(:,:,t) - V1_lead_(:,:,t))*J_lead_(:,:,t-1)';
end

for t = 2:(edgesize+1)
    Vc_lead_(:,:,t) = VT_lead_(:,:,t)*J_lead_(:,:,t-1)';
end


%% trail and post 1st order stats

x0_trail_ = zeros(hidsize, edgesize);
V0_trail_ = zeros(hidsize, hidsize, edgesize);
K_trail_ = zeros(hidsize, outsize, edgesize);
J_trail_ = zeros(hidsize, hidsize, edgesize);
x1_trail_ = zeros(hidsize, edgesize);
xT_trail_ = zeros(hidsize, edgesize);
V1_trail_ = zeros(hidsize, hidsize, edgesize);
VT_trail_ = zeros(hidsize, hidsize, edgesize);
Vc_trail_ = zeros(hidsize, hidsize, edgesize);

x1_post_ = zeros(hidsize, klag);
x0_post_ = zeros(hidsize, klag);
xT_post_ = zeros(hidsize, klag);


x1_post_(:,1) = zeros(hidsize,1);
for t = 1:klag
    if t > 1
        x1_post_(:,t) = A*x0_post_(:,t-1) + B*u_post_(:,t-1);
    end
    x0_post_(:,t) = x1_post_(:,t) + K*(y_post_(:,t) - C*x1_post_(:,t) - D*u_post_(:,t));
end

if edgesize > 0
    x1_trail_(:,1) = A*x0_post_(:,klag) + B*u_post_(:,klag);
    V1_trail_(:,:,1) = V_0_1;
    for t = 1:edgesize
        if t > 1
            x1_trail_(:,t) = A*x0_trail_(:,t-1) + B*u_trail_(:,t-1);
            V1_trail_(:,:,t) = A*V0_trail_(:,:,t-1)*A' + Q;
        end
        invS_t = inv(C*V1_trail_(:,:,t)*C' + R);
        K_trail_(:,:,t) = V1_trail_(:,:,t)*C'*invS_t;
        V0_trail_(:,:,t) = V1_trail_(:,:,t) - K_trail_(:,:,t)*C*V1_trail_(:,:,t);
        
        y_inov = y_trail_(:,t) - C*x1_trail_(:,t) - D*u_trail_(:,t);
        LL = LL - y_inov'*invS_t*y_inov / 2;
        LL = LL - outsize*log(2*pi)/2 + logdet(invS_t)/2;
        err = err + y_inov'*y_inov;
        
        x0_trail_(:,t) = x1_trail_(:,t) + K_trail_(:,:,t)*y_inov;
    end

    xT_trail_(:,edgesize) = x0_trail_(:,edgesize);
    VT_trail_(:,:,edgesize) = V0_trail_(:,:,edgesize);
    for t = edgesize:-1:2
        J_trail_(:,:,t-1) = V0_trail_(:,:,t-1)*A'*inv(V1_trail_(:,:,t));
        xT_trail_(:,t-1) = x0_trail_(:,t-1) + J_trail_(:,:,t-1)*(xT_trail_(:,t) - x1_trail_(:,t));
        VT_trail_(:,:,t-1) = V0_trail_(:,:,t-1) + J_trail_(:,:,t-1)*(VT_trail_(:,:,t) - V1_trail_(:,:,t))*J_trail_(:,:,t-1)';
    end

    
    %TEMP:
    xT_post_(:,klag) = x0_post_(:,klag) + (V_0_0*A'*inv(V1_trail_(:,:,1)))*(xT_trail_(:,1) - x1_trail_(:,1));
    %xT_post_(:,klag) = x0_post_(:,klag) + J*(xT_trail_(:,1) - x1_trail_(:,1));
    %xT_post_(:,klag) = x0_post_(:,klag);
    
    for t = 2:edgesize
        Vc_trail_(:,:,t) = VT_trail_(:,:,t)*J_trail_(:,:,t-1)';
    end
    Vc_trail_(:,:,1) = VT_trail_(:,:,1)*J';
    
else
    xT_post_(:,klag) = x0_post_(:,klag);
end
for t = klag:-1:2
    xT_post_(:,t-1) = x0_post_(:,t-1) + J*(xT_post_(:,t) - x1_post_(:,t));
end

%$$%have we used the incorrect assumption x0_post_end = xT_post_end anywhere?

%% ================================= approximate u_x_{klim+1} =================================
u_x_(:,:,F+klim+1) = ((u_y_(:,:,F+klim) - u_pre_(:,1+klim)*y_pre_(:,1)' )*K' + u_u_(:,:,F+klim+1)*L' + (u_u_(:,:,F+klim)  - u_pre_(:,1+klim)*u_pre_(:,1)')*nKD' + u_pre_(:,1+klim)*x0_pre_(:,1)') * inv(eye(hidsize) - H');
%u_x_(:,:,F+klim+1) = zeros(insize,hidsize);


%% ================================= compute u_x =================================
for k = klim:-1:0
    u_x_(:,:,F+k) = u_x_(:,:,F+k+1)*H' + (u_y_(:,:,F+k) - u_pre_(:,1+k)*y_pre_(:,1)')*K' + u_u_(:,:,F+k+1)*L' + (u_u_(:,:,F+k) - u_pre_(:,1+k)*u_pre_(:,1)')*nKD' + u_pre_(:,1+k)*x0_pre_(:,1)';
end

%% ================================= compute x_u =================================
x_u_(:,:,F+0) = u_x_(:,:,F+0)';
for k = 1:klim+1
    x_u_(:,:,F+k) = H*(x_u_(:,:,F+k-1) - x0_post_(:,end)*u_post_(:,end-k+1)') + K*y_u_(:,:,F+k) + L*(u_u_(:,:,F+k-1) - u_post_(:,end)*u_post_(:,end-k+1)') + nKD*u_u_(:,:,F+k);
end


%% ================================= compute temporary y_x =================================
y_x_(:,:,F+klim+1) = C*A*(-x0_post_(:,end)*x0_post_(:,end-klim)') + CB*(u_x_(:,:,F+klim) - u_post_(:,end)*x0_post_(:,end-klim)') + D*u_x_(:,:,F+klim+1);
for k = klim:-1:0
    y_x_(:,:,F+k) = y_x_(:,:,F+k+1)*H' + (y_y_{F+k})*K' - T*meanShift*(meanShift'*K');
   % y_x_(:,:,F+k) = y_x_(:,:,F+k+1)*H' + (y_y_{F+k} - y_pre_(:,1+k)*y_pre_(:,1)')*K' + y_u_(:,:,F+k+1)*L' + (y_u_(:,:,F+k) - y_pre_(:,1+k)*u_pre_(:,1)')*nKD' + y_pre_(:,1+k)*x0_pre_(:,1)';
end

%% ================================= compute temporary x_y =================================
x_y_(:,:,F+0) = y_x_(:,:,F+0)';
for k = 1:klim
%    x_y_(:,:,F+k) = H*(x_y_(:,:,F+k-1) - x0_post_(:,end)*y_post_(:,end-k+1)') + K*y_y_{F+k} + L*(u_y_(:,:,F+k-1) - u_post_(:,end)*y_post_(:,end-k+1)') + nKD*u_y_(:,:,F+k);
    x_y_(:,:,F+k) = H*(x_y_(:,:,F+k-1)) + K*y_y_{F+k} - T*(K*meanShift)*meanShift';

end


%% ================================= compute x_x_(:,:,F+klim) =================================

%x_x_(:,:,F+klim) = x_x_(:,:,F+klim+1)*H' + (x_y_(:,:,F+klim) - x0_(:,1+klim)*y_(:,1))*K' + x_u_(:,:,F+klim+1)*L' + (x_u_(:,:,F+klim) - x0_pre_(:,1+klim)*u_pre_(:,1)')*nKD' + x0_(:,1+klim)*x0_(:,1)';
%x_x_(:,:,F+klim+1) = H*(x_x_(:,:,F+klim) - x0_(:,T)*x0_(:,T-klim)') + K*y_x_(:,:,F+klim+1) + L*(u_x_(:,:,F+klim) - u_post_(:,end)*x0_post_(:,end-klim)') + nKD*u_x_(:,:,F+klim+1);
%approx: y_x_(:,:,F+klim+1) ~= C*A*(x_x_(:,:,F+klim) - x0_(:,T)*x0_(:,T-klim)') + CB*(u_x_(:,:,F+klim) - u_post_(:,end)*x0_post_(:,end-klim)') + D*u_x_(:,:,F+klim+1)
%x_x_(:,:,F+klim+1) = H*(x_x_(:,:,F+klim) - x0_(:,T)*x0_(:,T-klim)') + K*( C*A*(x_x_(:,:,F+klim) - x0_(:,T)*x0_(:,T-klim)') + CB*(u_x_(:,:,F+klim) - u_post_(:,end)*x0_post_(:,end-klim)') + D*u_x_(:,:,F+klim+1) ) + L*(u_x_(:,:,F+klim) - u_post_(:,end)*x0_post_(:,end-klim)') + nKD*u_x_(:,:,F+klim+1);
%x_x_(:,:,F+klim+1) = A*(x_x_(:,:,F+klim) - x0_(:,T)*x0_(:,T-klim)') + B*(u_x_(:,:,F+klim) - u_post_(:,end)*x0_post_(:,end-klim)');
%x_x_(:,:,F+klim) = ( A*(x_x_(:,:,F+klim) - x0_(:,T)*x0_(:,T-klim)') + B*(u_x_(:,:,F+klim) - u_post_(:,end)*x0_post_(:,end-klim)') )*H' + (x_y_(:,:,F+klim) - x0_(:,1+klim)*y_(:,1))*K' + x_u_(:,:,F+klim+1)*L' + (x_u_(:,:,F+klim) - x0_pre_(:,1+klim)*u_pre_(:,1)')*nKD' + x0_(:,1+klim)*x0_(:,1)';
const = ( A*( - x0_post_(:,end)*x0_post_(:,end-klim)') + B*(u_x_(:,:,F+klim) - u_post_(:,end)*x0_post_(:,end-klim)') )*H' + (x_y_(:,:,F+klim) - x0_pre_(:,1+klim)*y_pre_(:,1)')*K' + x_u_(:,:,F+klim+1)*L' + (x_u_(:,:,F+klim) - x0_pre_(:,1+klim)*u_pre_(:,1)')*nKD' + x0_pre_(:,1+klim)*x0_pre_(:,1)';
%norm(x0_post_(:,end)*x0_post_(:,end-klim)')
%fprintf('norms %f %f %f\n',(1/T)*norm(( A*( - x0_post_(:,end)*x0_post_(:,end-klim)') + B*(u_x_(:,:,F+klim) - u_post_(:,end)*x0_post_(:,end-klim)') )*H'),(1/T)*norm((x_y_(:,:,F+klim) - x0_pre_(:,1+klim)*y_pre_(:,1)')*K'),(1/T)*norm(x0_pre_(:,1+klim)*x0_pre_(:,1)'));
%x_x_(:,:,F+klim) = dlyap(A, H', const);
%x_x_(:,:,F+klim) = SylvDoubling(A, H', const, 50);


M = H^(2*klim + 1);
x_x_(:,:,F+klim) = zeros(hidsize);
tmp = const;
for i = 1:35 %35 seems like it should be enough
    %tmp = dlyap( A, H', tmp ); %you should use this instead of the below statement if you have access to the dlyap function (it comes in the control toolbox)
    tmp = SylvDoubling( A, H', tmp, 50 );
    x_x_(:,:,F+klim) = x_x_(:,:,F+klim) + tmp;
    if norm(vec(tmp)) < 1e-12*norm(vec(x_x_(:,:,F+klim)))
        break;
    end
    tmp = M*tmp'*A'*C'*K';
end
if(any(isnan(x_x_(:))))
    disp('here');
end
assert(~any(isnan(x_x_(:))));

%i
%{
for i = 1:500
    x_x_(:,:,F+klim) = A*x_x_(:,:,F+klim)*H' + M*x_x_(:,:,F+klim)'*A'*C'*K' + const;
end

%save for future use:
%o.x_x_klim = x_x_(:,:,F+klim);
%}
%norm(A)
%norm(H)
%norm(M)*norm(A'*C'*K')

%% ================================= update y_x =================================
M = eye(hidsize);
for k = klim:-1:0
    M = M*H;
    y_x_(:,:,F+k) = y_x_(:,:,F+k) + CA*x_x_(:,:,F+klim)*M';
end
%% ================================= update x_y =================================
for k = 0:klim
    x_y_(:,:,F+k) = x_y_(:,:,F+k) + M*x_x_(:,:,F+klim)'*CA';
    M = M*H;
end

%% ================================= compute x_x =================================

for k = (klim-1):-1:0
    x_x_(:,:,F+k) = x_x_(:,:,F+k+1)*H' + (x_y_(:,:,F+k) - x0_pre_(:,1+k)*y_pre_(:,1)')*K' + x_u_(:,:,F+k+1)*L' + (x_u_(:,:,F+k) - x0_pre_(:,1+k)*u_pre_(:,1)')*nKD' + x0_pre_(:,1+k)*x0_pre_(:,1)';
end

%xT_t = J xT_t+1 + P x_t - JB u_t

%% ================================= compute xT_x =================================
xT_x_(:,:,F+klim) = x_x_(:,:,F+klim);
for k = (klim-1):-1:0
    xT_x_(:,:,F+k) = J*xT_x_(:,:,F+k+1) + P*(x_x_(:,:,F+k) - x0_post_(:,end)*x0_post_(:,end-k)') + nJB*(u_x_(:,:,F+k)- u_post_(:,end)*x0_post_(:,end-k)') + xT_post_(:,end)*x0_post_(:,end-k)';
end

%% ================================= compute x_xT =================================
x_xT_(:,:,F+0) = xT_x_(:,:,F+0)';


%% ================================= compute xT_u =================================
xT_u_(:,:,F+klim) = x_u_(:,:,F+klim);
for k = (klim-1):-1:0
    xT_u_(:,:,F+k) = J*xT_u_(:,:,F+k+1) + P*(x_u_(:,:,F+k) - x0_post_(:,end)*u_post_(:,end-k)') + nJB*(u_u_(:,:,F+k) - u_post_(:,end)*u_post_(:,end-k)') + xT_post_(:,end)*u_post_(:,end-k)';
end

%% ================================= compute u_xT =================================
u_xT_(:,:,F+0) = xT_u_(:,:,F+0)';


%% ================================= compute xT_xT =================================
%xT_t = J xT_t+1 + P x_t - JB u_t

%xT_xT_(:,:,F+k) = J*xT_xT_(:,:,F+k+1) + P*(x_xT_(:,:,F+k) - x0_(:,T)*xT_(:,T-k)') + nJB*(u_xT_(:,:,F+k) - u_(:,T)*xT_(:,T-k)') + xT_(:,T)*xT_(:,T-k)';
%xT_xT_(:,:,F+k+1) = (xT_xT_(:,:,F+k) - xT_(:,k+1)*xT_(:,1)')*J' + xT_x_(:,:,F+k+1)*P' + xT_u_(:,:,F+k+1)*nJB';
%xT_xT_(:,:,F+k) = J*( (xT_xT_(:,:,F+k) - xT_(:,k+1)*xT_(:,1)')*J' + xT_x_(:,:,F+k+1)*P' + xT_u_(:,:,F+k+1)*nJB' ) + P*(x_xT_(:,:,F+k) - x0_(:,T)*xT_(:,T-k)') + nJB*(u_xT_(:,:,F+k) - u_(:,T)*xT_(:,T-k)') + xT_(:,T)*xT_(:,T-k)';
const = J*( (- xT_pre_(:,0+1)*xT_pre_(:,1)')*J' + xT_x_(:,:,F+0+1)*P' + xT_u_(:,:,F+0+1)*nJB' ) + P*(x_xT_(:,:,F+0) - x0_post_(:,end)*xT_post_(:,end-0)') + nJB*(u_xT_(:,:,F+0) - u_post_(:,end)*xT_post_(:,end-0)') + xT_post_(:,end)*xT_post_(:,end-0)';
%xT_xT_(:,:,F+0) = dlyap( J, const ); %you should use this instead of the below statement if you have access to the dlyap function (it comes in the control toolbox)

%todo: symmetrize const?
%assertSymmetric(xT_xT_(:,:,F+0));


xT_xT_(:,:,F+0) = LyapDoubling( J, const, 50 );

%Smooth out potential numerical problems:
%xT_xT_(:,:,F+0) = (xT_xT_(:,:,F+0) + xT_xT_(:,:,F+0)')/2;

%xT_xT_(:,:,F+k) = (xT_xT_(:,:,F+k-1) - xT_(:,1+k-1)*xT_(:,1)')*J' + xT_x_(:,:,F+k)*P' + xT_u_(:,:,F+k)*nJB';
xT_xT_(:,:,F+1) = (xT_xT_(:,:,F+1-1) - xT_pre_(:,1+1-1)*xT_pre_(:,1)')*J' + xT_x_(:,:,F+1)*P' + xT_u_(:,:,F+1)*nJB';

%% ================================= compute xT_y =================================
xT_y_(:,:,F+klim) = x_y_(:,:,F+klim);
for k = (klim-1):-1:0
     xT_y_(:,:,F+k) = J*xT_y_(:,:,F+k+1) + P*(x_y_(:,:,F+k));
%    xT_y_(:,:,F+k) = J*xT_y_(:,:,F+k+1) + P*(x_y_(:,:,F+k) - x0_post_(:,end)*y_post_(:,end-k)') + nJB*(u_y_(:,:,F+k) - u_post_(:,end)*y_post_(:,end-k)') + xT_post_(:,end)*y_post_(:,end-k)';
end


%% Compute likelihood only if you have R, which you don't get in the first iter if you initialized with SSID
if(computeLikelihood)
offsetTerm1 = -MIL_rightMultiply(invR,meanShift)'*meanShift;
firstTrace = MIL_traceOfProduct(invR,y_y_{F+0}) + T*offsetTerm1;
%secondTrace = trace((C'*(invR*y_y_{F+0}*invR)*C)*invS.termInInverse);

%secondTrace = matrixProductTrace((invS.Ctranspose_invR*(y_y_{F+0} - meanShift*meanShift')*invS.invR_C),invS.termInInverse);
offsetTerm2 = -matrixProductTrace((invS.Ctranspose_invR*meanShift)*(meanShift'*invS.invR_C),invS.termInInverse);
secondTrace = matrixProductTrace((invS.Ctranspose_invR*(y_y_{F+0})*invS.invR_C),invS.termInInverse) + T*offsetTerm2;

term1 = firstTrace - secondTrace;

% term2 = trace((A'*C')*invR*y_x_(:,:,F+1)) - trace((A'*(C'*invR*C))*invS.termInInverse*(C'*invR*y_x_(:,:,F+1)));

firstTrace = matrixProductTrace(A'*invS.Ctranspose_invR,y_x_(:,:,F+1));
secondTrace = matrixProductTrace((A'*invS.Ctranspose_invR_C)*invS.termInInverse,(invS.Ctranspose_invR*y_x_(:,:,F+1)));

term2 = firstTrace - secondTrace;

% term3 = trace((A'*(C'*invR*C)*A)*x_x_(:,:,F+0)) - trace((A'*(C'*invR*C))*invS.termInInverse*(C'*invR*C)*(A*x_x_(:,:,F+0)));
%here, there are no intermediate matrices that have have an outsize dimension size
firstTrace = trace((A'*(invS.Ctranspose_invR_C)*A)*x_x_(:,:,F+0));
secondTrace = trace((A'*(invS.Ctranspose_invR_C))*invS.termInInverse*(invS.Ctranspose_invR_C)*(A*x_x_(:,:,F+0)));

term3 =  firstTrace - secondTrace;

LLincr = -1/2*(term1 -2*term2 + term3);
LLincr = LLincr + (T-2*edgesize)*(-outsize*log(2*pi) + invS.logDet)/2;
LL = LL + LLincr;
end

if(testingMode && computeLikelihood)

% predmat = y_y_{F+0} - y_pre_(:,1)*y_pre_(:,1)' - 2*y_x_(:,:,F+1)*CA' + CA*(x_x_(:,:,F+0) - x0_post_(:,end)*x0_post_(:,end)')*CA' ...
%         - 2*y_u_(:,:,F+1)*CB' + 2*CB*(u_x_(:,:,F+0) - u_post_(:,end)*x0_post_(:,end)' )*CA' + CB*(u_u_(:,:,F+0) - u_post_(:,end)*u_post_(:,end)')*CB' ...
%         - 2*(y_u_(:,:,F+0) - y_pre_(:,1)*u_pre_(:,1)')*D' + 2*D*u_x_(:,:,F+1)*CA' + 2*D*u_u_(:,:,F+1)*CB' + D*(u_u_(:,:,F+0) - u_pre_(:,1)*u_pre_(:,1)')*D' ...
%         + (y_pre_(:,1) - C*x1_pre_(:,1) - D*u_pre_(:,1))*(y_pre_(:,1) - C*x1_pre_(:,1) - D*u_pre_(:,1))';     
    predmat = y_y_{F+0}  - T*meanShift*meanShift'- 2*y_x_(:,:,F+1)*CA' + CA*(x_x_(:,:,F+0))*CA'; %todo: is this not properly accounting for meanShift?
    
    LLincr_test = -1/2*vec(invS_test)'*vec(predmat);
    LLincr_test = LLincr_test + (T-2*edgesize)*(-outsize*log(2*pi) + logdet(invS_test))/2;
    fprintf('LL1 = %f and LL2 = %f\n',LLincr_test,LLincr);
    assert(abs((LLincr_test - LLincr)/LLincr) < 0.00001,['diff = ' num2str(abs(LLincr_test - LLincr))]);
    err = err + trace(predmat);
end

if(isnan(LL))
    disp('here');
end

Ex_x_0 = xT_xT_(:,:,F+0) + xT_lead_(:,1:edgesize)*xT_lead_(:,1:edgesize)' + xT_trail_*xT_trail_' + V_0_T*(T - 2*edgesize) + sum( VT_lead_(:,:,1:edgesize), 3 ) + sum( VT_trail_, 3);
Ex_x_1 = xT_xT_(:,:,F+1) + xT_lead_(:,2:edgesize+1)*xT_lead_(:,1:edgesize)' + xT_trail_(:,2:edgesize)*xT_trail_(:,1:edgesize-1)' + V_1_T*(T - 2*edgesize - 1) + sum(Vc_lead_, 3) + sum(Vc_trail_, 3);
if edgesize > 0
    Ex_x_1 = Ex_x_1 + xT_trail_(:,1)*xT_post_(:,klag)';
end
    
Ey_x_0 = xT_y_(:,:,F+0)' + y_lead_*xT_lead_(:,1:edgesize)' + y_trail_*xT_trail_';
Eu_x_0 = u_xT_(:,:,F+0) + u_lead_*xT_lead_(:,1:edgesize)' + u_trail_*xT_trail_';
Ex_u_1 = xT_u_(:,:,F+1) + xT_lead_(:,2:edgesize+1)*u_lead_' + xT_trail_(:,2:edgesize)*u_trail_(:,1:edgesize-1)';
if edgesize > 0
    Ex_u_1 = Ex_u_1 + xT_trail_(:,1)*u_post_(:,klag)';
end

Ex_.start = xT_lead_(:,1);
if edgesize > 0
    Ex_.end = xT_trail_(:,edgesize);
else
    Ex_.end = xT_post_(:,klag);
end

Exx_.start = Ex_.start*Ex_.start' + VT_lead_(:,:,1);
if edgesize > 0
    Exx_.end = Ex_.end*Ex_.end' + VT_trail_(:,:,edgesize);
else
    Exx_.end = Ex_.end*Ex_.end' + V_0_T;
end

expt = estruct( Ex_x_0, Ex_x_1, Ey_x_0, Eu_x_0, Ex_u_1, Ex_, Exx_, [] );

end

