function [params, LLhist] = word_ApproxEMforLDS( iparams, priors, numiters, updateTransMat, updateHiddenCovMat, updateObservedCovMat,updateInit,  firstIterToUpdateEmbeddings,firstIterToUpdateTransitions,ape, T, covStyle)
LLhist = zeros(1,numiters);
[A,B,C,D,Q,invR,insize,outsize,hidsize,pi_1,V_1] = word_pextract( iparams );

y_y_0 = ape.y_y_0;
y_u_0 = ape.y_u_0;
u_0_0 = ape.u_u_0;
u_u_0_2 = ape.u_u_0_2;
inv_u_u_0 = ape.inv_u_u_0;
inv_u_u_0_2 = ape.inv_u_u_0_2;
debugMode = size(y_y_0,1) < 20;

if ~isempty( priors )
    error( 'priors not implemented' );
    %[sigA, nQ, invVQ, nR, invVR] = rextract( priors );
end
%start of loop:
updateEmbeddings = @(iter) iter > firstIterToUpdateEmbeddings;
updateTransition = @(iter) iter > firstIterToUpdateTransitions;

%since the bohning covariance doesn't depend on inference, might as well set it early
if(strcmp(covStyle,'bohning'))
    %invR = mat2MIL(0.5*speye(outsize),ones(outsize,1),-0.5/(outsize+1),ones(1,outsize));
    invR = mat2MIL(2*speye(outsize),ones(outsize,1),1,ones(1,outsize)); 
    Psi = MIL_innerMultiply(invR,C',C);
    assertPSD(Psi); 
end
            
for iter = 1:numiters
    tic;    
    %% Compute approximate E-step    
    currentParams =  word_pstruct( A, B, C, D, Q, invR, pi_1, V_1 );
    [ape, expt, err, LL] = word_Step( ape, currentParams);
    LLhist(iter) = LL;
    disp( ['Iteration ' num2str(iter) ':'] );
    disp( ['Approximate LL: ' num2str(LL)]);
    checkLikelihood = (strcmp(covStyle,'full') || iter > 2) && (~strcmp(covStyle,'bohning'));
    if(iter > 1 && LLhist(iter - 1) > LLhist(iter) && checkLikelihood)
        disp('warning: likelihood went down');
        %assert(false);
        %break;
    end
    
    %% M-Step
    [Ex_x_0, Ex_x_1, Ey_x_0, Eu_x_0, Ex_u_1, Ex_, Exx_] = eextract( expt );
    
    newC = ( (1/T)*(Ey_x_0 - y_u_0*inv_u_u_0*Eu_x_0))*inv( (1/T)*(Ex_x_0 - Eu_x_0'*inv_u_u_0*Eu_x_0 ));
    
   % newD = (y_u_0 - newC*Eu_x_0')*inv_u_u_0;
    
    Eu_x_0x = Eu_x_0;%- u_(:,T)*Ex_.end';
    newA = (Ex_x_1 - Ex_u_1*inv_u_u_0_2*Eu_x_0x)/((Ex_x_0 - Exx_.end) - Eu_x_0x'*inv_u_u_0_2*Eu_x_0x);
    
    if(max(abs(eig(newA))) > 1)
        disp(['A is unstable' num2str(max(abs(eig(newA))))]);
    end
%    assert(max(abs(eig(newA))) < 1,['A is unstable' num2str(max(abs(eig(newA))))]);
    
    if(true)
        %%this is what you should do is A is unstable. Ideally we would
        %%never have this problem though, so we're throwing an exception
        %%above
        maxEig = 0.9999;
        [Aproj, denom] = projectSpectrum(newA,maxEig,'A');
        cob  = newA \ Aproj;
        newA = Aproj;
        newC = newC*inv(cob);
        assert(all(isreal(newA(:))));
        assert(all(isreal(newC(:))));
    end
    
    if updateTransMat && updateTransition(iter)
        A = newA;
    end
    
    %    B = newB;
    if(updateEmbeddings(iter))
        C = newC;
    end
    
    %    D = newD;
    
    if(updateObservedCovMat);
        if(strcmp(covStyle,'full'))
            Ex_x_0 = 0.5*(Ex_x_0 + Ex_x_0');
            C_M = [C (1/T)*Ey_x_0];
            C_M_2 = [C'; (1/T)*Ey_x_0'];
            Bm = [(1/T)*Ex_x_0 -eye(hidsize);-eye(hidsize) zeros(hidsize,hidsize)];
            inv_Sigma_0 = inv((1/T)*y_y_0);
            %inv_Sigma_0 = mat2MIL((1/T)*y_y_0,meanShift,-1,meanShift');
            invR =  mat2MIL(inv_Sigma_0,C_M,Bm,C_M_2,true);%%see writeup for why we ignore the meanshift term here.
           % invR_withoutShift =  mat2MIL(inv((1/T)*y_y_0),C_M,Bm,C_M_2,true);

            assertSymmetric((1/T)*Ex_x_0);
            assertSymmetric(Bm);
            Psi = MIL_innerMultiply(invR,C',C);
            assertPSD(Psi); 
            if(debugMode)
                debug_invR_without_shift = mat2MIL((1/T)*y_y_0,C_M,Bm,C_M_2);
                debug_m = MIL2mat(debug_invR_without_shift);
                m = MIL2mat(invR);
                assertPSD(m);
                assert(norm(debug_m  - m) < 0.000001);
            end
        elseif(strcmp(covStyle,'spherical'))
            %%this is a terrible hack; for if you're using 4-whitened data
            sqrtFreq = spdiags(inv((1/T)*y_y_0)).^(-1);
            s = sum(sqrtFreq.^2);
            assert(abs(-1 + s) < 0.0000001,num2str(s));
            freq = sqrtFreq.^2;
            %at every timestep, you observe 1/(f^4)
            diag_term = sum(freq.*(freq.^(-.25)));
            %%%%terrible hack ends here
            
            d = diag_term -  2*(1/T)*matrixProductTrace(C,Ey_x_0') + (1/T)*matrixProductTrace(Ex_x_0 , C'*C);
            assert(d> 0);
            d = d/outsize; 
            assert(d> 0);
            R = d*speye(outsize);
            fprintf('new spherical variance = %f; if not divided = %f\n',d,d*outsize);
            invR = inv(R);
        elseif(strcmp(covStyle,'diagonal'))
            Ex_x_0 = 0.5*(Ex_x_0 + Ex_x_0');
            C_M = [C (1/T)*Ey_x_0];
            C_M_2 = [C'; (1/T)*Ey_x_0'];
            Bm = [(1/T)*Ex_x_0 -eye(hidsize);-eye(hidsize) zeros(hidsize,hidsize)];
            V = size(C,1);
            Sigma_0 = (1/T)*spdiags(y_y_0) - ape.meanShift(:).^2; %%unlike for the full covariance case, we actually include the meanShift term. That's because we no longer are defining an inverse of a rank-deficient matrix
            diagTerm2 = diagOfMatrixProduct(C_M*Bm,C_M_2)';
            diagTerm = Sigma_0(:) - diagTerm2(:);                        
            invR = inv(spdiags(diagTerm,0,V,V));
        elseif(strcmp(covStyle,'bohning'))
            %invR = mat2MIL(0.5*speye(outsize),ones(outsize,1),-0.5*1/(outsize+1),ones(1,outsize));
            invR = mat2MIL(2*speye(outsize),ones(outsize,1),1,ones(1,outsize)); 
        else
           assert(false,'invalid cov style'); 
        end
    end
    if(updateHiddenCovMat); Q = 1/(T)*((Ex_x_0 - Exx_.start) - newA*Ex_x_1' - B*Ex_u_1');end

    
    if updateInit
        pi_1 = Ex_.start;
        V_1 = Exx_.start - Ex_.start*Ex_.start';
    end
    
    %% Stamp out some potential numerical problems
    
    V_1 = (V_1 + V_1')/2;
    Q = (Q + Q')/2;
    
    
    %% Sanity Checks

    if(~isstruct(invR) &&  any(eig(inv(invR)) < -0.001 ))
        min(sort(eig(inv(invR))))
                error('invR developed negative eigenvalues')
        [U S V ] = svd(inv(invR));
        ridge = 0.01;
        invR = inv(U*diag(max(diag(S),ridge))*V');
    end
    
    if any(eig(Q) < 0 )
        sort(eig(Q))
        error('Q developed negative eigenvalues')
    end
   
end

params = word_pstruct( A, B, C, D, Q, invR, pi_1, V_1 );
params.Ey_x_0 = Ey_x_0/T;
params.Ex_x_0 = Ex_x_0/T;
end
%%this is stuff for modeling with a spherical covariance and if you're
%%dealing with non-diagonal observation covariances
%         if(uniformCovariance)
%             diag_term = ape.wNormSquared/T;
%
%             %switch this when you're running on non-word data (i.e.
%             %observations aren't indicator vectors)
%             % diag_term = 1;
%
%             d = diag_term -  2*(1/T)*matrixProductTrace(Ey_x_0',C) + (1/T)*matrixProductTrace(Ex_x_0 , C'*C);
%             d = d/outsize;
%             assert( d> 0);
%             R = d*speye(outsize);
%             invR = inv(R);
%         elseif(false)
%             %             d = (1/T)*(frequencies - 2*sum(Ey_x_0.*C,2) + sum((C*Ex_x_0).*C,2)); %%note this covariance contains both the xbar-xbar' term and the V_0_T term
%             %             R = sparse(1:outsize,1:outsize,d);
%             if(~updatedSomething)
%                 cross = -Ey_x_0*C';
%                 R = 1/T*(y_y_0 + C*Ex_x_0*C'  + cross + cross');
%             else
%                 R = 1/T*(y_y_0 - C*Ey_x_0' - D*y_u_0');
%             end
%             R = (R + R')/2;
%             invR = inv(R);
%         elseif(false)
%             %ignore the cross terms, since we don't know how to handle
%             %them with an MILstruct
%             invR = mat2MIL((1/T)*y_y_0,C,Ex_x_0,C');
%         else
%


%% Computing the M-step objective (for debugging purposes)
function E = computeObj(y_y_0, u_u_0, u_u_0_2, y_u_0, u_, T, params, priors, expt)

[A,B,C,D,Q,R,insize,outsize,hidsize,pi_1,V_1] = pextract( params );
[Ex_x_0, Ex_x_1, Ey_x_0, Eu_x_0, Ex_u_1, Ex_, Exx_] = eextract( expt );
if ~isempty( priors )
    [sigA, nQ, invVQ, nR, invVR] = rextract( priors );
end

invR = inv(R);
invQ = inv(Q);
invV_1 = inv(V_1);

E = 0;

cross = -Ey_x_0*C' + D*Eu_x_0*C' - y_u_0*D';
E = E - vec(invR)'*vec(y_y_0 + C*Ex_x_0*C' + D*u_u_0*D' + 2*cross );

Eu_x_0x = Eu_x_0 - u_(:,T)*Ex_.end';
cross = -Ex_u_1*B' + B*Eu_x_0x*A' - Ex_x_1*A';
E = E - vec(invQ)'*vec((Ex_x_0 - Exx_.start) + B*u_u_0_2*B' + A*(Ex_x_0 - Exx_.end)*A' + 2*cross);

E = E - vec(invV_1)'*vec(Exx_.start - 2*pi_1*Ex_.start' + pi_1*pi_1');
E = E - T*logdet(R) - (T)*logdet(Q) - logdet(V_1);
E = E - T*(outsize + hidsize)*log(2*pi);

if ~isempty( priors )
    E = E - 1/(2*sigA^2)*vec(A - eye(hidsize))'*vec(A - eye(hidsize));
    E = E - vec(invR)'*vec(invVR) - (nR + outsize + 1)*logdet(R); %missing some constants
    E = E - vec(invQ)'*vec(invVQ) - (nQ + hidsize + 1)*logdet(Q); %missing some constants
    E = E - (2*hidsize^2*log(sigA) + hidsize^2*log(2*pi));
end

E = E/2;

if imag(E) ~= 0
    error( 'Log of negative number while computing E (bad)' );
end
end
