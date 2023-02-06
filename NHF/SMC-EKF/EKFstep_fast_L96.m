function [Pupdated, Xupdated, lW, PSIstep] = EKFstep_fast_L96(F,A,s2x,s2y,h,Tobs,K,y,P0,x0)
% function [Pu, Xu, lW, PSIstep] = EKFstep_fast_RK4(F,A,s2x,s2z,h,Tobs,Oobs_x,z,P0,X0,nosc)
%
% Inputs
% F, A      : F parameter, and coefficients A(1) and A(2)
% s2x   : process noise variances for the slow variables
% s2y   : observation noise variance
% h     : Euler integration step
% Tobs  : observations are collected every Tobs time steps
% K     : one out every K slow oscillators is observed 
% y     : observation
% P0    : prior covariance of the state
% X0    : prior mean of the state
% nosc  : number of oscillators
% 
% Outputs
% Pupdated    : updated covariance matrix
% Xupdated    : updated mean
% lW          : log likelihood
% PSIstep     : ansatz
%

%% Initialisation
% A = [0.07; 0.01];
nosc = size(P0,1);

% Predictive mean and covariance matrix
Ppredictive = zeros(nosc,nosc);
xpredictive = zeros(nosc,1);

% Observation transition matrix
dim_obs = length(1:K:nosc);
Hy_full = eye(nosc);
Hy = Hy_full(1:K:nosc,:);
clearvars Hy_full

%% --Prediction
for k = 1:Tobs   
    
    % Predictive mean
    PSI1 = A'*[x0.'; (x0.').^2];
    k1 = fLx_ansatz(x0,F,PSI1);
    PSI2 = A'*[(x0+1/2.*h.*k1).'; ((x0+1/2.*h.*k1).').^2];
    k2 = fLx_ansatz(x0+1/2.*h.*k1,F,PSI2);
    PSI3 = A'*[(x0+1/2.*h.*k2).'; ((x0+1/2.*h.*k2).').^2];
    k3 = fLx_ansatz(x0+1/2.*h.*k2,F,PSI3);
    PSI4 = A'*[(x0+h.*k3).'; ((x0+h.*k3).').^2];
    k4 = fLx_ansatz(x0+h.*k3,F,PSI4);
    xpredictive= x0 + 1/6*h*(k1+2*k2+2*k3+k4);
    PSIstep = (1/6)*(( PSI1 + (2*PSI2) + (2*PSI3) + PSI4 ));
    
    % Predictive covariance
    J = compute_jacobian2(xpredictive,A,h);
    
%     % Simplification of matrix calculations by blocks  
%     nblocks = nosc/10; n=nosc/nblocks; 
%     blocksH = nosc/20; nH = nosc/blocksH;
%     
%     init = (0:nblocks-1)*nosc/nblocks;
%     aux = zeros(nosc,nosc);
%     
%     % aux = J*P0;
%     aux(1:n,[1:nH,(blocksH-1)*nH+1:nosc])=J( 1:n ,[1:1+n,nosc-1,nosc])*P0([1:1+n,nosc-1,nosc],[1:nH,(blocksH-1)*nH+1:nosc]);
%     for i=2:nblocks-1
%        aux(1+init(i):n+init(i),floor((init(i)-1)/nH)*nH+1:ceil((n+1+init(i))/nH)*nH)=J( 1+init(i):n+init(i) ,init(i)-1:init(i)+1+n)*P0(init(i)-1:init(i)+1+n,floor((init(i)-1)/nH)*nH+1:ceil((n+1+init(i))/nH)*nH);
%     end
%     aux(1+init(nblocks):n+init(nblocks),[1:nH,(blocksH-1)*nH+1:nosc])=J( 1+init(nblocks):n+init(nblocks) ,[1,init(nblocks)-1:end])*P0([1,init(nblocks)-1:end],[1:nH,(blocksH-1)*nH+1:nosc]);
% 
%     % Pp = aux*J';
%     Ppredictive(:,1:n) = aux(:,[1:1+n,nosc-1,nosc])*J( 1:n ,[1:1+n,nosc-1,nosc])';
%     for i=2:nblocks-1
%        Ppredictive(:,1+init(i):n+init(i)) = aux(:,init(i)-1:init(i)+1+n)*J( 1+init(i):n+init(i) ,init(i)-1:init(i)+1+n)';
%     end
%     Ppredictive(:,1+init(nblocks):n+init(nblocks)) = aux(:,[1,init(nblocks)-1:end])*J( 1+init(nblocks):n+init(nblocks) ,[1,init(nblocks)-1:end])';
%     Ppredictive = Ppredictive + s2x*eye(nosc);
    
    Ppredictive = J*P0*J' + s2x*eye(nosc);

    % For the next step
    x0 = xpredictive;
    P0 = Ppredictive;

end %k
clearvars P0 J

%% --Update
% Kalman update eqs
S = Hy*Ppredictive*Hy' + s2y*eye(dim_obs);

if nosc>100
    numberofblocks = nosc/10;
else
     numberofblocks = 2;
end


try   
    S_eig = zeros(dim_obs,1);
    S_inv = zeros(dim_obs,dim_obs);
    for i = 1:numberofblocks
        part = ((i-1)*(dim_obs/numberofblocks))+1:i*dim_obs/numberofblocks;
        aux = S(part,part);
        S_inv(part,part) = inv(aux);
        S_eig(part,1)=eig(aux);
    end

    % S_inv = inv(S);
    % S_eig = eig(S);

    Kalman_Gain = Ppredictive*Hy'*S_inv;
    innov = y(1:K:nosc)' - Hy*xpredictive;
    
    Xupdated = xpredictive + Kalman_Gain*innov;
    Pupdated = (eye(nosc) - Kalman_Gain*Hy)*Ppredictive;

    % Likelihood
    lW = -0.5 * ( sum(log(S_eig)) + innov' * S_inv * innov );

catch
    % when errors in computing inv(S) or eig(S)
    lW = -1e5;
    Xupdated = xpredictive;
    Pupdated = Ppredictive;
    
end
