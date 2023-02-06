function [Xupdated, xmean_updated, lW, PSIstep] = EnKFstep_fast_RK4_SQMC( F, A, s2x, s2y, h, Tobs, K, y, X0, samples_noise)
% function [Xupdated, xmean_updated, lW, PSIstep] = EnKFstep_fast_RK4_SQMC( F, A, s2x, s2y, h, Tobs, K, y, X0, noise_step)
%
% Inputs
% F, A      : F parameter, and coefficients A(1) and A(2)
% s2x       : process noise variances for the slow variables
% s2y       : observation noise variance
% h         : Euler integration step
% Tobs      : observations are collected every Tobs time steps
% K         : one out every K slow oscillators is observed 
% y         : observation
% X0        : prior ensemble of the state
% samples_noise     : random (SMC) or quasirandom (SQMC) point set
% 
% Outputs
% Xupdated          : updated ensemble
% xmean_updated     : updated ensemble
% lW                : log likelihood
% PSIstep           : ansatz


%% Initialisation

[nosc,M_ensemble] = size(X0);
PSIstep_ensemble = zeros(nosc,M_ensemble);

% Observation transition matrix
dim_obs = length(1:K:nosc);
Hy_full = eye(nosc);
Hy = Hy_full(1:K:nosc,:);
clearvars Hy_full


%% -- Prediction    
for t=1:Tobs
    
    Xpredictive = zeros(nosc,M_ensemble);
    for j=1:M_ensemble       
        
        PSI1 = A'*[X0(:,j).'; (X0(:,j).').^2];
        k1 = fLx_ansatz(X0(:,j),F,PSI1);
        
        aux = X0(:,j)+1/2.*h.*k1+1/2.*sqrt(s2x).*sqrt(2).*erfinv(2.*(samples_noise(j*4-3,:)')-1);
        PSI2 = A'*[(aux).'; ((aux).').^2];
        k2 = fLx_ansatz(aux,F,PSI2);
        
        aux = X0(:,j)+1/2.*h.*k2+1/2.*sqrt(s2x).*sqrt(2).*erfinv(2.*(samples_noise(j*4-2,:)')-1);
        PSI3 = A'*[(aux).'; ((aux).').^2];
        k3 = fLx_ansatz(aux,F,PSI3);
        
        aux = X0(:,j)+h.*k3+sqrt(s2x).*sqrt(2).*erfinv(2.*(samples_noise(j*4-1,:)')-1);
        PSI4 = A'*[(aux).'; ((aux).').^2];
        k4 = fLx_ansatz(aux,F,PSI4);
        
        Xpredictive(:,j)= X0(:,j) + 1/6*( h*(k1+2*k2+2*k3+k4)+  sqrt(s2x).*sqrt(2).*erfinv(2.*(samples_noise(j*4,:)')-1) );
        PSIstep_ensemble(:,j) = (1/6)*(( PSI1 + (2*PSI2) + (2*PSI3) + PSI4 ));
        
     
    end
      
    X0 = Xpredictive;
    
end %t
clearvars noiseinstep

%% -- Update

PSIstep = (1/M_ensemble)*PSIstep_ensemble*ones(M_ensemble,1);
Xmean = (1/M_ensemble)*Xpredictive*ones(M_ensemble,1);
Y = Hy*Xpredictive;
Ymean = Hy*Xmean;


% Kalman estimation of state
Zhat = Hy*Xpredictive - (Ymean*ones(1,M_ensemble));
Xhat = Xpredictive - (Xmean*ones(1,M_ensemble));
Mmatrix = (Xhat*Zhat')/(M_ensemble-1);
R = eye(dim_obs)*s2y;
S = ((Zhat*Zhat')/(M_ensemble-1))+R;



% if nosc>100
%     numberofblocks = floor(nosc/10);
% else
%      numberofblocks = 2;
% end
% S_eig = zeros(dim_obs,1);
% S_inv = zeros(dim_obs,dim_obs);
% 
% for i = 1:numberofblocks
%     part = ((i-1)*(dim_obs/numberofblocks))+1:i*dim_obs/numberofblocks;
%     aux = S(part,part);
%     S_inv(part,part) = inv(aux);
%     S_eig(part,1)=eig(aux);
% end

S_inv = inv(S);
S_eig = eig(S);

Kalman_Gain = Mmatrix*S_inv;
innov = y(1:K:nosc)'*ones(1,M_ensemble)-Y;

Xupdated = Xpredictive + Kalman_Gain*innov;
xmean_updated = (1/M_ensemble)*Xupdated*ones(M_ensemble,1);

% Likelihood
innov_mean = y(1:K:nosc)' - Hy*xmean_updated;
lW = -0.5 * ( sum(log(S_eig)) + innov_mean' * S_inv * innov_mean );


end 



