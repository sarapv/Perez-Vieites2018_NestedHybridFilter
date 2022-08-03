function [Xparticles, Xest, lW] = NPF_statelayer_1step(F,A,s2x,s2y,h,Tobs,K,y,Xold)
% function [Xparticles, Xest, lW] =
% bf_fast_1step(Params,t,s2x,s2z,h,Tobs,Oobs_x,z,Xold) does one step in one
% of the inner PFs in the second layer of the NPF (state estimation layer)
%
% Inputs
% F, A  : two-scale Lorenz 96 with an ansatz for the fast vars. (A is a 2x1 vector)
% s2x   : process noise variances
% s2y   : observation noise variance for the slow variables
% h     : integration step
% Tobs  : observations are collected every Tobs time steps
% K     : one out every Oobs_x slow oscillators is observed
% y     : observations from the slow variables
% Xold  : particles for the slow variables
%
% Outputs
% Xparticles    : particles
% Xest  : estimates of the state
% lW    : log-weight
%

% recovers the no. of particles (N), the no. of slow oscillators (nosc) and
% the no. of fast oscillators (nosc_fast)
[nosc, N] = size(Xold);

%% initialisation
Xparticles = zeros(size(Xold));

%% sampling
for k = 1:Tobs
    
    PSI1 = A(1).*Xold + A(2).*(Xold).^2 ;
    k1 = fLX_ansatz(Xold,F*ones(1,N),PSI1,nosc,N);
    
    aux = (Xold+1/2.*h.*k1);
    PSI2 = A(1).*aux + A(2).*(aux).^2 ;
    k2 = fLX_ansatz(aux,F*ones(1,N),PSI2,nosc,N);
    
    aux = Xold+1/2.*h.*k2;
    PSI3 = A(1).*aux + A(2).*(aux).^2 ;
    k3 = fLX_ansatz(aux,F*ones(1,N),PSI3,nosc,N);
    
    aux = Xold+h.*k3;
    PSI4 = A(1).*aux + A(2).*(aux).^2 ;
    k4 = fLX_ansatz(aux,F*ones(1,N),PSI4,nosc,N);
    Xparticles = Xold + h*(1/6*(k1+2*k2+2*k3+k4)  );    

    Xold = Xparticles;
end %k

%% Likelihoods, weights & estimates

% likelihoods
llk = -(1/(2*s2y)).*sum( ( y(1:K:nosc)*ones(1,N) - Xparticles(1:K:nosc,1:N) ).^2 );

% catches errors: if a particle diverges it is given an extremely low
% weight and its value reset
if not( isempty( find( isnan(llk) | isinf(llk), 1 ) ) )
    idx_bad = find( isnan(llk) | isinf(llk) );
    llk( idx_bad ) = -1000;
    Xparticles(:,idx_bad) = rand(nosc,length(idx_bad)); 
end %if

% weights
wu = exp( llk - max(llk) );
w = wu./sum(wu);

% log outer weight
lW = log( mean( wu ) ) + max(llk); 

% state estimate
Xest = Xparticles*w';

% resampling
idx = randsample(1:N,N,true,w);
Xparticles = Xparticles(:,idx);
