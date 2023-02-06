function [Param_est, xest, Param_part] = NPF_parameterlayer(s2x,s2y,h,NT,Tobs,K,y,par_range,nosc,M,N,s2M,x_ref,A_ref)
%
% function [Pe, xe, Pp] = nbf_fast_FA_a(s2x,s2z,h,NT,Tobs,Oobs_x,z,par_range,nosc,M,N,s2M,x_ref,A_ref)
%
% Inputs
% s2x   : signal variance
% s2y   : variance of observations
% h     : integration period
% NT    : no. of discrete-time steps
% Tobs  : one observations every Tobs time steps
% K     : one slow variable observed every Oobs_x ones 
% y     : observations  
% par_range     : support of the unknown parameter F
% nosc  : no. of slow variables 
% M     : no. of particles, outer filter
% N     : no. of particles, inner filter
% s2M   : jittering variance
%
% Outputs
% Param_est     : estimates of F, A1 and A2
% xest          : estimates of x
% Param_part    : last particles of parameters

%% Estimates
xest = zeros(nosc,NT); 
Fest = zeros(1,NT);
Aest = zeros(2,NT);

%% Weights
W = ones([M 1])/M;

%% Prior sampling in the parameter space (uniform prior)
Fparticles = par_range(1,1) + ( par_range(1,2) - par_range(1,1) ) .* rand(1,M);
Aparticles = zeros(2,M);
Aparticles(1,1:M) = par_range(2,1) + ( par_range(2,2) - par_range(2,1) ) .* rand(1,M);
Aparticles(2,1:M) = par_range(3,1) + ( par_range(3,2) - par_range(3,1) ) .* rand(1,M);

%% Prior sampling in the state space
Xl = cell([M 1]); 
for m = 1:M
    Xl{m} = rand(nosc,N); 
end %n

%% Time loop
for t = 1+Tobs:Tobs:NT
    
%     t0 = clock;
    
    % Jittering
    Fparticles = rnd_tgaussian(Fparticles,s2M(1)*ones(1,M),par_range(1,1),par_range(1,2));
    Aparticles(1,1:M) = rnd_tgaussian(Aparticles(1,1:M),s2M(2)*ones(1,M),par_range(2,1),par_range(2,2));
    Aparticles(2,1:M) = rnd_tgaussian(Aparticles(2,1:M),s2M(2)*ones(1,M),par_range(3,1),par_range(3,2));
    
    % Weights (one particle at a time)
    Xep = zeros(nosc,M); % state estimations of each PF of the second layer
    lWu = ones(M,1);
    
   % This loop can be parallelised using a parfor
    for m = 1:M
        % One BF step. It returns the log-weight (non normalised), the new
        % particles in the state space and a conditional estimate of the
        % states.
        [Xl{m}, aux, lWu(m)] = NPF_statelayer_1step(Fparticles(m),Aparticles(:,m),s2x,s2y,h,Tobs,K,y(:,t),Xl{m});
        Xep(:,m) =aux; 
    end %n
    
    % Weight normalisation
    Wu = exp( lWu - max(lWu) );
    W = Wu ./ sum(Wu);

    % State estimation
    xest(1:nosc,t) = Xep(:,1:M)*W;
    
    % Parameter estimation
    Fest(t) = Fparticles*W;
    Aest(:,t) = Aparticles*W;
    
    % Resampling
    idx = randsample(1:M,M,true,W);
    Fparticles_old = Fparticles;
    Fparticles = Fparticles(idx);
    Aparticles = Aparticles(:,idx);
    Xl_old = Xl;
    for m = 1:M
        Xl{m} = Xl_old{idx(m)};
    end %m
    
    % Figuras
    if rem(t,Tobs)==1
        figure(3);
        
        for kk = 1:min(nosc,12)
            subplot(4,4,kk);
            plot(h*(1:t),x_ref(kk,1:t),'k-');
            hold on;
            plot(h*(1+Tobs:Tobs:t),xest(kk,1+Tobs:Tobs:t),'g--');
            hold off;
            xlabel('time');
            etiq = sprintf('slow var. %d', kk);
            ylabel(etiq);
        end %kk
        
        subplot(4,4,13);
        plot(h*(1:t),8*ones([1 t]),'k-');
        hold on;
        plot(h*(1+Tobs:Tobs:t),Fest(1+Tobs:Tobs:t),'g-');
        plot(h*t*ones([1 M]),Fparticles(1:M),'ro');
        hold off;
        axis([0 h*(t+20) 0 25]); 
        xlabel('time');
        ylabel('F');

        subplot(4,4,14);
        semilogy(Fparticles_old,Wu,'bs');
        ylabel('n.n. weights');
        xlabel('F');
        axis([par_range(1,1) par_range(1,2) 1e-3 1]);
        
        subplot(4,4,15);
        plot(h*(1:t),A_ref(1)*ones([1 t]),'k-');
        hold on;
        plot(h*(1:t),A_ref(2)*ones([1 t]),'k--');
        plot(h*(1+Tobs:Tobs:t),Aest(1,1+Tobs:Tobs:t),'g-');
        plot(h*(1+Tobs:Tobs:t),Aest(2,1+Tobs:Tobs:t),'g--');
        hold off;
        xlabel('time');
        ylabel('ansatz');
        axis([0 h*(t+20) par_range(2,1) par_range(2,2)]);
        
        subplot(4,4,16);
        plot(Aparticles(1,1:M),Aparticles(2,1:M),'co');
        hold on;
        plot(A_ref(1),A_ref(2),'kx');
        hold off;
        axis([par_range(2,1) par_range(2,2) par_range(2,1) par_range(2,2)]);
        xlabel('a(1)'); 
        ylabel('a(2)');
        
    end %if
    
    if rem(t,10*Tobs)==1
        fprintf(1,'NPF step %d\n', t);
    end
    
end %t

%% Outputs
Param_est = [Fest; Aest];
Param_part = [Fparticles; Aparticles];