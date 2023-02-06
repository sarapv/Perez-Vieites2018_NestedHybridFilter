function [Param_est, xest, Param_particles, MSEx, MSE_ansatz, xref] = SMCEnKF_Lorenz96(s2x,s2z,s2y,h,NT,Tobs,K,par_range,N,M_ensemble,nosc,F,C,B,H,fps,s2M)
% [Param_est, xest, Param_particles, MSEx, MSE_ansatz, xref] = SMCEnKF_Lorenz96(s2x,s2z,s2y,h,NT,Tobs,K,par_range,N,M_ensemble,nosc,F,C,B,H,fps,s2M)
% Algorithm with SMC for parameter estimation and EnKFs for state estimation
%
% Inputs
% s2x   : signal variance (slow state)
% s2z   : signal variance (fast state)
% s2y   : variance of observations
% h     : integration period
% NT    : no. of discrete-time steps
% Tobs  : one observations every Tobs time steps
% K     : one slow variable observed every Oobs_x ones
% par_range     : support of the unknown parameter F
% N             : no. of particles, first layer
% M_ensemble    : no. of samples in the second layer (EnKF)
% nosc          : no. of slow variables 
% F, C, H, B    : parameters
% fps           : no. of fast variables per slow variable 
% s2M           : jittering variance

%
% Outputs
% Param_est     : estimates of F, A1 and A2
% xest          : estimates of x
% Param_particles    : samples of parameters at different moments
% MSEx               : MSE of the slow state (x)
% MSE_ansatz         : MSE of the ansatz (polynomial of degree 2)
% xref               : true signals of x_{1,t} and x_{2,t}
%


%% Initialisations
nosc_fast = nosc*fps;

x0 = rand(nosc,1);
x_old = x0;
z0 = (1/(C*B))*rand(nosc_fast,1) - 1/(2*C*B);
z_old = z0;

saveparticles = round(NT/20);

% Estimates
Fest = zeros(1,NT);
Aest = zeros(2,NT);
Xparticles_layer1 = zeros(nosc,N);
Param_particles = [];
PSIparticles = zeros(nosc,N);
xest = zeros(nosc,NT);
MSEx = zeros(1,N);
MSE_ansatz = zeros(1,NT);

%% Covariance matrices and mean vectors
for n = 1:N
    Xensemble{n} = rand(nosc,M_ensemble);
    Xparticles_layer1(:,n) = mean(Xensemble{n},2);
end %n

%% Prior sampling in the parameter space (uniform prior)
Fparticles = par_range(1,1) + ( par_range(1,2) - par_range(1,1) ) .* rand(1,N);
Aparticles = zeros(2,N);
Aparticles(1,1:N) = par_range(2,1) + ( par_range(2,2) - par_range(2,1) ) .* rand(1,N);
Aparticles(2,1:N) = par_range(3,1) + ( par_range(3,2) - par_range(3,1) ) .* rand(1,N);

%% Generate ground truth and run algorithm
xref=zeros(2,NT); % save true values of x_{1,t} and x_{2,t}
xref(:,1)=[x_old(1); x_old(2)];

for t = 2:NT
    
    % Ground truth
    [x1, PSI1] = fLX(x_old,z_old,F,C,B,H,nosc,fps);
	z1 = fLZ(x_old,z_old,F,C,B,H,nosc,fps);
            
    [x2, PSI2] = fLX(x_old+1/2.*h.*x1+1/2.*sqrt(s2x)*randn(size(x_old)),z_old+1/2.*h.*z1+1/2.*sqrt(s2z).*randn(size(z_old)),F,C,B,H,nosc,fps);
    z2 = fLZ(x_old+1/2.*h.*x1+1/2.*sqrt(s2x)*randn(size(x_old)),z_old+1/2.*h.*z1+1/2.*sqrt(s2z).*randn(size(z_old)),F,C,B,H,nosc,fps);
            
    [x3, PSI3]  = fLX(x_old+1/2.*h.*x2+1/2.*sqrt(s2x)*randn(size(x_old)),z_old+1/2.*h.*z2+1/2.*sqrt(s2z).*randn(size(z_old)),F,C,B,H,nosc,fps);
	z3 = fLZ(x_old+1/2.*h.*x2+1/2.*sqrt(s2x)*randn(size(x_old)),z_old+1/2.*h.*z2+1/2.*sqrt(s2z).*randn(size(z_old)),F,C,B,H,nosc,fps);
            
    [x4, PSI4] = fLX(x_old+h.*x3+sqrt(s2x)*randn(size(x_old)),z_old+h.*z3+sqrt(s2z).*randn(size(z_old)),F,C,B,H,nosc,fps);
    z4 = fLZ(x_old+h.*x3+sqrt(s2x)*randn(size(x_old)),z_old+h.*z3+sqrt(s2z).*randn(size(z_old)),F,C,B,H,nosc,fps);
            
    z = z_old + (1/6)*(  h*( z1 + (2*z2) + (2*z3) + z4 ) + sqrt(s2z).*randn(size(z_old)) );%+ sqrt(s2y)*randn;
    x = x_old + (1/6)*(  h*( x1 + (2*x2) + (2*x3) + x4 ) + sqrt(s2x).*randn(size(x_old)) );%+ sqrt(s2x)*randn;
    PSItruth = (1/6)*(( PSI1 + (2*PSI2) + (2*PSI3) + PSI4 ));
    
    xref(:,t)=[x(1); x(2)];
    

    % Observations
    y = x' + sqrt(s2y)*randn(size(x'));   % for the slow variables
    
    % Run algorithm
    if (mod(t,Tobs)==1)
        
        % Jittering
        Fparticles = rnd_tgaussian(Fparticles,s2M(1)*ones(1,N),par_range(1,1),par_range(1,2));
        Aparticles(1,1:N) = rnd_tgaussian(Aparticles(1,1:N),s2M(2)*ones([1 N]),par_range(2,1),par_range(2,2));
        Aparticles(2,1:N) = rnd_tgaussian(Aparticles(2,1:N),s2M(2)*ones([1 N]),par_range(3,1),par_range(3,2));
        
        % quasirandom point set to generate noise later - needed when using SQMC
        %p=haltonset(nosc,'Skip',1e3,'Leap',1e2);
        %p = scramble(p,'RR2');
        %QMC_noise=net(p,M_ensemble*4); 
        
        % random point set to generate noise later - when using SMC
        just_noise = rand(M_ensemble*4,M_ensemble);


        %SECOND LAYER
        % Weights (one particle at a time)
        llw = ones(N,1);
        % this loop can be changed to a parfor in order to parallelise
        for n = 1:N %cambiar a parfor
            % One EnKF step. It returns the log-weight (non normalised), the new
            % covariance matrix and the new mean vector in the state space
            [Xaux, xmean, llw_aux, psi] = EnKFstep_fast_RK4_SQMC(Fparticles(n),Aparticles(:,n),s2x,s2y,h,Tobs,K,y,Xensemble{n},just_noise);                       
            Xensemble{n}=Xaux;
            Xparticles_layer1(:,n)=xmean;
            llw(n)=llw_aux;
            PSIparticles(:,n) = psi;

        end %n
        clearvars Xaux xmean llw_aux psi              
        
        % Weight normalisation
        wu = exp( llw - max(llw) );
        w = wu ./ sum(wu);

        %State and parameter estimation
        xest(:,t) = Xparticles_layer1(:,1:N)*w;
        Fest(t) = Fparticles*w;
        Aest(:,t) = Aparticles*w;
        PSIest = PSIparticles*w;
        
        %MSE
        MSEx(t) = sum((xest(:,t)-x).^2)/nosc;     
        MSE_ansatz(t) = sum((PSItruth - PSIest).^2)/nosc;

        % Resampling
        idx = randsample(1:N,N,true,w);
        Fparticles_old = Fparticles;
        Fparticles = Fparticles(idx);
        Aparticles = Aparticles(:,idx);
        Xparticles_layer1 = Xparticles_layer1(:,idx);
        Xensemble_old = Xensemble;
        
        for n = 1:N
            Xensemble{n} = Xensemble_old{idx(n)};
        end %n

        
        
        % Figures
        if rem(t,Tobs)==1
            
            figure(3);
            for kk = 1:2
                subplot(2,4,kk);
                plot(h*(1:t),xref(kk,1:t),'k-');
                hold on;
                plot(h*(1+Tobs:Tobs:t),xest(kk,1+Tobs:Tobs:t),'g--');
                hold off;
                xlabel('time');
                etiq = sprintf('SMC slow var. %d', kk);
                ylabel(etiq);
            end %kk
            
            subplot(2,4,5);
            plot(h*(1:t),8*ones(1,t),'k-');
            hold on;
            plot(h*(1+Tobs:Tobs:t),Fest(1+Tobs:Tobs:t),'g-');
            plot(h*t*ones(1,N),Fparticles(1:N),'ro');
            hold off;
            axis([0 h*(t+20) 0 25]); 
            xlabel('time'), ylabel('F');

            subplot(2,4,6);
            semilogy(Fparticles_old,wu,'bs');
            ylabel('n.n. weights');
            xlabel('F');
            axis([par_range(1,1) par_range(1,2) 1e-3 1]);

            subplot(2,4,7);
            %plot(h*(1:t),A_ref(1)*ones([1 t]),'k-');
            %plot(h*(1:t),A_ref(2)*ones([1 t]),'k--');
            plot(h*(1+Tobs:Tobs:t),Aest(1,1+Tobs:Tobs:t),'g-');
            hold on;
            plot(h*(1+Tobs:Tobs:t),Aest(2,1+Tobs:Tobs:t),'g--');
            hold off;
            xlabel('time');
            ylabel('ansatz');
            axis([0 h*(t+20) par_range(2,1) par_range(2,2)]);

            subplot(2,4,8);
            plot(Aparticles(1,1:N),Aparticles(2,1:N),'co');
            hold on;
            %plot(A_ref(1),A_ref(2),'kx');
            hold off;
            axis([par_range(2,1) par_range(2,2) par_range(2,1) par_range(2,2)]);
            xlabel('a(1)'), ylabel('a(2)');
           
            subplot(2,4,4)
            semilogy(h*(1+Tobs:Tobs:t),MSEx(1,1+Tobs:Tobs:t),'b');
            xlabel('MSE_x')
            
             
            pause(0.1);
            
        end %if rem
        
        
        if mod(t,saveparticles) == 1
            Param_particles = [Param_particles; Fparticles; Aparticles];
        end
        
    end %if mod
    
    x_old=x;
    z_old=z;
    
end %t

%% Outputs
Param_est = [Fest; Aest];

end