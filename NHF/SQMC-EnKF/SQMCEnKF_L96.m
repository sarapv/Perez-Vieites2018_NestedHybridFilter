function [ Param_est_sqmc, xest_sqmc, Param_particles_sqmc, MSEx_sqmc, MSE_ansatz, xref] = SQMCEnKF_L96(s2x,s2z,s2y,h,NT,Tobs,K,par_range,s2M,M,F,C,B,H,fps,nosc)
       
% function[Param_est, xest, Param_particles, MSEx, MSE_ansatz, x_ref] = = SQMCEnKF_L96(s2x,s2y,s2z,h,NT,Tobs,K,par_range,s2M,M,F,C,B,H,fps,nosc)
% Algorithm with SQMC for parameter estimation and EKFs for state estimation
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
% s2M           : jittering variance
% M             : no. of particles, first layer
% F, C, H, B    : parameters
% fps           : no. of fast variables per slow variable 
% nosc          : no. of slow variables 
%
% Outputs
% Param_est     : estimates of F, A1 and A2
% xest          : estimates of x
% Param_particles    : samples of parameters at different moments
% MSEx               : MSE of the slow state (x)
% MSE_ansatz         : MSE of the ansatz (polynomial of degree 2)
% xref              : true signals of x_{1,t} and x_{2,t}
%

%% Initializations
nosc_fast = fps*nosc;
M_ensembles = nosc;     % no. of particles or ensemble members of the EnKFs

% Slow state x0 
p=haltonset(1,'Skip',1e3,'Leap',1e2);
p = scramble(p,'RR2');
x0=net(p,nosc);

% Initial ensemble
p=haltonset(M_ensembles,'Skip',1e3,'Leap',1e2);
p = scramble(p,'RR2');
X_ensembles=net(p,nosc);

% Fast state z0
p=haltonset(1,'Skip',1e3,'Leap',1e2);
p = scramble(p,'RR2');
vz=net(p,nosc_fast);
z0 = (1/(C*B))*vz - 1/(2*C*B);

x_old=x0;
z_old=z0;

saveparticles = round(NT/20);


% Estimates
%xe = zeros([nosc NT]); %COMENTO ESTO PORQUE SOLO VOY A CONSIDERAR 2 OSCILADORES PARA PINTAR Y LO PONGO MAS ABAJO
Fest_sqmc = zeros([1 NT]);
Aest_sqmc = zeros([2 NT]);
xest_sqmc = zeros([nosc M]);

% Error
MSEx_sqmc = zeros([1 NT]);
MSE_ansatz = zeros([1 NT]);

% Particles
Xparticles_layer1_sqmc = zeros(nosc,M);
Param_particles_sqmc = [];
PSIparticles_sqmc = zeros([nosc M]);

%% Ensembles (2nd layer)
for m = 1:M
    p=haltonset(M_ensembles,'Skip',1e3,'Leap',1e2);
    p = scramble(p,'RR2');
    vx=net(p,nosc);
    
    Xensemble_layer2_sqmc{m} = vx;
    Xparticles_layer1_sqmc(:,m) = mean(vx,2);
end %m


%% Prior sampling in the parameter space (uniform prior)

% QMC noise
p=haltonset(3,'Skip',1e3,'Leap',1e2);
p = scramble(p,'RR2');
v_tao=net(p,M);  

% Prior parameters
Fparticles_sqmc = par_range(1,1) + ( par_range(1,2) - par_range(1,1) ) .* v_tao(:,1).';
Aparticles_sqmc = zeros([2 M]);
Aparticles_sqmc(1,1:M) = par_range(2,1) + ( par_range(2,2) - par_range(2,1) ) .* v_tao(:,2);
Aparticles_sqmc(2,1:M) = par_range(3,1) + ( par_range(3,2) - par_range(3,1) ) .* v_tao(:,3);

Fparticles_sqmc = rnd_tgaussian_un(Fparticles_sqmc,par_range(1,2),v_tao(:,1),par_range(1,1),par_range(1,2));
Aparticles_sqmc(1,1:M) = rnd_tgaussian_un(Aparticles_sqmc(1,1:M),par_range(2,2),v_tao(:,2),par_range(2,1),par_range(2,2));
Aparticles_sqmc(2,1:M) = rnd_tgaussian_un(Aparticles_sqmc(2,1:M),par_range(3,2),v_tao(:,3),par_range(3,1),par_range(3,2));


%% Generate ground truth and run algorithm
xref=zeros([2 NT]); % save true values of x_{1,t} and x_{2,t}
xref(:,1)=[x_old(1); x_old(2)];

for t = 2:NT
    
    % Ground truth (Runge-Kutta 4)
    [x1, PSI1] = fLX(x_old,z_old,F,C,B,H,nosc,fps);
	z1 = fLZ(x_old,z_old,F,C,B,H,nosc,fps);
            
    [x2, PSI2] = fLX(x_old+1/2.*h.*x1+1/2.*sqrt(s2x)*randn(size(x_old)),z_old+1/2.*h.*z1+1/2.*sqrt(s2z).*randn(size(z_old)),F,C,B,H,nosc,fps);
    z2 = fLZ(x_old+1/2.*h.*x1+1/2.*sqrt(s2x)*randn(size(x_old)),z_old+1/2.*h.*z1+1/2.*sqrt(s2z).*randn(size(z_old)),F,C,B,H,nosc,fps);
            
    [x3, PSI3]  = fLX(x_old+1/2.*h.*x2+1/2.*sqrt(s2x)*randn(size(x_old)),z_old+1/2.*h.*z2+1/2.*sqrt(s2z).*randn(size(z_old)),F,C,B,H,nosc,fps);
	z3 = fLZ(x_old+1/2.*h.*x2+1/2.*sqrt(s2x)*randn(size(x_old)),z_old+1/2.*h.*z2+1/2.*sqrt(s2z).*randn(size(z_old)),F,C,B,H,nosc,fps);
            
    [x4, PSI4] = fLX(x_old+h.*x3+sqrt(s2x)*randn(size(x_old)),z_old+h.*z3+sqrt(s2z).*randn(size(z_old)),F,C,B,H,nosc,fps);
    z4 = fLZ(x_old+h.*x3+sqrt(s2x)*randn(size(x_old)),z_old+h.*z3+sqrt(s2z).*randn(size(z_old)),F,C,B,H,nosc,fps);
            
    z = z_old + (1/6)*(  h*( z1 + (2*z2) + (2*z3) + z4 ) + sqrt(s2z).*randn(size(z_old)) );%+ sqrt(s2z)*randn;
    x = x_old + (1/6)*(  h*( x1 + (2*x2) + (2*x3) + x4 ) + sqrt(s2x).*randn(size(x_old)) );%+ sqrt(s2x)*randn;
    PSIreal = (1/6)*(( PSI1 + (2*PSI2) + (2*PSI3) + PSI4 ));
    
    xref(:,t)=[x(1); x(2)];
    

    % Observations
    y = x' + sqrt(s2y)*randn(size(x'));   % for the slow variables
    
    
    % Run algorithm
    if (mod(t,Tobs)==1)
  
        % QMC noise
        p=haltonset(nosc,'Skip',1e3,'Leap',1e2);
        p = scramble(p,'RR2');
        noise=net(p,M_ensembles*4);     
        
         %%%% 2nd layer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        % Weights (one particle at a time)
        llk_sqmc = ones([M 1]); %log-likelihood
        for m = 1:M %can be changed for a parfor
            [Xensemble_aux, xaux, llk_aux, psi] = EnKFstep_fast_L96_SQMC(Fparticles_sqmc(m),Aparticles_sqmc(:,m),s2x,s2y,h,Tobs,K,y,Xensemble_layer2_sqmc{m}, noise);
            Xensemble_layer2_sqmc{m}=Xensemble_aux;
            Xparticles_layer1_sqmc(:,m)=xaux;
            llk_sqmc(m)=llk_aux;
            PSIparticles_sqmc(:,m) = psi;     
        end %m
        clearvars Xensemble_aux xaux llk_aux psi noise
        
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Weight normalisation
        wu_sqmc = exp( llk_sqmc - max(llk_sqmc) );
        w_sqmc = wu_sqmc ./ sum(wu_sqmc);

        % Parameter and state estimation
        xest_sqmc(:,t) = Xparticles_layer1_sqmc(:,1:M)*w_sqmc;
        Fest_sqmc(t) = Fparticles_sqmc*w_sqmc;
        Aest_sqmc(:,t) = Aparticles_sqmc*w_sqmc;
        PSIest_sqmc = PSIparticles_sqmc*w_sqmc;
               
        % MSE
        MSEx_sqmc(t) = sum((xest_sqmc(:,t)-x).^2)/nosc; 
        MSE_ansatz(t) = sum((PSIreal - PSIest_sqmc).^2)/nosc;
        
        % Parameter estimation
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Resampling
        % QMC
        p=haltonset(3 + 1,'Skip',1e3,'Leap',1e2);
        p = scramble(p,'RR2');
        vu=net(p,M);  

        % Function phi
        p_sup = mean( [Fparticles_sqmc.' Aparticles_sqmc.']) + 2.* std( [Fparticles_sqmc.' Aparticles_sqmc.']);
        p_inf = mean( [Fparticles_sqmc.' Aparticles_sqmc.']) - 2.* std( [Fparticles_sqmc.' Aparticles_sqmc.']);
        phi = (1 + exp( - ([Fparticles_sqmc.' Aparticles_sqmc.'] - ones(M,1)*p_inf)./(ones(M,1)*p_sup - ones(M,1)*p_inf) ) ).^(-1);
               
        % Function sort
        m_hilbert = 256;
        dist = zeros(M,1);
        phi_norm = floor(phi*(m_hilbert-1)/max(max(phi)));
        for pp = 1:M
           dist(pp) = xyz2d(m_hilbert, phi_norm(pp,1), phi_norm(pp,2), phi_norm(pp,3)); 
        end
        [~, I_sigma] = sort(dist,'ascend');
        

        % Sort X and w
        Aparticles_sigma = Aparticles_sqmc(:,I_sigma);
        Fparticles_sigma = Fparticles_sqmc(I_sigma);
        W_sigma = w_sqmc(I_sigma);
        Xensemble_layer1_sigma = Xparticles_layer1_sqmc(:,I_sigma);
        Xensemble_layer2_old_sqmc = Xensemble_layer2_sqmc;
        for m = 1:M
            Xensemble_layer2_sigma{m} = Xensemble_layer2_old_sqmc{I_sigma(m)};
        end %m

        % Sort u
        [u_tao, Iu_tao] = sort(vu(:,1),'ascend');

        Alabel = InvTransformMethod(u_tao,W_sigma,M);
        Aparticles_sqmc = Aparticles_sigma(:,Alabel);
        Fparticles_old_sqmc = Fparticles_sigma;
        Fparticles_sqmc = Fparticles_sigma(Alabel);
        Xparticles_layer1_sqmc = Xensemble_layer1_sigma(:,Alabel);
        Xensemble_layer2_old_sqmc = Xensemble_layer2_sigma;
        for m = 1:M
            Xensemble_layer2_sqmc{m} = Xensemble_layer2_old_sqmc{Alabel(m)};
        end %m

        % Jittering
        v = vu(:,2:end);
        v_tao = v(Iu_tao,:);
 
       
        % Jittering
        Fparticles_sqmc = rnd_tgaussian_un(Fparticles_sqmc,s2M(1),v_tao(:,1),par_range(1,1),par_range(1,2));
        Aparticles_sqmc(1,1:M) = rnd_tgaussian_un(Aparticles_sqmc(1,1:M),s2M(2),v_tao(:,2),par_range(2,1),par_range(2,2));
        Aparticles_sqmc(2,1:M) = rnd_tgaussian_un(Aparticles_sqmc(2,1:M),s2M(2),v_tao(:,3),par_range(3,1),par_range(3,2));
     

        
       
        % Figures
        if rem(t,Tobs)==1
            figure(3);

            for kk = 1:2
                subplot(2,4,kk);
                plot(h*(1:t),xref(kk,1:t),'k-');
                hold on;
                plot(h*(1+Tobs:Tobs:t),xest_sqmc(kk,1+Tobs:Tobs:t),'g--');
                hold off;
                xlabel('time');
                etiq = sprintf('SQMC slow var. %d', kk);
                ylabel(etiq);
            end %kk

            
            subplot(2,4,5);
            plot(h*(1:t),8*ones([1 t]),'k-');
            hold on;
            plot(h*(1+Tobs:Tobs:t),Fest_sqmc(1+Tobs:Tobs:t),'g-');
            plot(h*t*ones([1 M]),Fparticles_sqmc(1:M),'ro');
            hold off;
            axis([0 h*(t+20) 0 25]); 
            xlabel('time');
            ylabel('F');

            subplot(2,4,6);
            semilogy(Fparticles_old_sqmc,w_sqmc,'bs');
            ylabel('n.n. weights');
            xlabel('F');
            axis([par_range(1,1) par_range(1,2) 1e-3 1]);

            subplot(2,4,7);
            %plot(h*(1:t),A_ref(1)*ones([1 t]),'k-');
            %plot(h*(1:t),A_ref(2)*ones([1 t]),'k--');
            plot(h*(1+Tobs:Tobs:t),Aest_sqmc(1,1+Tobs:Tobs:t),'g-');
            hold on;
            plot(h*(1+Tobs:Tobs:t),Aest_sqmc(2,1+Tobs:Tobs:t),'g--');
            hold off;
            xlabel('time');
            ylabel('ansatz');
            axis([0 h*(t+20) par_range(2,1) par_range(2,2)]);

            subplot(2,4,8);
            plot(Aparticles_sqmc(1,1:M),Aparticles_sqmc(2,1:M),'co');
            hold on;
            %plot(A_ref(1),A_ref(2),'kx');
            hold off;
            axis([par_range(2,1) par_range(2,2) par_range(2,1) par_range(2,2)]);
            xlabel('a(1)'); 
            ylabel('a(2)');
             
            subplot(2,4,4)
            semilogy(h*(1+Tobs:Tobs:t),MSEx_sqmc(1,1+Tobs:Tobs:t),'r');
            xlabel('MSE_x')
             
            pause(0.1);
        end %if rem
        
        
        if mod(t,saveparticles) == 1
            Param_particles_sqmc = [Param_particles_sqmc; Fparticles_sqmc; Aparticles_sqmc];
        end
    end %if mod

    x_old=x;
    z_old=z;
    
end %t

%% Outputs
Param_est_sqmc = [Fest_sqmc; Aest_sqmc];


end