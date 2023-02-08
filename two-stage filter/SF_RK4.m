function [MSEx, Aest, Fest, Xest, x_ref] = SF_RK4(s2x,s2z,s2y,h,NT,Tobs,K,par_range,N,M,F,C,B,H,nosc,fps,s2M)
%
% function [MSE, A_estim, F_estim, X_estim, x_ref] = SF_RK4(s2x,s2z,s2y,h,NT,Tobs,K,par_range,N,M,F,C,B,H,nosc,fps,s2M)
%
% Inputs
% s2x        : slow state noise variance
% s2z        : fast state noise variance
% s2y        : variance of observations
% h          : integration period
% NT         : no. of discrete-time steps
% Tobs       : one observations every Tobs time steps
% K          : one slow variable observed every Oobs_x ones 
% par_range  : support of the unknown parameter F
% N          : no. of particles, first layer
% M          : no. of ensembles
% F, C, H, B : true parameters
% nosc       : no. of oscillators or dimension of slow variables (x)
% fps        : no. of fast variables (z) per slow variable (x)
% s2M        : jittering variance
%
% Outputs
% MSE          : state MSE
% Aest         : estimates of A1 and A2
% Fest         : estimates of F
% Xest         : estimates of x
% x_ref        : true values of x_{1,t} and x_{2,t}
%



% fast variables dimension
nosc_fast = nosc*fps; 

% dimension of observations
n_obs = length(1:K:nosc);

% observation matrix
Hz_full = eye(nosc);
Hz = Hz_full(1:K:nosc,:);
clearvars Hz_full 

%% Initial conditions
x_old = rand(nosc,1);
z_old = (1/(C*B))*rand(nosc_fast,1) - 1/(2*C*B);


%% Estimates
Xest = zeros(nosc,NT);
Fest = zeros(1,NT);
Aest = zeros(2,NT);

Xensemble = randn(nosc,M); 
Xest_aux = randn(nosc,1);
MSEx = zeros(1,NT);

%% Prior sampling in the parameter space (uniform prior)
Fparticles = par_range(1,1) + ( par_range(1,2) - par_range(1,1) ).*rand(1,N);
Fest(1) = mean(Fparticles);

Aparticles = zeros(2,N);
Aparticles(1,1:N) = (par_range(2,1) + ( par_range(2,2) - par_range(2,1) ).*rand(1,N));
Aparticles(2,1:N) = (par_range(3,1) + ( par_range(3,2) - par_range(3,1) ) .*rand(1,N));
Aest(:,1) = mean(Aparticles,2);


%% Para generar las observaciones 
x_ref=zeros([2 NT]); %save x_{1,t} and x_{2,t} true values
x_ref(:,1)=[x_old(1); x_old(2)];
 
for t = 2:NT   
    %% Generating grount truth
    
    % States
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
        
    x_ref(:,t)=[x(1); x(2)];
    
    % Observations
    y = x' + sqrt(s2y)*randn(size(x'));   % for the slow variables
    
 %% Filter   
    if (mod(t,Tobs)==1)
        
        % Jittering      
        Fparticles = rnd_tgaussian(Fparticles,s2M(1)*ones([1 N]),par_range(1,1),par_range(1,2));
        Aparticles(1,1:N) = rnd_tgaussian(Aparticles(1,1:N),s2M(2)*ones(1,N),par_range(2,1),par_range(2,2));
        Aparticles(2,1:N) = rnd_tgaussian(Aparticles(2,1:N),s2M(2)*ones(1,N),par_range(3,1),par_range(3,2));
             
        % Prediction (N particles)
        Xpredictive = zeros(nosc,N);
        Xprev = Xest_aux*ones(1,N);
        for ii = 1:Tobs
            for i = 1:N
                Faux = Fparticles(i);
                Aaux = Aparticles(:,i);

                PSI1 = Aaux'*[Xprev(:,i).'; (Xprev(:,i).').^2];
                k1 = fLx_ansatz(Xprev(:,i),Faux,PSI1,nosc);

                aux = Xprev(:,i)+1/2.*h.*k1;
                PSI2 = Aaux'*[(aux).'; ((aux).').^2];
                k2 = fLx_ansatz(aux,Faux,PSI2,nosc);

                aux = Xprev(:,i)+1/2.*h.*k2;
                PSI3 = Aaux'*[(aux).'; ((aux).').^2];
                k3 = fLx_ansatz(aux,Faux,PSI3,nosc);

                aux = Xprev(:,i)+h.*k3;
                PSI4 = Aaux'*[(aux).'; ((aux).').^2];
                k4 = fLx_ansatz(aux,Faux,PSI4,nosc);

                Xpredictive(:,i)= Xprev(:,i) + 1/6*( h*(k1+2*k2+2*k3+k4)) ;

            end
            
            Xprev = Xpredictive;
           
        end
       
        % Likelihood
        y_estim = Xpredictive + sqrt(s2y)*randn(nosc,N);    
        llk = -(1/(2*s2y)).*sum( ( y(1:K:nosc).'*ones([1 N]) - y_estim(1:K:nosc,1:N) ).^2 );

        % weights
        wu = exp( llk - max(llk) );
        w = wu./sum(wu);     
        if sum(isnan(w)) > 0
            w = (1/N).*ones(1,N);
        end
            
        % Compute point estimates
        Aest(:,t) = Aparticles*w.';
        Fest(t)= Fparticles*w.';
               
        
        % Resampling
        idx_mc = randsample(1:N,N,true,w);
        Fparticles_old = Fparticles;
        Aparticles_old = Aparticles;
        Fparticles = Fparticles_old(idx_mc);    
        Aparticles = Aparticles_old(:,idx_mc);

        
        % EnKF stage
        X0 = Xensemble;
        Aaux = Aest(:,t);

        for k=1:Tobs
            PSI1 = Aaux(1).*X0 + Aaux(2).*(X0.^2);
            k1 = fLX_EnKF_ansatz(X0,Fest(t),PSI1);

            aux = X0+1/2.*h.*k1;
            PSI2 = Aaux(1).*aux + Aaux(2).*(aux.^2);
            k2 = fLX_EnKF_ansatz(aux,Fest(t),PSI2);
            clearvars  aux 

            aux = X0+1/2.*h.*k2;
            PSI3 = Aaux(1).*aux + Aaux(2).*(aux.^2);
            k3 = fLX_EnKF_ansatz(aux,Fest(t),PSI3);
            clearvars  aux 

            aux = X0+h.*k3;
            PSI4 = Aaux(1).*aux + Aaux(2).*(aux.^2);
            k4 = fLX_EnKF_ansatz(aux,Fest(t),PSI4);
            clearvars  aux 

            Xensemble = X0 + 1/6*( h*(k1+2*k2+2*k3+k4));

            if isnan(Xensemble)
               fprintf(1,'nan'); 
            end
            
            X0 = Xensemble + 0.5*randn(nosc,M);

        end %t

        Xmean = (1/M)*Xensemble*ones(M,1);
        Y = Hz*Xensemble;% +  sqrt(s2z)*randn(size(n_obs,M));
        Ymean = (1/M)*Y*ones(M,1);


        % Kalman update
        Z_hat = Y - (Ymean*ones(1,M));
        X_hat = Xensemble - (Xmean*ones(1,M));
        Matrix = (X_hat*transpose(Z_hat))/(M-1);
        R = eye(n_obs)*s2y;
        S = (Z_hat*transpose(Z_hat)/(M-1))+R;

        KalmanGain = Matrix*inv(S);

        innov = y(1:K:nosc)'*ones(1,M)-(Y+sqrt(s2y)*randn(n_obs,M));
        Xupdated = Xensemble + KalmanGain*innov;
        Xensemble = Xupdated;
        Xest_aux = (1/M)*Xupdated*ones(M,1);
        
        
        if sum(isnan(Xest_aux))>1 || max(Xest_aux)>50 || min(Xest_aux)<-50
            Xest_aux = randn(nosc,1);
            Xensemble = randn(nosc,M);
            Fparticles = par_range(1,1) + ( par_range(1,2) - par_range(1,1) ).*rand(1,N);
            Aparticles(1,1:N) = (par_range(2,1) + ( par_range(2,2) - par_range(2,1) ).*rand(1,N));
            Aparticles(2,1:N) = (par_range(3,1) + ( par_range(3,2) - par_range(3,1) ) .*rand(1,N));
            fprintf(1,'Reset filter\n');
        end
        
        % Estimates and MSE
        Xest(:,t) = Xest_aux;
        MSEx(t) = real(sum((Xest_aux-x).^2)/nosc); 
      
        
        % Figures
        if rem(t,Tobs)==1
            
            figure(3);
            for k = 1:2
                subplot(2,4,k);
                plot(h*(1:t),x_ref(k,1:t),'k-');
                hold on;
                plot(h*(1+Tobs:Tobs:t),Xest(k,1+Tobs:Tobs:t),'g--');
                hold off;
                xlabel('time');
                etiq = sprintf('slow var. %d',k);
                ylabel(etiq);
            end %kk
            
            subplot(2,4,5);
            plot(h*(1:t),F*ones([1 t]),'k-');
            hold on;
            plot(h*(1+Tobs:Tobs:t),Fest(1+Tobs:Tobs:t),'g-');
            hold off;
            axis([0 h*(t+20) 0 25]); 
            xlabel('time');
            ylabel('F');

            subplot(2,4,7);
            %plot(h*(1:t),A_ref(1)*ones([1 t]),'k-');
            %plot(h*(1:t),A_ref(2)*ones([1 t]),'k--');
            plot(h*(1+Tobs:Tobs:t),Aest(1,1+Tobs:Tobs:t),'g-');
            hold on;
            plot(h*(1+Tobs:Tobs:t),Aest(2,1+Tobs:Tobs:t),'g--');
            hold off;
            xlabel('time');
            %ylabel('ansatz');
            axis([0 h*(t+20) par_range(2,1) par_range(2,2)]);

            subplot(2,4,8);
            plot(Aparticles_old(1,:),Aparticles_old(2,:),'co');
%             hold on;
%             plot(A_ref(1),A_ref(2),'kx');
            hold off;
            axis([par_range(2,1) par_range(2,2) par_range(2,1) par_range(2,2)]);
            xlabel('a(1)'), ylabel('a(2)')
             
            subplot(2,4,6);
            semilogy(Fparticles_old,w,'bs')
            ylabel('n.n. weights'), xlabel('F')
            axis([par_range(1,1) par_range(1,2) 1e-4 1]);
            
            subplot(2,4,4)
            semilogy(h*(1+Tobs:Tobs:t),MSEx(1,1+Tobs:Tobs:t),'b');   
            ylabel('MSE_x')
                   
            pause(0.1);
        end %if rem        
        
    end %if mod
    
    % Prepare true state variables for next loop
    x_old=x;
    z_old=z;
    
end %t

end
