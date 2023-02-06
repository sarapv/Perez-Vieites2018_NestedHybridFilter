function [Output_FA] = NPF_Lorenz96(nosc,N,M,K,t_final,iter)
%
%   [Output_FA] = NPF_Lorenz96(nosc,N,M,K,t_final,iter)
%
% Inputs:
%   nosc    : no. of oscillators or dimension of slow variables (x)
%   N       : no. particles of first layer
%   M       : no. particles of second layer
%   K       : we observe one every K slow oscillators
%   t_final : duration of the simulation in natural time units
%   iter    : no. to label the experiment
%
% Outputs:
%   Output_FA   : structure that contains FAe (estimates of param. F A1
%   A2), Xe_FA (estimates of x) and FAp (particles of param. F A1 A2)


%% LORENZ 96 parameters
F = 8;                  % forcing parameter
H = 0.75;               % coupling between slow and fast variables
C = 10;                 % time scale of variables y
B = 15;                 % inverse amplitude of the fast variables

% time variables
h = 5e-3;               % integration period in natural units
NT = fix(t_final/h);    % no. of discrete time steps
Tobs = fix(0.05/h);     % signals observed every Tobs time steps

% fast variables (z)
fps = 10;
nosc_fast = fps*nosc;   % no. of fast variables per slow variable 

% noise parameters
s2y = 4;                % variance of the observations: slow variables
s2x = h/4;              % variance of the signal noise (slow variables)
s2z = h/16;              % variance of the signal noise (fast variables)


% Initial conditions
x = zeros(nosc, NT);
x(1:nosc,1) = rand(nosc,1);
z = zeros(nosc_fast,NT);
z(1:nosc_fast,1) = (1/(C*B))*rand(nosc_fast,1) - 1/(2*C*B);
x_old = x(1:nosc,1);
y_old = z(1:nosc_fast,1);
PSI = zeros(nosc,NT);

% Particle filter parameters
% N = 800;                % no. of particles in the inner PF
% M = 800;                % no. of particles in the outer PF
fprintf(1,'M=%d, N=%d\n',M,N);
F_range = [2 30];
A_range = [0 0.2; 0 0.2];
FA_range = [F_range; A_range];
s2M = [20 20 1/5 20 0.04 0.04]./(M^1.5); % [F C B H A1 A2] - jittering variance


%% RK4 integration
    t0 = clock;
    ok = 0;
    while not(ok)
        for n = 2:NT
                      
            [x1, PSI1] = fLX(x_old,y_old,F,C,B,H,nosc,fps);
            y1 = fLZ(x_old,y_old,F,C,B,H,nosc,fps);

            [x2, PSI2] = fLX(x_old+1/2.*h.*x1+1/2.*sqrt(s2x)*randn(size(x_old)),y_old+1/2.*h.*y1+1/2.*sqrt(s2z).*randn(size(y_old)),F,C,B,H,nosc,fps);
            y2 = fLZ(x_old+1/2.*h.*x1+1/2.*sqrt(s2x)*randn(size(x_old)),y_old+1/2.*h.*y1+1/2.*sqrt(s2z).*randn(size(y_old)),F,C,B,H,nosc,fps);

            [x3, PSI3]  = fLX(x_old+1/2.*h.*x2+1/2.*sqrt(s2x)*randn(size(x_old)),y_old+1/2.*h.*y2+1/2.*sqrt(s2z).*randn(size(y_old)),F,C,B,H,nosc,fps);
            y3 = fLZ(x_old+1/2.*h.*x2+1/2.*sqrt(s2x)*randn(size(x_old)),y_old+1/2.*h.*y2+1/2.*sqrt(s2z).*randn(size(y_old)),F,C,B,H,nosc,fps);

            [x4, PSI4] = fLX(x_old+h.*x3+sqrt(s2x)*randn(size(x_old)),y_old+h.*y3+sqrt(s2z).*randn(size(y_old)),F,C,B,H,nosc,fps);
            y4 = fLZ(x_old+h.*x3+sqrt(s2x)*randn(size(x_old)),y_old+h.*y3+sqrt(s2z).*randn(size(y_old)),F,C,B,H,nosc,fps);

            z(:,n) = y_old + (1/6)*(  h*( y1 + (2*y2) + (2*y3) + y4 ) + sqrt(s2z).*randn(size(y_old)) );%+ sqrt(s2y)*randn;
            x(:,n) = x_old + (1/6)*(  h*( x1 + (2*x2) + (2*x3) + x4 ) + sqrt(s2x).*randn(size(x_old)) );%+ sqrt(s2x)*randn;
            PSI(:,n) = 1/6 * (PSI1 + 2*PSI2 + 2*PSI3 + PSI4);
           
            y_old = z(:,n);
            x_old = x(:,n);
            
        end %n
        ok = isempty( find( isnan(x(1,:)) | isinf(x(1,:)), 1 ) );
    end %while

    % Observations
    y = x + sqrt(s2y)*randn(size(x));   % for the slow variables

    %% Ansatz (estimation of coefficients A1 and A2 with least squares)
    A = lsestimate(x,PSI);
    fprintf(1,'2nd order ansatz: %5.4f, %5.4f.\n', A(1), A(2));
    fprintf(1,'time: %5.4f s \n', etime(clock,t0));
    
     %% Nested BF, unknown F with unknown ansatz
    t0 = clock;
    
    [FAest, Xest_FA, FApart] = NPF_parameterlayer(s2x,s2y,h,NT,Tobs,K,y,FA_range,nosc,M,N,s2M([1 5]),x,A); 
    MSE_FA_slow = sum( sum( (Xest_FA(:,1+Tobs:Tobs:NT) - x(:,1+Tobs:Tobs:NT)).^2 ) ) / sum( sum( x(:,1+Tobs:Tobs:NT).^2 ) );

    Output_FA = struct('FAest',FAest,'Xest',Xest_FA,'FApart',FApart,'MSE',MSE_FA_slow);
    ttotal = etime(clock,t0)/60;
    MSE = sum( (Xest_FA(:,1+Tobs:Tobs:NT) - x(:,1+Tobs:Tobs:NT)).^2 )/nosc;
    fprintf(1,'NB+PF (F unk., ansatz): MSE = %7.3f\n',mean(MSE));
    fprintf(1, 'time = %7.3f min\n',ttotal);
   
    
     %% Saves data
     etiq_save = sprintf('data/lorenz96fast_ansatz_NPF_FA_nosc%d_N%d_M%d_iter%d.mat',nosc,N,M,iter);
     save(etiq_save);

end