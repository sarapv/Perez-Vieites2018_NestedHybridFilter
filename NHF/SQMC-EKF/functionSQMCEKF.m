function [ SQMC_EKF ] = functionSQMCEKF( nosc , M, gap, iteration )


rng(100*sum(clock)+100*iteration,'twister')


%% LORENZ 96 parameters
F = 8;                  % forcing parameter
H = 0.75; %0.75               % coupling between slow and fast variables
C = 10;                 % time scale of variables y
B = 15;%15;                 % inverse amplitude of the fast variables
h = 5e-3;               % integration period in natural units, antes era 2e-4
t_final = 40;           % duration of the simulation in natural time units
NT = fix(t_final/h);    % no. of discrete time steps
% nosc = 300;           % no. of oscillators
fps = 10;
nosc_fast = fps*nosc;   % no. of fast variables per slow variable 
s2z = 4;                % variance of the observations: slow variables
s2u = 1/10;             % variance of the observations: fast variables
s2x = h/4; %1/2, h/2;              % variance of the signal noise (slow variables)
s2y = h/16; %1              % variance of the signal noise (fast variables)
Tobs = fix(gap/h);             % signals observed every Tobs time steps
Oobs_x = 2;             % we observe one every Oobs_x slow oscillators
Oobs_y = 5;             % we observe one every Oobs_y fast oscillators

fprintf(1,'Parameters: F=%2.2f - H=%2.2f - B=%2.2f - C=%2.2f \n', F,H,B,C);

fprintf(1,'gap: %2.2f, nosc: %2.0f, M: %2.0f  \n', gap, nosc, M);
% Initial conditions
% x0 = zeros([nosc 1]);
% x0 = rand([nosc 1]);
% y0 = zeros([nosc_fast 1]);
% y0 = (1/(C*B))*rand([nosc_fast 1]) - 1/(2*C*B);


p=haltonset(1,'Skip',1e3,'Leap',1e2);
p = scramble(p,'RR2');
vx=net(p,nosc);
x0 = vx;

p=haltonset(1,'Skip',1e3,'Leap',1e2);
p = scramble(p,'RR2');
vy=net(p,nosc_fast);
y0 = (1/(C*B))*vy - 1/(2*C*B);

% Particle filter parameters
% M = 100;                % no. of particles in the outer PF
F_range = [2 30];
A_range = [0 0.2; 0 0.2];
FA_range = [F_range; A_range];
s2M = [20 20 1/5 20 0.04 0.04]./(M^1.5); %[20 20 1/5 20 0.04 0.04]./(M^1.5); % [F C B H A1 A2]


fprintf(1,'iteration= %d \n', iteration);
t0 = clock;

%%
[FAe_sqmc, Xe_sqmc, FAp_sqmc, MSE_sqmc, MSEpsi_sqmc, xref] = SQMCEKF_RK4(s2x,s2y,s2z,h,NT,Tobs,Oobs_x,FA_range,nosc,M,10*eye(nosc),x0,y0,F,C,B,H,fps,nosc_fast,s2M([1 5])); %TB TENDRE QUE GENERAR LA x, la z y la A AQUI
       
ttotal = etime(clock,t0)/60;
fprintf(1, 'time to finish = %7.3f min\n',ttotal);       
  
SQMC_EKF = struct('FAe',FAe_sqmc,'Xe',Xe_sqmc,'FAp',FAp_sqmc,'MSE',MSE_sqmc, 'MSEpsi', MSEpsi_sqmc,'xref',xref);  
mit = floor(NT/Tobs/2);
mit2 = floor(NT/Tobs/4*3);
fprintf(1,'SQMC+EKF (F unk., ansatz): MSE (all) = %7.3f, (last half) %7.3f, (last quarter) %7.3f \n', mean(MSE_sqmc(1+Tobs:Tobs:NT)), mean(MSE_sqmc(1+mit*Tobs:Tobs:NT)), mean(MSE_sqmc(1+mit2*Tobs:Tobs:NT)));
fprintf(1,'\n ------------------------------------------------------------------------------------------------- \n \n');

etiq_save = sprintf('/export/usuarios01/spvieites/private_html/PAPER/DATA/EKF_SQMC_FA_nosc%d_M%d_Tobs%d_iter%d.mat', nosc, M, Tobs, iteration);
save(etiq_save);


end