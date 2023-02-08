function [  Output_SQMCEnKF ] = NHF_SQMCEnKF_Lorenz96( nosc, M, gap, t_final, iter )
%
%   [Output_SQMCEnKF] = NHF_SQMCEnKF_Lorenz96( nosc, M, gap, t_final, iter )
%
% Inputs:
%   nosc    : no. of oscillators or dimension of slow variables (x)
%   M       : no. particles of first layer
%   gap     : time (in continuous-time units) between consecutive
%   observations
%   t_final : duration of the simulation in natural time units
%   iter    : no. to label the experiment
%
% Outputs:
%   Output_SQMCEKF   : structure that contains 
%                       - FAest (estimates of param. F A1 and A2), 
%                       - Xest (estimates of x),
%                       - FAparticles (particles of param. F A1 and A2),
%                       - MSEx (MSE of state x),
%                       - MSE_ansatz, and
%                       - xref (true values of state variables x_{1,t} and
%                       x_{2,t}

rng(100*sum(clock)+100*iter,'twister');

%% LORENZ 96 parameters
F = 8;                  % forcing parameter
H = 0.75;               % coupling between slow and fast variables
C = 10;                 % time scale of variables y
B = 15;                 % inverse amplitude of the fast variables

h = 5e-3;               % integration period in natural units, antes era 2e-4
% t_final = 40;           % duration of the simulation in natural time units
Tobs = fix(gap/h);             % signals observed every Tobs time steps
NT = fix(t_final/h);    % no. of discrete time steps

fps = 10;               % no. of fast variables per slow variable 
nosc_fast = fps*nosc;   % no. of fast variables 

s2y = 4;                % variance of the observations: slow variables
s2x = h/4;              % variance of the signal noise (slow variables)
s2z = h/16;             % variance of the signal noise (fast variables)

K = 2;             % we observe one every Oobs_x slow oscillators

%% Initializations

% Particle filter parameters
F_range = [2 30];
A_range = [0 0.2; 0 0.2];
FA_range = [F_range; A_range];
s2M = [20 20 1/5 20 0.04 0.04]./(M^1.5); %[20 20 1/5 20 0.04 0.04]./(M^1.5); % [F C B H A1 A2]



%% Run algorithm

fprintf(1,'Parameters: F=%2.2f - H=%2.2f - B=%2.2f - C=%2.2f \n', F,H,B,C);
fprintf(1,'gap: %2.2f, nosc: %2.0f, M: %2.0f  \n', gap, nosc, M);
fprintf(1,'iteration= %d \n', iter);
t0 = clock;


[FAest_sqmc, Xest_sqmc, FAparticles_sqmc, MSEx_sqmc, MSE_ansatz, xref]= SQMCEnKF_L96(s2x,s2z,s2y,h,NT,Tobs,K,FA_range,s2M([1 5]),M,F,C,B,H,fps,nosc); %TB TENDRE QUE GENERAR LA x, la z y la A AQUI

ttotal = etime(clock,t0)/60;

fprintf(1, 'time to finish = %7.3f min\n',ttotal);
fprintf(1,'SQMC+EnKF (F unk., ansatz): MSEx  = %7.3f\n', mean(MSEx_sqmc(1+Tobs:Tobs:NT)));
fprintf(1,'\n ---------------------------------------------------- \n \n');
% 
% etiq_save = sprintf('/export/usuarios01/spvieites/private_html/PAPER/DATA/EnKF_SQMC_FA_nosc%d_M%d_Tobs%d_iter%d.mat', nosc, M, Tobs, iteration);

% Save data
Output_SQMCEnKF = struct('FAest',FAest_sqmc,'Xest',Xest_sqmc,'FAparticles',FAparticles_sqmc,'MSEx',MSEx_sqmc, 'MSE_ansatz', MSE_ansatz,'xref',xref);  
etiq_save = sprintf('data/SQMCEnKF_FA_nosc%d_M%d_Tobs%d_iter%d.mat', nosc, M, Tobs, iter);
save(etiq_save);


end