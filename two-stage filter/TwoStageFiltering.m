function [ Output_2SF ] = TwoStageFiltering(nosc,N,M,K,t_final,iter)
%
%   [Output_2SF] = TwoStageFiltering(nosc,N,M,K,t_final,iter)
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
%   Output_2SF   : structure that contains Aest (estimates of param. A1 and
%    A2), Fest (estimates of param. F), Xest (estimates of x) and MSE



rng(100*sum(clock)+100*iter,'twister')

%% LORENZ 96 parameters
F = 8;                  % forcing parameter
H = 0.75;               % coupling between slow and fast variables
C = 10;                 % time scale of variables y
B = 15;                 % inverse amplitude of the fast variables

%time variables
h = 5e-3;               % integration period in natural units
NT = fix(t_final/h);    % no. of discrete time steps
gap = 0.05;
Tobs = fix(gap/h);      % signals observed every Tobs time steps

% fast variables (z)
fps = 10;               % no. of fast variables per slow variable 

% noise parameters
s2z = 4;                    % variance of the observations: slow variables
s2x = h/4; %1/2, h/2;       % variance of the signal noise (slow variables)
s2y = h/16; %1              % variance of the signal noise (fast variables)


% Particle filter parameters
F_range = [2 30];
A_range = [0 0.2; 0 0.2];
FA_range = [F_range; A_range];
s2M = [20 20 1/5 20 0.04 0.04]./(N^1.5); % [F C B H A1 A2] - jittering variables


fprintf(1,'Parameters: F=%2.2f - H=%2.2f - B=%2.2f - C=%2.2f \n', F,H,B,C);
fprintf(1,'gap: %2.2f, nosc: %2.0f, M: %2.0f  \n', gap, nosc, N);

%%
fprintf(1,'iteration= %d \n', iter);

t0 = clock;

[MSEx, Aest, Fest, Xest, xref]= SF_RK4(s2x,s2y,s2z,h,NT,Tobs,K,FA_range,N,M,F,C,B,H,nosc,fps,s2M([1 5])); 

ttotal = etime(clock,t0)/60;

fprintf(1, 'time to finish = %7.3f min\n',ttotal);
fprintf(1,'2 stage filter (F unk., ansatz): MSEx = %7.3f \n', mean(MSEx(1+Tobs:Tobs:NT)));
fprintf(1,'\n ------------------------------------------------------------------------------------------------- \n \n');

Output_2SF = struct('Aest',Aest,'Fest',Fest,'Xest',Xest,'MSE',MSEx,'xref',xref);
etiq_save = sprintf('data/2SF_FA_nosc%d_N%d_M%d_Tobs%d_iter%d.mat', nosc, N,M, Tobs, iter);
save(etiq_save);


end

