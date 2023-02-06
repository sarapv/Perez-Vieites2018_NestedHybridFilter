%% Initializations

nosc = 10;      % nosc = d_x = no. of oscillators or dimension of slow variables (x)
N = 500;        % no. particles of first layer
M = 500;        % no. particles of second layer
K = 2;          % we observe one every K slow oscillators, e.g., [x_1, x_3, ..., x_9]
t_final = 0.6;   % duration of the simulation in natural time units

total_iterations = 1;

%% Run

for iter = total_iterations

    [Output_FA] = NPF_Lorenz96(nosc,N,M,K,t_final,iter);
    

end