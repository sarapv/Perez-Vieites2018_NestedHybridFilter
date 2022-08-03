function Xnew = fLX_ansatz(X,F,PSI,nosc,N)
% fLX_ansatz integrates the slow state using both x and z variables
%   
% Inputs
%   X   : set of N previous slow state vectors
%   F   : model parameter
%   PSI     : contribution of fast variables (polynomial with 2 coefficients)
%   nosc    : dimension of the slow state
%   N   : number of particles
%
% Outputs
%   Xnew    : N new state vectors (matrix of size [nosc, N])
%

% initialisation
Xnew = zeros(nosc,N);

% 3 problematic cases: indexes 1, 2 and N
Xnew(1,:)= -X(nosc,:).*(X(nosc-1,:)-X(2,:))-X(1,:)+F-PSI(1,:);
Xnew(2,:)= -X(1,:).*(X(nosc,:)-X(3,:))-X(2,:)+F-PSI(2,:);
Xnew(nosc,:)=-X(nosc-1,:).*(X(nosc-2,:)-X(1,:))-X(nosc,:)+F-PSI(nosc,:);

% general case
for i=3:nosc-1
 Xnew(i,:)= -X(i-1,:).*(X(i-2,:)-X(i+1,:))-X(i,:)+F-PSI(i,:);
end


end