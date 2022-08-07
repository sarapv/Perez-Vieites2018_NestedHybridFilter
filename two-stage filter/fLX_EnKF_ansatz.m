function Xnew = fLX_EnKF_ansatz(X,F,PSI)
% fLEnKF integrates the slow state using previous x variables 
%   
% Inputs
%   X   : set of samples of previous slow state vectors
%   F   : model parameter
%   PSI     : contribution of fast variables (polynomial with 2 coefficients)
%
% Outputs
%   Xnew    : N new state vectors (matrix of size [nosc, N])
%

[nosc, n_ensembles] = size(X);
Xnew = zeros(nosc,n_ensembles);
Fvector = F*ones(1,n_ensembles);

% 3 problematic cases: 1,2,N
Xnew(1,:)= -X(nosc,:).*(X(nosc-1,:)-X(2,:))-X(1,:)+Fvector-PSI(1,:);
Xnew(2,:)= -X(1,:).*(X(nosc,:)-X(3,:))-X(2,:)+Fvector-PSI(2,:);
Xnew(nosc,:)=-X(nosc-1,:).*(X(nosc-2,:)-X(1,:))-X(nosc,:)+Fvector-PSI(nosc,:);

% general case
for i=3:nosc-1
 Xnew(i,:)= -X(i-1,:).*(X(i-2,:)-X(i+1,:))-X(i,:)+Fvector-PSI(i,:);
end


end