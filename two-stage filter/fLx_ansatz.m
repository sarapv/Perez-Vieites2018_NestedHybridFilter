function xnew = fLx_ansatz(x,F,PSI,nosc)
% fLx_ansatz integrates the slow state using previous x variables
% 
%   
% Inputs
%   x   : set of N previous slow state vectors
%   F   : model parameter
%   PSI     : contribution of fast variables (polynomial with 2 coefficients)
%   nosc    : dimension of the slow state
%
% Outputs
%   xnew    : new state vectors 



xnew = zeros(nosc,1);

% 3 problematic cases: 1,2,N
xnew(1)= -x(nosc).*(x(nosc-1)-x(2))-x(1)+F-PSI(1);
xnew(2)= -x(1).*(x(nosc)-x(3))-x(2)+F-PSI(2);
xnew(nosc)=-x(nosc-1).*(x(nosc-2)-x(1))-x(nosc)+F-PSI(nosc);

% general case
for i=3:nosc-1
 xnew(i)= -x(i-1).*(x(i-2)-x(i+1))-x(i)+F-PSI(i);
end


end