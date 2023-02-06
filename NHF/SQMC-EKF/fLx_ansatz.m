function xnew = fLx_ansatz(x0,F,PSI)
% fLx_ansatz integrates the slow state using the previous x variable and
% the ansatz (polynomial of degree 2 with coefficients A(1) and A(2))
%  
% Inputs
%   x0   :  slow state vector
%   F    : model parameter
%   PSI     : contribution of fast variables (polynomial with 2 coefficients)
%
% Outputs
%   xnew    : new state vectors
%

nosc = size(x0,1);
xnew = zeros(nosc,1);

% 3 problematic cases: 1,2,N
xnew(1)= -x0(nosc).*(x0(nosc-1)-x0(2))-x0(1)+F-PSI(1);
xnew(2)= -x0(1).*(x0(nosc)-x0(3))-x0(2)+F-PSI(2);
xnew(nosc)=-x0(nosc-1).*(x0(nosc-2)-x0(1))-x0(nosc)+F-PSI(nosc);

%general case
for i=3:nosc-1
 xnew(i)= -x0(i-1).*(x0(i-2)-x0(i+1))-x0(i)+F-PSI(i);
end


end