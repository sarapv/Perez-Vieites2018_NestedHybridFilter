function [Xnew, PSI] = fLX(x,z,F,C,B,H,nosc,fps)
% fLX integrates the slow state using both x and z variables
%   
% Inputs
%   x   : previous slow state vector
%   z   : previous fast state vector
%   F   : model parameter
%   C   : model parameter
%   B   : model parameter
%   H   : model parameter
%   nosc    : dimension of the slow state
%   fps     : number of fast states per slow state
%
% Outputs
%   Xnew    : new state vector
%   PSI     : contribution of fast variables in this integration step 
%
%

% initialisation
Xnew = zeros(nosc,1);
PSI = zeros(nosc,1);

% 3 problematic cases: indexes 1, 2 and N
Xnew(1)= -x(nosc)*(x(nosc-1) - x(2)) - x(1) + F - ((H*C)/B)*sum(z(1:fps));
PSI(1) = ((H*C)/B)*sum(z(1:fps));

Xnew(2)= -x(1)*(x(nosc) - x(3)) - x(2) + F - ((H*C)/B)*sum(z(fps+(1:fps)));
PSI(2) = ((H*C)/B)*sum(z(fps+(1:fps)));

Xnew(nosc)=-x(nosc-1)*(x(nosc-2) - x(1)) - x(nosc) + F - ((H*C)/B)*sum(z((nosc-1)*fps+(1:fps)));
PSI(nosc) = ((H*C)/B)*sum(z((nosc-1)*fps+(1:fps)));


% general case
for i=3:nosc-1
 Xnew(i)= -x(i-1)*(x(i-2) - x(i+1)) - x(i) + F - ((H*C)/B)*sum(z((i-1)*fps+(1:fps)));
 PSI(i) = ((H*C)/B)*sum(z((i-1)*fps+(1:fps)));
end


end
