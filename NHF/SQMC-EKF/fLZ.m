function znew = fLZ(x,z,F,C,B,H,nosc,fps)
% fLZ integrates the fast state using both x and z variables
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
%   znew    : new state vector


% Initialisations
nosc_fast = fps*nosc;
znew = zeros(nosc_fast,1);

% 3 problematic cases: 1,N-1,N
id_x = 1;
znew(1)=C*B*z(2)*(z(nosc_fast) - z(3)) - C*z(1) + (C*H/B)*x(id_x)+ F*C/B;
id_x = nosc;
znew(nosc_fast-1)=  C*B*z(nosc_fast)*(z(nosc_fast-2) - C*z(1)) - z(nosc_fast-1) + (C*H/B)*x(id_x) + F*C/B;
znew(nosc_fast)= C*B*z(1)*(z(nosc_fast-1) - z(2)) - C*z(nosc_fast) + (C*H/B)*x(id_x) + F*C/B;

% the general case
for j=2:nosc_fast-2
    id_x = 1 + floor(j/fps);
    znew(j)= C*B*z(j+1)*(z(j-1) - z(j+2)) - C*z(j)  + (C*H/B)*x(id_x) + F*C/B;
end


end