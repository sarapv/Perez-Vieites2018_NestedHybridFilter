function A = lsestimate(x,PSI)
%
% function A = lsestimate(x,PSI)
%

[nosc,NT] = size(x);

% autocorrelation matrix
Rx2 = zeros([2 2]);
Cx = zeros([2 1]);
for jj = 1:nosc
    x2 = [ x(jj,1:NT-1); x(jj,1:NT-1).^2 ];
    Rx2 = Rx2 + x2*x2'./NT;
    
    Cx = Cx + sum( ( ones([2 1])*PSI(jj,1:NT-1) ).*x2, 2)./NT;
end %jj

A = inv(Rx2)*Cx;


