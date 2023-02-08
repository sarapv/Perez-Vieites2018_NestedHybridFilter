function [X] = rnd_tgaussian_un(m,s2,v,t0,t1)
%
%
% Sampling from a Gaussian with mean m and variance s2, truncated outside the interval [t0,t1], by rejection sampling.
%
% 'm' and 's2' are vectors of the same size, containing the means and variances, respectively.
% 't0' and 't1' are the bounds; they can be either vectors of the same size as 'm' and 's2' or scalars.
% 'v' is a quasirandom point set



% First step
Fa = 0.5.*(1+erf((t0-m)./(sqrt(s2).*sqrt(2))));
Fb = 0.5.*(1+erf((t1-m)./(sqrt(s2).*sqrt(2))));

% Draw new point
X = m + sqrt(2).*sqrt(s2).*erfinv(2.*(Fa+(v.'.*(Fb-Fa)))-1);

