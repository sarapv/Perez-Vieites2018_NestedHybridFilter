function [ a ] = InvTransformMethod( u, W , N)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

a = zeros(N,1);

s = W(1); m = 1;

for n = 1 : N
    while (s < u(n) )
       m = m+1;
       s = s + W(m);
    end
    
    a(n) = m;
end

