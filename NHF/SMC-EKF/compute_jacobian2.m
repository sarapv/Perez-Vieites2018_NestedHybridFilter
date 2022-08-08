function J = compute_jacobian2(X,A,h)
% J = compute_jacobian(X,A)
% X : linearisation point
% A : fixed (ansatz) parameters
% h : integration step

L = length(X);

% l = 0
J(1,L) = -h*(X(L-1) - X(2));
J(1,L-1) =- h* X(L);
J(1,1) = 1 - h*( 1 + A(1) + 2*A(2)*X(1) );
J(1,2) = h*X(L);

% l = 1
J(2,1) = -h*( X(L) - X(3) );
J(2,L) = -h*X(1);
J(2,2) = 1 - h*( 1+A(1)+2*A(2)*X(2) );
J(2,3) = h*X(1);

% l = L-1
J(L,L-1) = -h*( X(L-2) - X(1) );
J(L,L-2) = -h*X(L-1);
J(L,L) = 1 - h*( 1 + A(1) + 2*A(2)*X(L) );
J(L,1) = h*X(L-1);

% l = 2, ..., L-2
for idl = 3:L-1
    J(idl,idl-1) = - h.*( X(idl-2)-X(idl+1) );
    J(idl,idl-2) = - h.*X(idl-1);
    J(idl,idl) = 1 - h.*(1+A(1)+2.*A(2).*X(idl));
    J(idl,idl+1) = h*X(idl-1);
end; %for


end

