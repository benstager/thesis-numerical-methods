function [f, L, N] = burgersOperators(Nx,epsilon, sparse)
% BURGERS OPERATORS: Returns finite difference matrices U_xx, U_x , and L, N operators

if(nargin == 2)
    sparse = true;
end

xspan = [-1,1];
dx    = diff(xspan) / (Nx+1);
A     = Uxx(Nx,dx, sparse);
B     = Ux(Nx,dx, sparse);
L     = epsilon*A;
N     = @(U) - B*U.^2/2;

f = @(t,U) epsilon*A*U - B*U.^2/2;


end

