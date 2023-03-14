function [A, f] = heatOperators(Nx, sparse)
% returns right hand side for the heat equation
% Nx - number of internal grid points

if(nargin == 1)
    sparse = true;
end

xspan = [-1,1];
dx    = diff(xspan) / (Nx+1);
A     = Uxx(Nx,dx, sparse);
f     = @(t,U) A * U;

end