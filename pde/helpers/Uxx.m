function [A] = Uxx(Nx, dx, sparse)
% returns second-order finite difference matrix for the second derivative
% assuming Dirichlet boundary conditions
% -- Parameters
% Nx - number of internal grid points
% dx - grid spacing
% sparse - if true, returns sparse matrix, if false returns dense matrix

if(sparse) % create sparse operators
    o = ones(Nx, 1) / dx^2;
    A = spdiags([o, -2*o, o], [-1 0 1], Nx, Nx);
else % create non-sparse operators
    A = zeros(Nx,Nx);
    A(1:1+Nx:Nx*Nx) = -2/dx^2;
    A(Nx+1:1+Nx:Nx*Nx) = 1/dx^2;
    A(2:1+Nx:Nx*Nx-Nx) = 1/dx^2;
end

end