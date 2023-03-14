function [A] = Ux(Nx,dx, sparse)
% returns second-order finite difference matrix for the first derivative 
% assuming Dirichlet boundary conditions
% -- Parameters
% Nx - number of internal grid points
% dx - grid spacing
% sparse - if true, returns sparse matrix, if false returns dense matrix

if(sparse) % create sparse operators
    o = ones(Nx, 1) / (2 * dx);
    A = spdiags([-o, o], [-1 1], Nx, Nx);
else % create non-sparse operators
    A = zeros(Nx,Nx);
    A(Nx+1:1+Nx:Nx*Nx) = 1/(2*dx);
    A(2:1+Nx:Nx*Nx-Nx) = -1/(2*dx);
end

end