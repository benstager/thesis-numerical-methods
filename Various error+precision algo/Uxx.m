function [A] = Uxx(Nx,dx)

% Creating a separate file for scalar matrix

A = zeros(Nx,Nx);
A(1:1+Nx:Nx*Nx) = -2/dx^2;
A(Nx+1:1+Nx:Nx*Nx) = 1/dx^2;
A(2:1+Nx:Nx*Nx-Nx) = 1/dx^2;




% Heat eq: U' = DxxU
% Burgers: U' = epsilon*Dxx*U - Dx*U^2/2