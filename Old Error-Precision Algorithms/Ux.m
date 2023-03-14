function [A] = Ux(Nx,dx)

% Creating a separate file for scalar matrix

A = zeros(Nx,Nx);
A(1:1+Nx:Nx*Nx) = 0;
A(Nx+1:1+Nx:Nx*Nx) = 1/(2*dx);
A(2:1+Nx:Nx*Nx-Nx) = -1/(2*dx);