function [U] = burgersNumericalSolution(epsilon,Nt,Nx,tspan)

% computes Burgers solution with epsilon, and Nt, Nx mesh, on tspan

y0 = burgersParameters(Nx);
[f,L,N] = burgersOperators(Nx, epsilon);
U = expRK(L,N,tspan,y0,Nt);

end