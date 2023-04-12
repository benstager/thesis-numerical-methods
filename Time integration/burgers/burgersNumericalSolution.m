function [U,cpu_time] = burgersNumericalSolution(epsilon,Nt,Nx,tspan)

% computes Burgers solution with epsilon, and Nt, Nx mesh, on tspan

y0 = burgersParameters(Nx);
[f,L,N] = burgersOperators(Nx, epsilon);
U = IMEXdirk(L,N,tspan,y0,Nt);
% [U,cpu_time] = IMEXSemiLinearEulerBackslash(f,tspan,y0,Nt);

end