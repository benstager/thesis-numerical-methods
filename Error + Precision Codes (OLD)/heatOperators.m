function [A,f] = heatOperators(Nx)

xspan = [-1,1];
dx = diff(xspan)/(Nx+1);
A = Uxx(Nx,dx);
f = @(t,U) A*U;

end

% referenceBurgers(.02,10^4,128,[0,.4])