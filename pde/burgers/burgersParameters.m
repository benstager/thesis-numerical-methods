function [y0,xs,tspan] = burgersParameters(Nx)
% BURGERS PARAMETERS: returns initial condition, xs vector, and tspan array

xspan = [-1,1];
tspan = [0,.4];
xs = linspace(xspan(1),xspan(2),Nx+2);
f = @(x) -sin(pi*x);
y0 = f(xs(2:Nx+1))';

end

