function [y0,xs,tspan] = heatParameters(Nx)
% Returns problem parameters for the heat equation, namely the initial 
% condition, the spatial grid and the integration interval 
% Nx - number of internal grid points

tspan = [0 1];
xspan = [-1,1];
xs = linspace(xspan(1),xspan(2),Nx+2).';
f = @(x) (8*pi/3)*sin(pi*x);
y0 = f(xs(2:Nx+1));

end