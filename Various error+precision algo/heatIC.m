function [y0,xs] = heatIC(Nx)
% Returns Burgers IC, for Nx is internal gridpoints

xspan = [-1,1];
xs = linspace(xspan(1),xspan(2),Nx+2);
f = @(x) (8*pi/3)*sin(pi*x);
y0 = f(xs(2:Nx+1))';
end