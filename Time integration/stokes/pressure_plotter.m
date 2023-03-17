
% Pressure plotter plots p = f(x,y), where p is evaluated at x,y in mesh


% Plotting pressure as a surface u = p(x,y)
nx = 200;
ny = 200;

[xks,fks,xs,ys] = stokes_parameters(nx,ny);
p = zeros(ny,nx);

% F = @ (t,X) velocity_field(X,xks,fks) X' = F(X)
% surface vx and vy and streamlines

for i = 1:nx
    for j = 1:ny
        p(j,i) = pressure([xs(i),ys(j)],xks,fks);
    end
end

surf(xs, ys, p);
shading interp 
colormap jet