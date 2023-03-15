
% Plotting pressure as a surface u = p(x,y)
x = [-1,1];
y = [-1,1];
nx = 201;
ny = 201;
xs = linspace(x(1), x(2), nx);
ys = linspace(y(1), y(2), ny);

[xks,fks] = stokes_parameters();
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