
% Surfs the vx and vy components of velocity, can not be done
% simultaneously

% Mesh parameters
nx = 200;
ny = 200;

vx = zeros(ny,nx);
vy = zeros(ny,nx);
[xks,fks,xs,ys] = stokes_parameters(nx,ny);

% Vx surface
for i = 1:nx
    for j = 1:ny
        v = velocity_regularized([xs(i),ys(j)],xks,fks);
        vx(j,i) = v(1,1);
    end
end

% Vy surface
for i = 1:nx
    for j = 1:ny
        v = velocity_regularized([xs(i),ys(j)],xks,fks);
        vy(j,i) = v(2,1);
    end
end

figure(1)
surf(xs,ys,vx);
title('Vx');
shading interp
colormap jet

figure(2)
surf(xs,ys,vy);
title('Vy');
shading interp
colormap jet
