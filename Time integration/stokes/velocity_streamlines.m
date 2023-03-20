

% Function should produce streamline curves in the parametrized x-y plane

% Setting initial values
ny = 5;
nx = 5;
[xks,fks,xs,ys] = stokes_parameters(ny,nx);


% Producing a velocity vector at certain X = (x;y) position

y0 = [.3;.5];
% F = @(t,X) velocity(X,xks,fks);

for i = 1:nx
    for j = 1:ny
        v = velocity([xs(i),ys(j)],xks,fks);
        [ysx,cpux] = eulerLin(v(1),[xs(1),xs(end)],y0(1),nx);
        [ysy,cpuy] = eulerLin(v(2),[ys(1),ys(end)],y0(2),ny);
        plot(ysx,ysy);
        hold on;
    end
end

% Producing garbage currently
