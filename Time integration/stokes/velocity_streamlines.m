

% Function should produce streamline curves in the parametrized x-y plane

% Setting initial values
ny = 20;
nx = 20;
[xks,fks,xs,ys] = stokes_parameters(nx,ny);


% Producing a velocity vector at certain X = (x;y) position
tspan = [0,1];
f = @(t,X) velocity_regularized(X,xks,fks)';
Nt = 1000;

for i = 1:nx
    for j = 1:ny
        [X,cpu] = euler(f,tspan,[xs(i),ys(j)],Nt);
        plot(X(1,:),X(2,:),color = 'blue');
        hold on;
    end
end

% Producing garbage currently
