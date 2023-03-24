

% Function should produce streamline curves in the parametrized x-y plane

% Setting initial values
ny = 10;
nx = 10;
[xks,fks,xs,ys] = stokes_parameters(nx,ny);


% Producing a velocity vector at certain X = (x;y) position
tspan = [0,5];
f = @(t,X)  velocity(X,xks,fks)';
Nt = 1000;

for i = 1:nx
    for j = 1:ny
        [X,cpu] = euler(f,tspan,[xs(i),ys(j)],Nt);
        plot(X(1,:),X(2,:),color = 'blue');
        hold on;
    end
end

% Producing garbage currently
