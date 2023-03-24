
% Plots each blob at n amount of (x0,y0) points

nx = 40;
ny = 40;
n = 6;
epsilon = .5;

[xks,fks,xs,ys] = stokes_parameters(nx,ny);

X_unshaped = init_blob(n);
X = reshape(X_unshaped,[2*n,1]);


for i = 1:n
    z = blob_meshes(X(:,i),epsilon,xs,ys);
    surf(xs,ys,z); shading interp; colormap jet;
    hold on;
end

% fks = @(t) [2*sin(t)+1;0];
tspan = [0,.5];
f = @(t,X)  velocity(X,[0;0],fks(t))';
Nt = 1000;

for i = 1
    [Xs,cpu] = euler(f,tspan,X(:,i)',Nt);
end



