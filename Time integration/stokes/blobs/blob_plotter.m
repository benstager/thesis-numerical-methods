

nx = 20;
ny = 20;
n = 3;
epsilon = .5;

[xks,fks,xs,ys] = stokes_parameters(nx,ny);
X = init_blob(n);

for i = 1:n
    z = blob_meshes(X(:,i),epsilon,xs,ys);
    surf(z); shading interp; colormap jet;
    hold on;
end
