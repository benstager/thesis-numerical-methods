
% Plots movement of blob centers for certain force

nx = 40;
ny = 40;
n = 30;
epsilon = .5/10;

[xks,fks,xs,ys] = stokes_parameters(nx,ny);
X = init_blob(n);

% for i = 1:n
% %     z = blob_meshes(X(:,i),epsilon,xs,ys);
% % %     surf(xs,ys,z); shading interp; colormap jet;
% % %     hold on;
% end

Xv = reshape(X,[2*n,1]);
tspan = [0,1];
Nt = 1000;
[Xvs,cpu] = heun(@f,tspan,Xv,Nt);

figure()
for i = 1:Nt
    Xvi = Xvs(:,i);
    Xvi = reshape(Xvi,[2,length(Xvi)/2]);
    plot(Xvi(1,:),Xvi(2,:),'ko-');
    xlim([-1.5,1.5]); ylim([-1.5,1.5]);
    drawnow()
end
