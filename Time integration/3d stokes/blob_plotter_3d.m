
% Plots movement of blob centers for certain force
close all; clear all;

nx = 40;
ny = 40;
n = 15;
epsilon = .5/10;

[xks,fks,xs,ys] = stokes_parameters(nx,ny);
X = init_blob_3d(n);

Xv = reshape(X,[3*n,1]);
tspan = [0,2];
Nt = 10^4;
[Xvs,cpu] = euler(@f_3d,tspan,Xv,Nt);

figure()
for i = 1:Nt
    Xvi = Xvs(:,i);
    Xvi = reshape(Xvi,[3,length(Xvi)/3]);
    plot3(Xvi(1,:),Xvi(2,:),Xvi(3,:),'ko-');
    xlim([-1.5,1.5]); ylim([-1.5,1.5]), zlim([-1.5,1.5]);
    drawnow()
end
