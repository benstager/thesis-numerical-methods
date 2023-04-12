
% Plots movement of blob centers for certain force
close all; clear all;

nx = 40;
ny = 40;
n = 10;
epsilon = .5/10;

[xks,fks,xs,ys] = stokes_parameters(nx,ny);
X = init_blob(n);

Xv = reshape(X,[2*n,1]);
tspan = [0,1];
Nt = 250;
[Xvs,cpu] = IMEXSemiLinearEulerGMRES(@f,tspan,Xv,Nt);

figure()
for i = 1:Nt
    Xvi = Xvs(:,i);
    Xvi = reshape(Xvi,[2,length(Xvi)/2]);
    plot(Xvi(1,:),Xvi(2,:),'ko-');
    xlim([-1.5,1.5]); ylim([-1.5,1.5]);
    drawnow()
end
