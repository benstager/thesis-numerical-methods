

x = referenceBurgers(.02,10^4,63,[0,.4]);
xs = linspace(-1,1,63);
ts = linspace(0,.4,10^4+1);
% surf(ts,xs,x); shading interp; colormap jet;
% xlabel('t'); ylabel('x'); zlabel('u(x,t)');
% title('Numerical solution of Viscous Burgers Equation in 1-D');

for i = 1:10
    plot(x(:,i*1000));
    hold on;
 
end
