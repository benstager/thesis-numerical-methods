function final = pdeSolverOLD(Nx)


[T,X] = meshgrid(0:.0001:1,-1:2/Nx:.99);


U = (8*pi/3)*sin(pi*X).*exp(-pi^2.*T);
final = U(:,end);
surf(T,X,U)
xlabel('t')
ylabel('x')
zlabel('u(x,t)')
colormap jet
shading interp
