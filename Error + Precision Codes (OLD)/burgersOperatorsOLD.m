function [f,L,Nl] = burgersOperatorsOLD(epsilon,Nx)

xspan = [-1,1];
dx = diff(xspan)/(Nx+1);
f = @(t,U) epsilon*Uxx(Nx,dx)*U - Ux(Nx,dx)*(U.^2/2);
L =  epsilon*Uxx(Nx,dx);
Nl = @(U) - Ux(Nx,dx)*(U.^2/2);

end

% referenceBurgers(.02,10^4,128,[0,.4])