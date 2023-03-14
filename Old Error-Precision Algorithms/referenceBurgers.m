function [U] = referenceBurgers(epsilon,Nt,Nx,tspan)


y0 = burgersIC(Nx);
[f,L,Nl] = burgersOperatorsOLD(epsilon,Nx);
U = IMEXeuler(L,Nl,tspan,y0,Nt);

end