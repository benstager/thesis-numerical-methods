function [U] = referenceHeat(Nt,Nx,tspan)


y0 = heatIC(Nx);
[A,f] = heatOperators(Nx);
U = impMidpointLin(A,0,tspan,y0,Nt,'heat',1);

end
