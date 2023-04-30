function [ys,cpu_time] = IMEXSemiLinearEulerGMRES(F,tspan,y0,N)

% Integrator that transforms F into Ly + N(y) and uses Jacobian products
% to trick integrator into thinking it is of form y' = Ly + N(y)

ys = zeros(length(y0),N+1);
ys(:,1) = y0;
y = y0;
dt = diff(tspan)/N;
t = tspan(1);
dim = length(y0);
tol = 10^-4;
h = 10^-8;

tic
for i = 1:N
    Jv = jac_ProdFD(F,y,y,h)';
    [y,flag] = gmres(@(X) X-dt*jac_ProdFD(F,y,X,h)', y + dt*(F(t,y) - Jv),...
        dim,tol,dim);
    w = warning('query','last');
    warning('off', w.identifier);
    ys(:,i+1) = y;
    t = t + dt;
end
cpu_time = toc;

