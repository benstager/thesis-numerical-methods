function [ys,cpu_time] = IMEXSemiLinearEulerGMRES(F,tspan,y0,N)

% Integrator that transforms F into Ly + N(y) and uses Jacobian products
% to trick integrator into thinking it is of form y' = Ly + N(y)

ys = zeros(length(y0),N+1);
ys(:,1) = y0;
y = y0;
dt = diff(tspan)/N;
t = tspan(1);
I = speye(length(y0));
tol = 10^-2;

tic
for i = 1:N
    J_yn = sparse(jac_FD(F,y,dt));
    Jv = jac_ProdFD(F,y,y,dt)';
    [y,flag] = gmres(I-dt*J_yn,y + dt*(F(t,y)-Jv),...
        size(I,1),tol);
    ys(:,i+1) = y;
    t = t + dt;
end
cpu_time = toc;

