function [ys,cpu_time] = IMEXSemiLinearEulerGMRES(F,tspan,y0,N)

% Integrator that transforms F into Ly + N(y) and uses Jacobian products
% to trick integrator into thinking it is of form y' = Ly + N(y)

ys = zeros(length(y0),N+1);
ys(:,1) = y0;
y = y0;
dt = diff(tspan)/N;
t = tspan(1);
I = eye(length(y0));
tol = 10^-2;
h = 10^-8;

tic
for i = 1:N
    J_yn = jac_FD(F,y,h);
    Jv = jac_ProdFD(F,y,y,h)';
    A = I - J_yn;
    [y,flag] = gmres(@(y) (I - J_yn)*y, y + dt*(F(t,y) - Jv),...
        length(A),tol);
    ys(:,i+1) = y;
    t = t + dt;
end
cpu_time = toc;

