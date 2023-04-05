function [ys,cpu_time] = IMEXSemiLinearEulerBackslash(F,tspan,y0,N)

% Integrator that transforms F into Ly + N(y) and uses Jacobian

ys = zeros(length(y0),N+1);
ys(:,1) = y0;
y = y0;
dt = diff(tspan)/N;
t = tspan(1);
I = eye(length(y0));

tic
for i = 1:N
    y = (I-dt*jac_FD(F,y,dt))\(y + dt*(F(t,y)-jac_ProdFD(F,y,y,dt)'));
    ys(:,i+1) = y;
    t = t + dt;
end
cpu_time = toc;

