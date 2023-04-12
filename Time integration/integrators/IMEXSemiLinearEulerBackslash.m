function [ys,cpu_time] = IMEXSemiLinearEulerBackslash(F,tspan,y0,N)

% Integrator that transforms F into Ly + N(y) and uses Jacobian

ys = zeros(length(y0),N+1);
ys(:,1) = y0;
y = y0;
dt = diff(tspan)/N;
t = tspan(1);
I = eye(length(y0));
h = 10^-7;

tic
for i = 1:N
    J_yn = jac_FD(F,y,h);
    y = (I-dt*J_yn)\(y + dt*(F(t,y)-J_yn*y));
    ys(:,i+1) = y;
    t = t + dt;
end
cpu_time = toc;

