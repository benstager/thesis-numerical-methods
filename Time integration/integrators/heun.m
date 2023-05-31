function [ys,cpu_time] = heun(f,tspan,y0,N)
% f = f(t,y)
% tspan = [start,end]
% Nt is number times
ys = zeros(length(y0),N+1);
ys(:,1) = y0;
y = y0;
dt = diff(tspan)/N;
t = tspan(1);

tic
for i = 1:N
    k1 = f(t,y);
    k2 = f(t+dt,y+dt*k1);
    y = y +.5*dt*(k1+k2);
    ys(:,i+1) = y;
    t = t+ dt;
end
cpu_time = toc;


