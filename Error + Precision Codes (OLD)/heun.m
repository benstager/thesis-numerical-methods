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
    y_hat = y + dt*f(t,y);
    y = y + dt/2*(f(t,y)+ f(t,y_hat));
    ys(:,i+1) = y;
    t = t+ dt;
end
cpu_time = toc;


