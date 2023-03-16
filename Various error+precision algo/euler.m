function [ys,cpu_time] = euler(f,tspan,y0,N)
% A matrix or scalar
% f = f(t,y)
% tspan = [start,end]
% Nt is number times
ys = zeros(length(y0),N+1);
ys(:,1) = y0;
y = y0;
dt = diff(tspan)/N;
t = tspan(1);
% Heat" 
tic
for i = 1:N
    y = y + dt*f(t,y);
    ys(:,i+1) = y;
    t = t+dt;
end
cpu_time = toc;




