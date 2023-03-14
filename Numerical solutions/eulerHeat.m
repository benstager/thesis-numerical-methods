function [ys,cpu_time] = eulerHeat(A,tspan,y0,N)
% A matrix or scalar 
% f = f(t,y)
% tspan = [start,end]
% Nt is number times
ys = zeros(length(y0),N+1);
ys(:,1) = y0;
y = y0;
dt = diff(tspan)/N;
tic
for i = 1:N
    y = y + dt*(A*y);
    ys(:,i+1) = y;
end
cpu_time = toc;
end

