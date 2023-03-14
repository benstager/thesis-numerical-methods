function [ys,cpu_time] = heunLin(A,tspan,y0,N)
% Huen method for solving the autonomous linear differential equation
%
%       y' = A y        y(tspan(1)) = y0
%
% Parameters
%   A       - matrix that defines linear system
%   tspan   - [start,end]
%   y0      - initial condition
%   N       - number of timesteps

ys = zeros(length(y0),N+1);
ys(:,1) = y0;

y = y0;
dt = diff(tspan) / N;

tic
for i = 1:N
    y = y + dt * A * y + (dt^2) * ( A * ( A * ( y / 2 )));
    ys(:,i+1) = y;
end
cpu_time = toc;

end