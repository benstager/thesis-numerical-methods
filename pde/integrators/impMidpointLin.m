function [ys,cpu_time] = impMidpointLin(A,tspan,y0,N)
% Implicit midpoint method for solving the autonomous linear differential equation
%
%       y' = A y        y(tspan(1)) = y0
%
% Parameters
%   A       - matrix that defines linear system
%   tspan   - [start,end]
%   y0      - initial condition
%   N       - number of timesteps

ys = zeros(length(y0),N+1);
y = y0;
ys(:,1) = y0;
dt = diff(tspan)/N;

if(issparse(A))
    I = speye(length(y0));
else
    I = eye(length(y0));
end

tic
for i = 1:N
    Y1 = (I - (dt/2) * A) \ y;
    y  = y + dt * A * Y1;
    ys(:,i+1) = y;
end
cpu_time = toc;

end

