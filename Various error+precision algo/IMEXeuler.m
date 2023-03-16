function [ys,cpu_time] = IMEXeuler(L,Nl,tspan,y0,N)

% Imex scheme for Burgers
% f = Ly + N(y);
% How to scheme generally?

ys = zeros(length(y0),N+1);
ys(:,1) = y0;
y = y0;
dt = diff(tspan)/N;
I = speye(length(y0));
% Burgers" 
tic
for i = 1:N
    y = (I-dt*L)\(y + dt*Nl(y));
    ys(:,i+1) = y;
end
cpu_time = toc;



