function [ys,cpu_time] = expIntegrator(L,Nl,tspan,y0,N)

% Imex scheme for Burgers
% f = Ly + N(y);
% How to scheme generally?

ys = zeros(length(y0),N+1);
ys(:,1) = y0;
y = y0;
dt = diff(tspan)/N;
I = speye(length(y0));
% Burgers" 
p0 = expm(dt*L);
p1 = L\(I-expm(dt*L));
tic
for i = 1:N
    y = p0*y - p1*Nl(y);
    ys(:,i+1) = y;
end
cpu_time = toc;