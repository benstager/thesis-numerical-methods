function [ys,cpu_time] = expRK(L,Nl,tspan,y0,N)

% exponential runge kutta for y' = Ly + N(y)

ys = zeros(length(y0),N+1);
ys(:,1) = y0;
y = y0;
dt = diff(tspan)/N;

if(issparse(L))
    I = speye(length(y0));
else
    I = eye(length(y0));
end


tic

p0 = expm(dt*L);
p1 = L\(expm(dt*L)-I);
p2 = (L\(L\((expm(dt*L)-I-dt*L))))/dt;

for i = 1:N
    a_n = p0*y + p1*Nl(y);
    y = a_n + p2*(Nl(a_n)-Nl(y));
    ys(:,i+1) = y;
end
cpu_time = toc;