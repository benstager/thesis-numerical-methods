function [ys,cpu_time] = expRK(L,Nl,tspan,y0,N)



ys = zeros(length(y0),N+1);
ys(:,1) = y0;
y = y0;
dt = diff(tspan)/N;
I = speye(length(y0)); 

tic
for i = 1:N
    a_n = expm(dt*L)*y + (expm(dt*L)-I)*L^-1*Nl(y);
    y = a_n + (expm(dt*L)-I-dt*L)*1/dt*mpower(L,-2)*(Nl(a_n)-Nl(y));
    ys(:,i+1) = y;
end
cpu_time = toc;