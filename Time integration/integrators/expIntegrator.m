function [ys,cpu_time] = expIntegrator(L,Nl,tspan,y0,N)


ys = zeros(length(y0),N+1);
ys(:,1) = y0;
y = y0;
dt = diff(tspan)/N;

if(issparse(L))
    I = speye(length(y0));
else
    I = eye(length(y0));
end


p0 = expm(dt*L);
p1 = L\(I-expm(dt*L));
tic
for i = 1:N
    y = p0*y - p1*Nl(y);
    ys(:,i+1) = y;
end
cpu_time = toc;