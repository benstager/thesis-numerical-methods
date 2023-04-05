function [ys,cpu_time] = IMEXSemiLinearMidpointGMRES(F,tspan,y0,N)

% partioning y' = F(y) as y' = J(y)y + F(y) - J(y) y

ys = zeros(length(y0),N+1);
ys(:,1) = y0;
y = y0;
dt = diff(tspan)/N;
tol = 10^-2;
h = 10^-5;
t = tspan(1);

I = sparse(length(y0));

tic
for i = 1:N
    k1_hat = F(t,y) - jac_ProdFD(F,y,y,h)';
    k1 = jac_FD(F,y,h)*((I-dt/2*(jac_FD(F,y,h)))\(y+dt/2*k1_hat));
    k2_hat = F(t,y+dt/2*(k1+k1_hat)) + ...
    jac_ProdFD(F,y+dt/2*(k1+k1_hat),y+dt/2*(k1+k1_hat),h)';
    y = y+dt*(k1+k2_hat);
    ys(:,i+1) = y;
    t = t+h;
end
cpu_time = toc;
