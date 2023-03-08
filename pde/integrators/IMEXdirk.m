function [ys,cpu_time] = IMEXdirk(L,Nl,tspan,y0,N)

ys = zeros(length(y0),N+1);
ys(:,1) = y0;
y = y0;
dt = diff(tspan)/N;

if(issparse(L))
    I = speye(length(y0));
else
    I = eye(length(y0));
end
delta = (-2*sqrt(2))/3;
gamma = (2-sqrt(2))/2;

tic
for i = 1:N
    k1_hat = Nl(y);
    k1 = (I-dt*gamma*L)\(L*(y+dt*gamma*k1_hat));
    u1 = y + dt*gamma*k1 + dt*gamma*k1_hat;
    k2_hat = Nl(u1);
    k2 = (I-dt*gamma*L)\(L*(y+dt*((1-gamma)*k1)+dt*((delta)*Nl(y)+delta*k2_hat)));
    u2 = y+dt*((1-gamma)*k1+gamma*k2)+dt*(delta*k1_hat+(1-delta)*k2_hat);
    k3_hat = Nl(u2);
    y = y+dt*((1-gamma)*k1+gamma*k2) + dt*((1-gamma)*k1_hat+gamma*k2_hat);
    ys(:,i+1) = y;
end

cpu_time = toc;

