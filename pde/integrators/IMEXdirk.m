function [ys,cpu_time] = IMEXdirk(L,Nl,tspan,y0,N)

ys = zeros(length(y0),N+1);
ys(:,1) = y0;
y = y0;
dt = diff(tspan)/N;
I = speye(length(y0));
delta = 1;
%(-2*sqrt(2))/3;
gamma = 1;
%(2-sqrt(2))/2;

tic
for i = 1:N
    k1_hat = Nl(y);
    k1 = (1-dt*gamma*L)^-1*L*(y+dt*gamma*Nl(y));
    u1 = y + dt*gamma*k1 + dt*gamma*k1_hat;
    k2_hat = Nl(u1);
    k2 = (I-dt*gamma*L)^-1*L*(y+dt*((1-gamma)*k1)+dt*((delta)*Nl(y)+delta*k2_hat));
    u2 = y+dt*((1-gamma)*k1+gamma*k2)+dt*(delta*k1_hat+(1-delta)*k2_hat);
    k3_hat = Nl(u2);
    y = y+dt*((1-gamma)*k1+gamma*k2) + dt*((1-gamma)*k2_hat+gamma*k3_hat);
    ys(:,i+1) = y;
end

cpu_time = toc;

