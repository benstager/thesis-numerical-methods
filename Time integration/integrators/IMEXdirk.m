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


ImL = (I - dt * gamma * L);

tic
for i = 1:N
    
    % First Stage
    Y1 = y;
    k1_hat = Nl(y);    
    k1     = L * Y1;  % not needed - can be removed, just included for clarity  
    
    % Second Stage
    Y2     = ImL \ (y + dt * gamma * k1_hat);
    k2_hat = Nl(Y2);
    k2     = L * Y2;
    
    % Third Stage
    Y3     = ImL \ (y + dt * (delta * k1_hat + (1 - delta) * k2_hat + (1 - gamma) * k2));
    k3_hat = Nl(Y3);
    k3     = L * Y3;

    % Output
    y = y + dt * ((1 - gamma) * (k2_hat + k2) + gamma * (k3_hat + k3));
    ys(:,i+1) = y;

end

cpu_time = toc;
