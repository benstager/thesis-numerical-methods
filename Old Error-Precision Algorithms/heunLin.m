function [ys,cpu_time] = heunLin(A,B,tspan,y0,N,pde,epsilon)
% A matrix or scalar
% f = f(t,y)
% tspan = [start,end]
% Nt is number times
ys = zeros(length(y0),N);
y = y0;
dt = diff(tspan)/N;
switch pde
    case {'dalquist','heat'}
        tic
        for i = 1:N
            y = y+dt*A*y + (dt^2)*(A*(A*(.5)*y));
            ys(:,i+1) = y;
        end
        cpu_time = toc;
    case 'burgers'
        tic
        for i = 1:N
            k1 = epsilon*A*y + B*y.^2/2;
            k2 = epsilon*A*(y + dt*k1) + B*(y + dt*k1).^2/2;
            y = y + dt*(.5*k1 + .5*k2);
            ys(:,i+1) = y;
        end
        cpu_time = toc;

end

