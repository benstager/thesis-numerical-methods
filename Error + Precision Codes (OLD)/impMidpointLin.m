function [ys,cpu_time] = impMidpointLin(A,B,tspan,y0,N,pde,epsilon)
% A matrix or scalar
% f = f(t,y)
% tspan = [start,end]
% Nt is number times
switch pde
    case {'dalquist','heat'}
        ys = zeros(length(y0),N+1);
        y = y0;
        ys(:,1) = y0;
        dt = diff(tspan)/N;
        I = speye(length(y0));
        tic
        for i = 1:N
            K1 = (I-dt*A/2)\(A*y);
            %     y = y + dt*((A*y)/(1-A*(dt/2)));
            y = y + dt*K1;
            ys(:,i+1) = y;
        end
        cpu_time = toc;
    case 'burgers'
        %stub
end

