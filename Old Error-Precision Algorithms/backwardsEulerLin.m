function [ys,cpu_time] = backwardsEulerLin(A,B,tspan,y0,N,pde,epsilon)
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
            y = (I-dt*A)\y;
            ys(:,i+1) = y;
        end
        cpu_time = toc;
    case 'burgers'
        %stub
end

