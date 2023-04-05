function [ys,cpu_time] = IMEXeuler(L,Nl,tspan,y0,N)

% Imex scheme for y' = Ly + N(y)

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
for i = 1:N
    y = (I-dt*L)\(y + dt*Nl(y));
    ys(:,i+1) = y;
end
cpu_time = toc;



