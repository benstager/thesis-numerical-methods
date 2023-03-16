
% Function produces a phase plane in (x(t),y(t), solving dynamical system

% How do we find the solution to [vx;vy] = [f(x,y); g(x,y)]
% Y' = AY
% Can we return a functional handle for each vx and vy and go that way?




[xks,fks] = stokes_parameters();
F = @(t,X) velocity(X,xks,fks);

