function [UN,x,t] = oneDTransport()

% Inital data
nx = 100;
nt = 150;
a = 0;
b = 5;
t0 = 0;
t1 = 1;

% Time and spatial steps and wave number
dx = (b-a)/(nx - 1);
dt = (t1 - t0)/(nt - 1); 
c = .01;

% Each x and t vectors to evaluate u(x,t) at
x = a:dx:b;
t = t0:dt:t1;

% Coeffecient of finite approximation
alpha = -c*dt/dx;


% NUMERICAL SOLUTION


% Grid
UN = zeros(nx,nt);


% Initial conditions

UN(:,1) = sin(pi*x);


% Boundary conditions

UN(a+1,:) = 0;

% Assigning values to each u(xi,ti)
for j = 1:nt-1
    for i = 2:nx-1
        UN(i,j+1) = alpha*(UN(i+1,j)-UN(i,j)) + UN(i,j);
    end  
end

figure()
contourf(UN,200,'linecolor','non')
title('1D Heat eq')
xlabel('t')
ylabel('x')
colormap(jet(256))
colorbar
surf(t,x,UN);
shading interp
