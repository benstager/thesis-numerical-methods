

function [UN,x,t] = oneDWave()

% Inital data
nx = 100;
nt = 100;
a = 0;
b = 5;
t0 = 0;
t1 = 1;

% Time and spatial steps
dx = (b-a)/(nx - 1);
dt = (t1 - t0)/(nt - 1); 

% Each x and t vectors to evaluate u(x,t) at
x = a:dx:b;
t = t0:dt:t1;

% Coeffecient of finite approximation
alpha = dt^2/dx^2;


% NUMERICAL SOLUTION


% Grid
UN = zeros(nx,nt);


% Initial conditions

UN(:,1) = sin(pi*x);
UN(:,2) = sin(pi*x)*(1+dt);

% Boundary conditions

UN(a+1,:) = 0;
UN(nx,:) = 0;

% Assigning values to each u(xi,ti)
for j = 2:nt-1
    for i = 2:nx-1
        UN(i,j+1) = alpha*(UN(i+1,j)-2*UN(i,j)+ UN(i-1,j)) + 2*UN(i,j) - UN(i,j-1);
    end  
end

figure()
contourf(UN,200,'linecolor','non')
title('1D Wave eq')
xlabel('t')
ylabel('x')
colormap(jet(256))
colorbar
surf(t,x,UN);
shading interp