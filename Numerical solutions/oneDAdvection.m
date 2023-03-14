function [UN,x,t] = oneDAdvection()

% Inital data
nx = 10;
nt = 1000;
a = 0;
b = 1;
t0 = 0;
t1 = 2;
alpha = .00001;
beta = .1;

% Time and spatial steps
dx = (b-a)/(nx - 1);
dt = (t1 - t0)/(nt - 1); 

% Each x and t vectors to evaluate u(x,t) at
x = a:dx:b;
t = t0:dt:t1;

% Coeffecient of finite approximation
phi = (alpha*dt^2)/dx^2;
psi = (beta*dt^2)/dx;


% NUMERICAL SOLUTION


% Grid
UN = zeros(nx,nt);


% Initial conditions

UN(:,1) = sin(pi*x);
UN(:,2) = sin(pi*x)*(1+dt);

% Boundary conditions

UN(a+1,:) = 0;
UN(b,:) = 0;

% Assigning values to each u(xi,ti)
for j = 2:nt-1
    for i = 2:nx-1
        UN(i,j+1) = phi*(UN(i+1,j)-2*UN(i,j)+ UN(i-1,j)) + psi*(UN(i+1,j)-UN(i,j)) + 2*UN(i,j) - UN(i,j-1);
    end  
end

figure()
contourf(UN,200,'linecolor','non')
title('1D Ad-Diff eq')
xlabel('t')
ylabel('x')
colormap(jet(256))
colorbar
surf(t,x,UN);
