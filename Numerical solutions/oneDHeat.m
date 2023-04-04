
clear; clc; 
% Inital data
Nx = 200;
Nt = 20000;
xspan = [-1,1];
tspan = [0,.1];

% Time and spatial steps
dx = diff(xspan)/(Nx-1);
dt = diff(tspan)/(Nt-1); 

% Each x and t vectors to evaluate u(x,t) at
xs = linspace(xspan(1),xspan(end),Nx);
ts = linspace(xspan(1),xspan(end),Nt);

% Coeffecient of finite approximation
alpha = dt/dx^2;

% Grid
u = zeros(Nx,Nt);

% Initial and boundary conditions
u(:,1) = sin(pi*xs);
u(1,:) = 0;
u(Nx,:) = 0;

% Assigning values to each u(xi,ti)
for j = 1:Nt-1
    for i = 2:Nx-1
        u(i,j+1) = alpha*(u(i+1,j)-2*u(i,j)+ u(i-1,j))+ u(i,j);
    end  
end

figure()
contourf(u,200,'linecolor','non')
surf(ts,xs,u)
title('1D Heat eq')
xlabel('t')
ylabel('x')
colormap(jet(256))
colorbar
shading interp 
hold on

