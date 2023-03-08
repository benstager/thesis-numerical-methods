function [UN,x,t] = oneDHeat()
clear; clc; 
% Inital data
nx = 200;
nt = 20000;
a = -1;
b = 1;
t0 = 0;
t1 = 1;

% Time and spatial steps
dx = (b-a)/(nx - 1);
dt = (t1 - t0)/(nt - 1); 

% Each x and t vectors to evaluate u(x,t) at
x = a:dx:b;
t = t0:dt:t1;

% Coeffecient of finite approximation
alpha = dt/dx^2;


% NUMERICAL SOLUTION


% Grid
UN = zeros(nx,nt);


% Initial conditions

UN(:,1) = exp(-pi*x.^2);


% Boundary conditions

UN(1,:) = 0;
UN(nx,:) = 0;

% Assigning values to each u(xi,ti)
for j = 1:nt-1
    for i = 2:nx-1
        UN(i,j+1) = alpha*(UN(i+1,j)-2*UN(i,j)+ UN(i-1,j))+ UN(i,j);
    end  
end

figure()
contourf(UN,200,'linecolor','non')
surf(t,x,UN)
title('1D Heat eq')
xlabel('t')
ylabel('x')
colormap(jet(256))
colorbar
shading interp 
hold on

