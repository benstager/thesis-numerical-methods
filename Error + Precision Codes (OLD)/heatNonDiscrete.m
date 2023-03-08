function [U0] = heatNonDiscrete()
clear; clc; %close all;

% Attempt to solve the non-discretized heat diffusion problem

% Initial data
xleft = -1; % left length
xright = 1; % right length
Nx = 32; % number of interior grid points
tfinal = .01;
Nt = 10^5;
dt = tfinal/Nt; % time step
dx = (xright-xleft)/(Nx+1); % why NX+1?
f = @(x) sin(pi*x); % initial condition
method = 3;
epsilon = .7; % for inviscid Burger's

% Creating column matrices that are each time and space step
xs = linspace(xleft,xright,Nx+2);   
ts = linspace(0,tfinal,Nt+1);

% finite difference matrix
A = Uxx(Nx,dx);
B = Ux(Nx,dx);
% Initial condition at interior points?


X0 = f(xs(2:Nx+1))';

% Solution 
U0 = zeros(Nx+2,Nt+1);
% Only accessed at interior points
U0(2:Nx+1,1) = X0;


% This is A, our scalar matrix


% Each column of U, each Ui,j at each position i 
switch(method)
    % Forward Euler
    case 0
        for i = 1:Nt
            X0 = X0 + dt*(epsilon*A*X0 + B*X0.^2/2);
            U0(2:Nx+1,i+1) = X0;
        end
    % Backward Euler
    case 1
       I = eye(Nx);
       B = (I-dt*A);
       for i = 1:Nt
           X0 = B\X0;
           U0(2:Nx+1,i+1) = X0;
       end
    % RK4
    case 2
       for i = 1:Nt
           K1 = dt*A*X0;
           K2 = dt*A*(X0 + K1/2);
           K3 = dt*A*(X0 + K2/2);
           K4 = dt*A*(X0 + K3);
           X0 = X0 + K1/6 + K2/3 + K3/3 + K4/6;
           U0(2:Nx+1,i+1) = X0;
       end
    % Heun
    case 3
        for i = 1:Nt
%             X0 = X0 + dt*(.5*(epsilon*A*X0+B*X0.^2/2) + ...
%                 .5*(epsilon*A*(X0+dt*epsilon*A+B*X0.^2/2)+ ...
%                 .5*B*(X0+dt*epsilon*A*X0+B*X0.^2/2).^2/2));
            % U ' = epsilon*A*U + BU^2/2
            k1 = epsilon*A*X0 + B*X0.^2/2;
            k2 = epsilon*A*(X0 + dt*k1) + B*(X0 + dt*k1).^2/2;
            X0 = X0 + dt*(.5*k1 + .5*k2);
            U0(2:Nx+1,i+1) = X0;
        end
    % Explicit Midpoint (bad stability)
    case 4
        for i = 1:Nt
            X0 = X0 + dt*A*X0 + ((dt^2*A^2)/2)*X0;
            U0(2:Nx+1,i+1) = X0;
        end
     % Implicit Midpoint (great stability region)
    case 5
        for i = 1:Nt
            I = eye(Nx);
            B = (I-dt*.5*A);
            X0 = B\((I+.5*dt*A)*X0);
            U0(2:Nx+1,i+1) = X0;
        end
     % Crank-Nicholson (A-stability)
    case 6 
        for i = 1:Nt
            I = eye(Nx);
            K1 = A*X0;
            K2 = ((A+(A^2*dt)/2)/(I-A*dt/2))*X0;
            X0 = X0 + dt*(.5*K1+.5*K2);
            U0(2:Nx+1,i+1) = X0;
        end
     % DIRK method
    case 7
       for i = 1:Nt
           I = eye(Nx);
           K1 = (I-dt/4)\(A*X0);
           K2 = (A*X0+A*(X0*.5)*K1)/(I-(dt/4)*A);
           X0 = X0 + dt*(.5*K1+.5*K2);
       end

end

% Second way of solving
% Find Eigenvalues then evaluate each time function at each
% positional step
% set(gca,'YDir','reverse')
% view([-51,22])
end