function [xks,fks,xs,ys] = stokes_parameters(nx,ny)

% Empty argument function that returns forces and their respective
% locations, as well as x and y domains


% Forces and locations
xks = [[0;0],[.2;.4],[.1;.7]];
fks = [[1;1],[-2;-2],[.5;.3]];

% xks = [0;0];
% fks = [1;1];

% Domain and meshes
x = [-1,1];
y = [-1,1];
xs = linspace(x(1), x(2), nx);
ys = linspace(y(1), y(2), ny);



end

