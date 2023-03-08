Nx = 1000;
xspan = [-1,1];
dxs = diff(xspan)/(Nx+1);
xs = linspace(xspan(1),xspan(2),Nx+2)';
y0 = f(xs(2:Nx+1));
A = Uxx(Nx,dxs);
As = sparse(A);
trials = 20000;

tic
for i = 1:trials
%     A*y0;
end
toc

tic
for i = 1:trials
    As*y0;
end
toc

disp('Solve')

tic
for i = 1:trials
%     A\y0;
end
toc

tic
for i = 1:trials
    As\y0;
end
toc
% % crude matrix multiplication
% for i = 1:n
%     for j = i-1:i+1
%         x(i) = x(i) + A(i,j)*y0(j);
