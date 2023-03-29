
% Script plots the temporal error of centers for sequence of time steps
% taken to solve highly viscous stokes flow
clear all;

nx = 40;
ny = 40;
n = 6;
epsilon = .5;

[xks,fks,xs,ys] = stokes_parameters(nx,ny);
X = init_blob(n);

for i = 1:n
    z = blob_meshes(X(:,i),epsilon,xs,ys);
    % surf(xs,ys,z); shading interp; colormap jet;
    % hold on;
end

% reference fine solution
Xv = reshape(X,[2*n,1]);
tspan = [0,1];
Nt_ref = 2^5;
[Xvs,cpu] = euler(@f,tspan,Xv,Nt_ref);
Xvi = Xvs(:,end);
Xvi = reshape(Xvi,[2,length(Xvi)/2]);
x_ref = Xvi(1,:);
y_ref = Xvi(2,:);

% computing error for methods
methods = {@euler, @heun};
Nts = 2.^(1:4);
dts = diff(tspan)./Nts;

error = zeros(length(Nts),length(methods));
time = zeros(length(Nts),length(methods));
for i = 1:length(methods)
    for j = 1:length(Nts)
        [Xvs,cpu] = methods{i}(@f,tspan,Xv,Nts(j));
        Xvi = Xvs(:,end);
        Xvi = reshape(Xvi,[2,length(Xvi)/2]);
        x_coarse = Xvi(1,:);
        y_coarse = Xvi(2,:);
        error(j,i) = norm([norm(x_ref-x_coarse,2);norm(y_ref-y_coarse,2)],2);
        time(j,i) = cpu;
    end
end

getMethodName = @(f) functions(f).function;
legend_entries = cellfun(getMethodName, methods, 'UniformOutput', false);

% convergence diagram
figure(1)
loglog(dts,error,LineWidth=2.0); hold on;
xlabel('step size (Nt)');
ylabel('error');
title('Convergence Diagram for Stokes Flow');
legend(legend_entries); legend box on;


% precision diagram
figure(2)
loglog(time,error,LineWidth=2.0);
xlabel('time (sec)');
ylabel('error');
title('Precision Diagram for Stokes Flow');
legend(legend_entries); legend box off;

        
