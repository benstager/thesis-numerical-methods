
Nx = 2.^(5:8);
%Nt = @(N) (Nt0)/(exp((q/p)*log(Nx0/N)));

tspan = [0,1];
xspan = [-1,1];
dxs = @(N) diff(xspan)/(N+1);
%dts = diff(tspan)./Nt;
f = @(x) 8*pi/3*sin(pi*x); %(x+1)*(x-1)
% xs = linspace(xspan(1),xspan(2),Nx+2)';
% ts = linspace(tspan(1),tspan(2),50);
% Change for matrix2
epsilon = 1;
% A = sparse(Uxx(Nx,dxs));
% B = Ux(Nx,dxs);
% f = @(t,U) A*U;
% Two methods - 1. solve using fourier series, 2. solve with hight
% Nt and Rk4 method
% exact = pdeSolverNEW(xs(2:Nx+1),tspan(end));
% spaceError = norm(exact-expm(tspan(end)*A)*y0,'inf');
errorNorm = @(approx) norm(approx-expm(tspan(end)*A)*y0);
methods = {@eulerLin,@heunLin,@backwardsEulerLin,@impMidpointLin};
p = [1,2,1,2];
q = 2;
Nx0 = 15;
Nt0 = 600;
%Nt = @(N) (Nt0)/(exp((q/p)*log(Nx0/N)));
error = zeros(length(Nx),length(methods));
time = zeros(length(Nx),length(methods));

for i = 1:length(methods)
    for j = 1:length(Nx)
        Nt = @(N) (Nt0)/(Nx0/N)^(q/p(i));
        xs = linspace(xspan(1),xspan(2),Nx(j)+2)';
        y0 = f(xs(2:Nx(j)+1));
        exact = pdeSolverNEW(xs(2:Nx(j)+1),tspan(end));
        A = sparse(Uxx(Nx(j),dxs(Nx(j))));
        B = Ux(Nx(j),dxs(Nx(j)));
        [ys,cpu_time] = methods{i}(A,0,tspan,y0,ceil(Nt(Nx(j))),'heat',1);
        errorNorm = @(approx) norm(approx-expm(tspan(end)*A)*y0,'inf');
        error(j,i) = errorNorm(ys(:,end));
        time(j,i) = cpu_time;
    end
end

getMethodName = @(f) functions(f).function;
legend_entries = cellfun(getMethodName, methods, 'UniformOutput', false);

figure(1)
loglog(Nx,error,LineWidth=2.0);
xlabel('grid size (Nx)'); ylabel('error'); title('test1');
legend(legend_entries);
ylim([1e-16,1e2])


figure(2)
loglog(time,error,LineWidth=2.0);
xlabel('time (sec)'); ylabel('error'); title('test2');
legend(legend_entries);
ylim([1e-16,1e2])
