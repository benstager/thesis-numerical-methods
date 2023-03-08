clear all;

eq = 'burgers'; % heat dalquist type PDE ODE
switch eq
    case 'heat'
        Nx = 50;
        Nt = 10.^(1:5);
        tspan = [0,1];
        xspan = [-1,1];
        dxs = diff(xspan)/(Nx+1);
        dts = diff(tspan)./Nt;
        f = @(x) 8*pi/3*sin(pi*x); %(x+1)*(x-1)
        xs = linspace(xspan(1),xspan(2),Nx+2)';
        ts = linspace(tspan(1),tspan(2),50);
        y0 = f(xs(2:Nx+1));
        % Change for matrix2
        epsilon = 1;
        A = sparse(Uxx(Nx,dxs));
        B = Ux(Nx,dxs);
        f = @(t,U) A*U;
        % Two methods - 1. solve using fourier series, 2. solve with hight
        % Nt and Rk4 method
        exact = pdeSolverNEW(xs(2:Nx+1),tspan(end));
        spaceError = norm(exact-expm(tspan(end)*A)*y0,'inf');
        errorNorm = @(approx) norm(approx-expm(tspan(end)*A)*y0);
    case 'burgers'
        Nx = 31;
        Nt = 10.^(1:4);
        tspan = [0,.4];
        xspan = [-1,1];
        dxs = diff(xspan)/(Nx+1);
        dts = diff(tspan)./Nt;
        ts = linspace(tspan(1),tspan(2),50);
        [y0,xs] = burgersIC(Nx);
        epsilon = .02;
        [f,L,N] = burgersOperators(epsilon,Nx);
        % can be extended to autonomous problems
        % HELP
        exact = referenceBurgers(epsilon,10^5,31,tspan);
        errorNormSpace = @(approx) spatialError(approx,exact(:,end));
        errorNormTime = @(approx) norm((exact(:,end) - approx));
    case 'dalquist'
        tspan = [0,5];
        y0 = 1;
        Nx = 'N/A';
        % Change for matrix
        A = -1;
        B = 1;
        epsilon = 1;
        f = @(t,y) A*y;
        Nt = 2.^(1:5);
        dts = diff(tspan)./Nt;
        errorNorm = @(approx) abs(approx-exp(tspan(end)*A));
       
end

methodsIMEX = {@IMEXeuler, @expIntegrator};
methodsClassical = {@euler, @heun};
Nmi = length(methodsIMEX);
Nm = length(methodsClassical);
time = zeros(length(Nt),Nm + Nmi);
error = zeros(length(Nt),Nm+Nmi);

for i = 1:Nm
    for j = 1:length(Nt)
        [ys,cpu_time] = methodsClassical{i}(f,tspan,y0,Nt(j));
        error(j,i) = errorNormTime(ys(:,end));
        time(j,i) = cpu_time;
    end
end

for i = Nm+1:Nmi+Nm
    for j = 1:length(Nt)
        [ys,cpu_time] = methodsIMEX{i-2}(L,N,tspan,y0,Nt(j));
        error(j,i) = errorNormTime(ys(:,end));
        time(j,i) = cpu_time;
    end
end

% for i = 1:Nm
%     for j = 1:length(Nt)
%         [ys,cpu_time] = methodsClassical{i}(f,tspan,y0,Nt(j));
%         error(j,i+Nmi) = errorNorm(ys(:,end));
%         time(j,i+Nmi) = cpu_time;
%     end
% end
% loglog(dts,error,LineWidth = 2.0);
% hold on;
% loglog(dts,dts.^2,LineWidth = 2.0);
% xlabel('dt')
% ylabel('error')
% title("Error Convergence of Heun's Method (O(h^2))")
% legend('1');
% legend('2');

methods = [methodsIMEX,methodsClassical];
getMethodName = @(f) functions(f).function;
legend_entries = cellfun(getMethodName, methods, 'UniformOutput', false);

figure(1)
loglog(Nt,error,LineWidth=2.0);
xlabel('stepsize (dt)'); ylabel('error'); title('Error Convergence (Nt = 10^n, Nx = '+string(Nx) + ')');
legend(legend_entries);
ylim([1e-16,1e2])


figure(2)
loglog(time,error,LineWidth=2.0);
xlabel('time (sec)'); ylabel('error'); title('Precision Diagram (Nt = 10^n, Nx = ' + string(Nx) + ')');
legend(legend_entries);
ylim([1e-16,1e2])