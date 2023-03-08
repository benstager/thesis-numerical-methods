
close all;
clear all; %#ok<CLALL> 

eq = 'heat'; % heat dalquist type PDE ODE
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
        % exact = referenceHeat(10^6,Nx,tspan);
        spaceError = norm(exact(:,end)-expm(tspan(end)*A)*y0,'inf');
        errorNorm = @(approx) norm(approx-expm(tspan(end)*A)*y0,'inf');
        % errorNorm1 = @(approx) spatialError(approx,exact(:,end));
        methods = {@eulerLin,@heunLin,@backwardsEulerLin,@impMidpointLin};
    case 'burgers'
        Nx = 100;
        Nt = 10.^(1:4);
        tspan = [0,.01];
        xspan = [-1,1];
        dxs = diff(xspan)/(Nx+1);
        dts = diff(tspan)./Nt;
        f = @(x) sin(pi*x); %(x+1)*(x-1)
        xs = linspace(xspan(1),xspan(2),Nx+2);
        ts = linspace(tspan(1),tspan(2),50);
        y0 = f(xs(2:Nx+1))';
        epsilon = .7;
        A = Uxx(Nx,dxs);
        % Not working
        B = Ux(Nx,dxs);
        f = @(t,U) epsilon*A*U + B*(U.^2)/2;
        % HELP
        exact = heatNonDiscrete();
        errorNorm = @(approx) norm(approx-exact(2:Nx+1,end));
        methods = {@eulerLin,@heunLin,@backwardsEulerLin,@impMidpointLin};
    case 'dalquist'
        tspan = [0,1];
        y0 = 1;
        Nx = 'N/A';
        % Change for matrix
        A = 13;
        B = 1;
        epsilon = 1;
        f = @(t,y) A*y;
        Nt = 10.^(1:5);
        dts = diff(tspan)./Nt;
        errorNorm = @(approx) abs(approx-exp(tspan(end)*A));
        methods = {@eulerLin,@heunLin,@backwardsEulerLin,@impMidpointLin};
       
end

error = zeros(length(Nt),length(methods));
time = zeros(length(Nt),length(methods));

for i = 1:length(methods)
    for j = 1:length(Nt)
        [ys,cpu_time] = methods{i}(A,0,tspan,y0,Nt(j),'heat',1);
        error(j,i) = errorNorm(ys(:,end));
        time(j,i) = cpu_time;
    end
end
% loglog(dts,error,LineWidth = 2.0);
% hold on;
% loglog(dts,dts.^2,LineWidth = 2.0);
% xlabel('dt')
% ylabel('error')
% title("Error Convergence of Heun's Method (O(h^2))")
% legend('1');
% legend('2');
figure(1)
loglog(dts,error,LineWidth=2.0);
xlabel('time (sec)'); ylabel('error'); title('Error Convergence (Nt = 10^n, Nx = '+string(Nx) + ')');
%legend(legend_entries);


  figure(2)
  loglog(time,error,LineWidth=2.0);
  xlabel('time (sec)'); ylabel('error'); title('Precision Diagram (Nt = 10^n, Nx = ' + string(Nx) + ')');
  yline(spaceError, 'green', LineWidth = 2.0);
  %legend(legend_entries);
  ylim([1e-16,1e2])