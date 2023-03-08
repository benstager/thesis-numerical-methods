
% Euler: usual formula
% Heun: usual formula
% Backwards: Nt0 = 10^4, Nx0 = 31
% Imp Midpoint: Nt0 = 10^2, Nx0 = 15

% # of interior spatial points
Ns = 2.^(2:6) - 1;
Ns_ref = 2^7-1;
Nt_ref = 2*10^4;
% amount of temporal points based on spatial grid
methods = {@eulerLin,@heunLin,@backwardsEulerLin,@impMidpointLin};
% My old parameters
% dts = {@(n) 2/(n+1)^2, @(n) (2^2(2/n)^4)/(8*(2/n)^2+6)}
% epsilon = .7;
% tspan = [0,.01];
% IMEX: referenceBurgers(epsilon,10^4,127,tspan), Nt0 = 100, Nx0 = 63
% ExpInt: referenceBurgers(epsilon,10^4,127,tspan), Nt0 = 100, Nx0 = 63
% the same???
% euler = @(n) 2/((n+1)^2);
% heun = @(n) (4*(2/n)^4)/((8*(2/n)^2+6));
% backward = @(n) diff(tspan)/(10^4/(15/n)^2);
% imp = @(n) diff(tspan)/(10^2/(31/n));
tspan = [0,1];
dts = {@(n) 2/((n+1)^2), @(n) 2/((n+1)^2), @(n) diff(tspan)/(10^4/(15/(n+1))^2), @(n) diff(tspan)/(10^2/(31/(n+1)))};
tspan = [0,1];
xspan = [-1,1];
xs = linspace(xspan(1),xspan(2),Ns_ref+2)';
ts = linspace(tspan(1),tspan(2),Nt_ref+1);
dxs = diff(xspan)./Ns;
exact = referenceHeat(Nt_ref,Ns_ref,tspan);
% exact = pdeSolverNEW(xs(2:Ns_ref+1),ts);
error = zeros(length(Ns),length(methods));
time = zeros(length(Ns),length(methods));


for i = 1:length(methods)
    for j = 1:length(Ns)
        Nt = ceil(diff(tspan)/dts{i}(Ns(j)));
        y0 = heatIC(Ns(j));
        [A,f] = heatOperators(Ns(j));
        [ys,cpu_time] = methods{i}(A,0,tspan,y0,Nt,'heat',1);
        error(j,i) = spatialError(ys(:,end),exact(:,end));
        time(j,i) = cpu_time;
    end
end

methods = {@euler,@heun,@backward,@imp};
getMethodName = @(f) functions(f).function;
legend_entries = cellfun(getMethodName, methods, 'UniformOutput', false);
% 
% figure(1)
% loglog(Ns,error,LineWidth=2.0);
% xlabel('grid size(Nx)'); ylabel('error'); title('Convergence Diagram for Heat Equation');
% legend(legend_entries);
% ylim([10^-12,10^2])

figure(2)
loglog(time,error,LineWidth=2.0);
xlabel('time (sec)'); ylabel('error'); title('Precision Diagram for Heat Equation');
legend(legend_entries);
ylim([10^-12,10^2])