

% # of interior spatial points
Ns = 2.^(5:7) - 1;
Ns_ref = 2^9-1;
Nt_ref = 5000;
% amount of temporal points based on spatial grid
methods = {@euler,@heun};
% My old parameters
% dts = {@(n) 2/(n+1)^2, @(n) (2^2(2/n)^4)/(8*(2/n)^2+6)}
% epsilon = .7;
% tspan = [0,.01];
epsilon = .02;
% IMEX: referenceBurgers(epsilon,10^4,127,tspan), Nt0 = 100, Nx0 = 63
% ExpInt: referenceBurgers(epsilon,10^4,127,tspan), Nt0 = 100, Nx0 = 63
% the same???
dts = {@(n) 2/(epsilon*(n+1)^2), @(n) (4*(2/n)^4)/(epsilon*(8*(2/n)^2+6))};
tspan = [0,.4];
xspan = [-1,1];
dxs = diff(xspan)./Ns;
exact = referenceBurgers(epsilon,Nt_ref,Ns_ref,tspan);

error = zeros(length(Ns),length(methods));
time = zeros(length(Ns),length(methods));


for i = 1
    for j = 1:length(Ns)
        Nt = ceil(diff(tspan)/dts{i}(Ns(j)));
        y0 = burgersIC(Ns(j));
        f = burgersOperatorsOLD(epsilon,Ns(j));
        [ys,cpu_time] = methods{i}(f,tspan,y0,Nt);
        error(j,i) = spatialError(ys(:,end),exact(:,end));
        time(j,i) = cpu_time;
    end
end

methods = {@euler,@heun};
getMethodName = @(f) functions(f).function;
legend_entries = cellfun(getMethodName, methods, 'UniformOutput', false);

figure(1)
loglog(dxs,error,LineWidth=2.0);
xlabel('stepsize (dx)'); ylabel('error'); title('test');
legend(legend_entries);
ylim([10^-12,10^2])

figure(2)
loglog(time,error,LineWidth=2.0);
xlabel('time (sec)'); ylabel('error'); title('test');
legend(legend_entries);
ylim([10^-12,10^2])


