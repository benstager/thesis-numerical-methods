% =========================================================================
% Burgers equation convergence and precision diagrams for space / time
% refinement.
% =========================================================================

clear;

% Test Spatial Grid Points 2^exps - 1
exps = 4:10;

% compute Burgers fine numerical solution
Nx_ref = 2^(exps(end)+2) - 1;
Nt_ref = 5000;
epsilon = .02;
[y0,xs_ref,tspan] = burgersParameters(Nx_ref);
exact_numerical = burgersNumericalSolution(epsilon,Nt_ref,Nx_ref,tspan);

% methods{i} contains cell with { @integrator, @StepsizeSelector }
methods_classical = {
    { @euler,            @(n) 2 * diff(tspan) / ( epsilon * (n + 1)^2 ) } % stability bound
    { @heun,             @(n) 2 * diff(tspan) / ( epsilon * (n + 1)^2 ) } % stability bound
    };

methods_IMEX = {
    { @IMEXeuler,   @(n) stepsizeSelector(tspan, n, 50, 100, 2, 1) }   % obtained using Nx = 50 and [Userx,Usery] = ginput(1)
    { @expIntegrator,      @(n) stepsizeSelector(tspan, n, 50, 100, 2, 1) }     % obtained using Nx = 50 and [Userx,Usery] = ginput(1)
    { @IMEXmidpoint,    @(n) stepsizeSelector(tspan, n, 25, 60, 2, 2)}
    { @IMEXdirk,    @(n) stepsizeSelector(tspan, n, 25, 40, 2, 2)}
    { @expRK,   @(n) stepsizeSelector(tspan, n, 25, 40, 2, 2)}
    };

Nc = length(methods_classical);
Ni = length(methods_IMEX);
%% run numerical experiments

Nx = 2.^exps - 1;
error = zeros(length(Nx),Nc+Ni);
time = zeros(length(Nx),Nc+Ni);

for i = 1:Nc
    for j = 1:length(Nx)

        integrator = methods_classical{i}{1};
        stepsizeSelector = methods_classical{i}{2};

        Nt = ceil(diff(tspan)/stepsizeSelector(Nx(j)));
        [y0] = burgersParameters(Nx(j));
        [f,L,N] = burgersOperators(Nx(j),epsilon);
        [ys,cpu_time] = integrator(f,tspan,y0,Nt);
        error(j,i) = spatialError(ys(:,end),exact_numerical(:,end));
        time(j,i) = cpu_time;

    end
end


for i = 1:Ni
    for j = 1:length(Nx)

        integrator = methods_IMEX{i}{1};
        stepsizeSelector = methods_IMEX{i}{2};

        Nt = ceil(diff(tspan)/stepsizeSelector(Nx(j)));
        [y0] = burgersParameters(Nx(j));
        [f,L,N] = burgersOperators(Nx(j),epsilon);
        [ys,cpu_time] = integrator(L,N,tspan,y0,Nt);
        error(j,i+Nc) = spatialError(ys(:,end),exact_numerical(:,end));
        time(j,i+Nc) = cpu_time;

    end
end

%% Generate Plots

% legend names for methods
methods = {@euler,@heun,@IMEXeuler,@expIntegrator,@IMEXmidpoint,@IMEXdirk,@expRK};
getMethodName = @(f) functions(f).function;
legend_entries = cellfun(getMethodName, methods, 'UniformOutput', false);
set(0,'defaultAxesFontSize',13)

% convergence diagram
figure(1)
loglog(Nx,error,'-*',LineWidth=4.0); hold on;
loglog(Nx, Nx.^(-2), 'k--'); hold off;
xlabel('grid size (Nx)');
ylabel('error');
title('Convergence Diagram for Burgers Equation');
legend(legend_entries); legend box off;
%ylim([10^-12,10^2]); xlim('tight');

% precision diagram
figure(2)
loglog(time,error,'-*,LineWidth=4.0);
xlabel('time (sec)');
ylabel('error');
title('Precision Diagram for Burgers Equation');
legend(legend_entries); legend box off;
%ylim([10^-16,10^2])